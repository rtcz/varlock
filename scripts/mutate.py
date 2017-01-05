import pysam
import struct
import random
import time
import numpy as np

CIGAR_MATCH = 0  # M
CIGAR_INS = 1  # I
CIGAR_DEL = 2  # D
CIGAR_REF_SKIP = 3  # N
CIGAR_SOFT_CLIP = 4  # S
CIGAR_HARD_CLIP = 5  # H
CIGAR_PAD = 6  # P
CIGAR_EQUAL = 7  # E
CIGAR_DIFF = 8  # X
NM_TAG = 9

BASES = ["A", "T", "G", "C"]


def multi_random(rnd, p_dist):
    """
    Secure multinomial random.
    :param rnd: (secure) random generator
    :param p_dist: Array of probabilities. When all probabilities are zero, each outcome has equal probability.
    Each value represents probability of one outcome.
    :return: Always returns index of outcome in p_dist array.
    """
    if len(p_dist) == 1:
        # no other choice
        return 0
    
    if sum(p_dist) == 0:
        # zero probabilities, return random index
        return rnd.randint(0, len(p_dist) - 1)
    
    p_level = 0
    rnd_value = rnd.random()
    p_value = rnd_value * sum(p_dist)
    
    for i in range(len(p_dist)):
        p_level += p_dist[i]
        if p_value < p_level:
            return i
    
    raise ValueError(
        "sum of probability distribution %d must be greater then the probability value %d" % (sum(p_dist), p_value))


def strip_prefix(value, prefix):
    if value.startswith(prefix):
        return value[len(prefix):]
    else:
        return value


def strip_chr(value):
    return strip_prefix(value, 'chr')


def index2pos(index, fai_list):
    """
    :param index:
    :param fai_list:
    :return: reference name and 0-based position
    """
    ref_pos = index
    for ref in fai_list:
        new_ref_pos = ref_pos - ref["length"]
        if new_ref_pos < 0:
            return strip_chr(ref["name"]), ref_pos
        else:
            ref_pos = new_ref_pos
    raise ValueError("reference position for index %d not found" % index)


def pos2index(ref_name, ref_pos, fai_list):
    sum_ref_pos = 0
    for ref in fai_list:
        if strip_chr(ref["name"]) == strip_chr(ref_name):
            return int(ref_pos) + sum_ref_pos
        else:
            sum_ref_pos += ref["length"]
    raise ValueError("sequence name %s not found" % ref_name)


def parse_fai(fai_filename):
    data = []
    ref_keys = ["name", "length", "offset", "linebases"]
    with open(fai_filename, "r") as fai_file:
        for line in fai_file:
            raw_record = line.rstrip().split()
            
            ref_name = raw_record[0]
            ref_data = map(int, raw_record[1:])
            
            ref = [ref_name] + list(ref_data)
            record_map = dict(zip(ref_keys, ref))
            
            data.append(record_map)
    
    return data


def is_before_snv(alignment, snv, fai_list):
    """
    :return: True if alignment end is mapped before snv
    """
    # reference_end points to one past the last aligned residue
    alignment_end = pos2index(alignment.reference_name, alignment.reference_end - 1, fai_list)
    return alignment_end < snv['index']


def is_after_snv(alignment, snv, fai_list):
    """
    :return: True if alignment start is mapped after snv
    """
    alignment_start = pos2index(alignment.reference_name, alignment.reference_start, fai_list)
    return alignment_start > snv['index']


def create_mut_map(pileup_col, ref_ac, rnd=random.SystemRandom()):
    """
    Creates mutation map for base pileup column.
    pileup base -> mutated base
    :param pileup_col: list containing A, T, G, C, N letters
    :param ref_ac:
    :param rnd: random number generator
    :return: mutation map
    """
    # count allele count in pileup base column
    alt_ac = [0] * 4
    for base in pileup_col:
        if base != 'N':
            alt_ac[BASES.index(base)] += 1
    
    ref_bases = BASES[:]
    alt_bases = BASES[:]
    
    # add random value to distinguish tied values
    alt_ac = [ac + rnd.random() for ac in alt_ac]
    
    # unknown base is always mapped to itself
    mut_map = {'A': None, 'T': None, 'G': None, 'C': None, 'N': 'N'}
    while len(ref_bases) > 0:
        # draw base with reference probability
        ref_base_id = multi_random(rnd, ref_ac)
        
        # draw most abundant base from alt pileup column
        alt_base_id = np.argmax(alt_ac)
        mut_map[alt_bases[alt_base_id]] = ref_bases[ref_base_id]
        
        # delete processed bases
        
        del ref_bases[ref_base_id]
        del ref_ac[ref_base_id]
        
        del alt_bases[alt_base_id]
        del alt_ac[alt_base_id]
    
    return mut_map


def next_item(iterable):
    try:
        return next(iter(iterable))
    except StopIteration:
        return None


def read_vac_snv(vac_file, fai_list):
    byte_string = vac_file.read(12)
    if byte_string == "":
        return None
    data_list = list(struct.unpack('<IHHHH', byte_string))
    ref_name, ref_pos = index2pos(data_list[0], fai_list)
    return {'index': data_list[0], 'ref_name': ref_name, 'ref_pos': ref_pos, 'ac': data_list[1:]}


def ref_pos2seq_pos(alignment, ref_pos):
    """
    Retrieve alignment position of reference base position
    :param alignment: pysam.AlignedSegment
    :param ref_pos: reference position in alignment to find sequence position
    :return: sequence position at reference position
    """
    for seq_pos, current_ref_pos in alignment.get_aligned_pairs(matches_only=True, with_seq=False):
        # search for base in snv position
        if current_ref_pos == ref_pos:
            return seq_pos
    return None


def the_method(vac_filename, bam_filename, out_filename, rnd, verbose=False):
    """
    :param vac_filename: variant allele count filename
    :param bam_filename: input bam filename
    :param out_filename: output bam filename
    """
    with pysam.AlignmentFile(bam_filename, "rb") as sam_file, \
            open(vac_filename, "rb") as vac_file, \
            pysam.AlignmentFile(out_filename, "wb", template=sam_file) as out_file:
        
        fai_list = parse_fai(FAI_FILENAME)
        snv_alignments = []
        snv_alignment_map = {}
        alignment = next_item(sam_file)
        snv = read_vac_snv(vac_file, fai_list)
        counter = 1
        while True:
            if verbose and counter % 1000 == 0:
                print("%d alignments processed" % counter)
            
            if snv is None or alignment is None:
                if verbose:
                    if snv is None:
                        print("EOF VAC")
                    else:
                        print("EOF BAM")
                
                for snv_alignment in snv_alignments:
                    out_file.write(snv_alignment)
                    counter += 1
                while alignment is not None:
                    out_file.write(snv_alignment)
                    counter += 1
                    alignment = next_item(sam_file)
                break
            
            elif alignment.is_unmapped:
                out_file.write(alignment)
                counter += 1
            
            elif is_before_snv(alignment, snv, fai_list):
                # print("alignment is before snv")
                # time.sleep(0.1)
                out_file.write(alignment)
                counter += 1
                alignment = next_item(sam_file)
            
            elif is_after_snv(alignment, snv, fai_list):
                
                # print("alignment is after snv")
                # time.sleep(0.1)
                
                # get pileup of bases from alignments at snv mapping position
                pileup_col = []
                for snv_alignment in snv_alignments:
                    try:
                        snv_seq_pos = snv_alignment_map[snv_alignment]
                    except KeyError:
                        print("alignment id: %d" % id(snv_alignment))
                        exit(1)
                    
                    if snv_seq_pos is not None:
                        # alignment is mapped at snv position
                        # not sure if query_sequence or query_alignment_sequence
                        base = snv_alignment.query_alignment_sequence[snv_seq_pos]
                        pileup_col.append(base)
                        
                        # print("pileup base %s, snv_seq_pos %d, seq %s" % (base, snv_seq_pos, snv_alignment.query_alignment_sequence))
                        # time.sleep(0.1)
                        
                        # print(pileup_col)
                        # print(snv['ac'])
                
                # mutation mapping
                mut_map = create_mut_map(pileup_col=pileup_col, ref_ac=snv['ac'], rnd=rnd)
                # print("mut_map %s" % str(mut_map))
                
                # mutate alignments
                for snv_alignment in snv_alignments:
                    # snv_seq_pos is sequence snv position
                    snv_seq_pos = snv_alignment_map[snv_alignment]
                    # snv_seq_pos must erased after use
                    snv_alignment_map[snv_alignment] = None
                    if snv_seq_pos is not None:
                        # alignment is mapped at snv position
                        # not sure if query_sequence or query_alignment_sequence
                        base = snv_alignment.query_alignment_sequence[snv_seq_pos]
                        mut_base = mut_map[base]
                        
                        # print("BASE '%s' MUTATED TO '%s'" % (base, mut_base))
                        # time.sleep(1)
                        
                        if base != mut_base:
                            # base has been mutated to another base
                            # print("base '%s' mut_base '%s'" % (base, mut_base))
                            
                            # print(type(snv_alignment.query_sequence))
                            
                            mut_pos = snv_seq_pos + snv_alignment.query_alignment_start
                            query_sequence = str(snv_alignment.query_sequence[:mut_pos])
                            query_sequence += str(mut_base)
                            query_sequence += str(snv_alignment.query_sequence[mut_pos + 1:])
                            
                            # print("query_sequence before %s" % snv_alignment.query_sequence)
                            # print("query_sequence after  %s" % query_sequence)
                            
                            # print("base pos %d" % mut_pos)
                            # print("seq: %s" % snv_alignment.query_sequence)
                            # print("seq: %s" % query_sequence)
                            # print("qualities: %s" % snv_alignment.query_qualities)
                            
                            snv_alignment.query_sequence = query_sequence
                            
                            # start at -1 to work with 0-based position
                            seq_pos = -1
                            for cig_op_id, cig_op_len in snv_alignment.cigartuples:
                                if cig_op_id == CIGAR_MATCH or cig_op_id == CIGAR_INS:
                                    # match and insert are present in the query_alignment_sequence
                                    curr_seq_pos = seq_pos + cig_op_len
                                elif cig_op_id == CIGAR_DEL:
                                    # deletion is not present in the query_alignment_sequence
                                    continue
                                else:
                                    raise ValueError("Unsupported CIGAR operation %d" % cig_op_id)

                                if curr_seq_pos < snv_seq_pos:
                                    # curr_seq_pos has not reached snv_seq_pos
                                    seq_pos = curr_seq_pos
                                else:
                                    # cigar segment covers snv position, it must be a match
                                    if cig_op_id != CIGAR_MATCH:
                                        # snv is not matched
                                        # TODO update mapping quality
                                        # TODO update cigar string
                                        
                                        print("snv_seq_pos %d" % snv_seq_pos)
                                        print("seq_pos %d" % seq_pos)
                                        print(snv_alignment.cigarstring)
                                        print(snv_alignment.get_aligned_pairs(matches_only=False, with_seq=False))
                                        
                                        raise ValueError("Unexpected CIGAR operation %d" % cig_op_id)
                                    break
                
                # done with this snv, read next
                snv = read_vac_snv(vac_file, fai_list)
                new_snv_alignments = snv_alignments[:]
                for snv_alignment in snv_alignments:
                    # snv alignment is either before or overlapping next snv
                    if is_before_snv(snv_alignment, snv, fai_list):
                        out_file.write(snv_alignment)
                        counter += 1
                        
                        # print("deleted alignment id: %d" % id(snv_alignment))
                        
                        del snv_alignment_map[snv_alignment]
                        new_snv_alignments.remove(snv_alignment)
                    else:  # is overlapping snv
                        # find sequence position of snv
                        snv_seq_pos = ref_pos2seq_pos(snv_alignment, snv['ref_pos'])
                        if snv_seq_pos is not None:
                            # alignment is mapped at another snv position
                            snv_alignment_map[snv_alignment] = snv_seq_pos
                
                # print("snv_alignments size before %d" % len(snv_alignments))
                
                snv_alignments = new_snv_alignments
                
                # print("snv_alignments size after %d" % len(snv_alignments))
                # if len(snv_alignments) > 0:
                #     print([id(align) for align in snv_alignments])
            
            else:  # alignment is overlapping snv
                
                #                 print("alignment is overlapping snv")
                #                 time.sleep(0.1)
                
                # alignment dont have to be mapped to the snv (indel)
                
                # find sequence position of snv
                snv_seq_pos = ref_pos2seq_pos(alignment, snv['ref_pos'])
                
                snv_alignments.append(alignment)
                snv_alignment_map[alignment] = snv_seq_pos
                
                alignment = next_item(sam_file)


if __name__ == "__main__":
    BAM_FILENAME = "/data/projects/varlock/mapping/1020b_S23_chr22.bam"
    FAI_FILENAME = "/data/genome/human/hg19/hg19.fa.fai"
    VAC_FILENAME = "/data/projects/varlock/scripts/in/chr22.vac"
    OUT_FILENAME = "/data/projects/varlock/mapping/1020b_S23_chr22_mut.bam"
    
    # 517362 alignments in BAM_FILENAME
    # python3 /data/projects/varlock/scripts/mutate.py
    
    the_method(vac_filename=VAC_FILENAME,
               bam_filename=BAM_FILENAME,
               out_filename=OUT_FILENAME,
               rnd=random.SystemRandom(),
               verbose=True)
