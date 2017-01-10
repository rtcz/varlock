import random
import struct

import numpy as np
import pysam


class Mutator:
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
    
    UNKNOWN_BASE = "N"
    BASES = ["A", "T", "G", "C"]
    
    def __init__(self, fai_filename, rnd=random.SystemRandom(), verbose=False):
        """
        :param fai_filename:
        :param rnd: (secure) random generator
        :param verbose:
        """
        self.fai_list = self.parse_fai(fai_filename)
        self.rnd = rnd
        self.verbose = verbose
        
        self.alignment_counter = None
        self.mut_counter = None
        self.ns_mut_counter = None
        self.snv_counter = None
        self.unmapped_counter = None
        self.overlapping_counter = None
        self.max_snv_alignments = None
    
    @staticmethod
    def multi_random(p_dist, rnd):
        """
        Draw index from multinomial probability distribution.
        :param p_dist: Array of probabilities. When all probabilities are zero, each outcome has equal probability.
        Each value represents probability of one outcome.
        :param rnd: random number generator
        :return: Always returns index of outcome in p_dist array.
        """
        if len(p_dist) == 1:
            # only one choice
            return 0
        
        if sum(p_dist) == 0:
            # each outcome has the equal probability
            return rnd.randint(0, len(p_dist) - 1)
        
        p_level = 0
        rnd_value = rnd.random()
        
        # make relative
        p_value = rnd_value * sum(p_dist)
        
        for i in range(len(p_dist)):
            p_level += p_dist[i]
            if p_level > p_value:
                # random outcome has been reached
                return i
        
        raise ValueError(
            "sum of probability distribution %d must be greater then the probability value %d" % (sum(p_dist), p_value))
    
    @staticmethod
    def strip_prefix(value, prefix):
        if value.startswith(prefix):
            return value[len(prefix):]
        else:
            return value
    
    @classmethod
    def strip_chr(cls, value):
        """
        Strip 'chr' prefix from string.
        :param value: string ot strip
        :return: stripped string
        """
        return cls.strip_prefix(value, 'chr')
    
    def index2pos(self, index):
        """
        Convert absolute position (index) on genome to reference position.
        :param index: absolute position on genome
        :return: reference position as tuple (<reference name>, <0-based position>)
        """
        ref_pos = index
        for ref in self.fai_list:
            new_ref_pos = ref_pos - ref["length"]
            if new_ref_pos < 0:
                return self.strip_chr(ref["name"]), ref_pos
            else:
                ref_pos = new_ref_pos
        raise ValueError("reference position for index %d not found" % index)
    
    def pos2index(self, ref_name, ref_pos):
        """
        Convert reference position to absolute position (index) on genome.
        :param ref_name: reference name
        :param ref_pos: 0-based reference position
        :return: absolute position on genome
        """
        sum_ref_pos = 0
        for ref in self.fai_list:
            if self.strip_chr(ref["name"]) == self.strip_chr(ref_name):
                return int(ref_pos) + sum_ref_pos
            else:
                sum_ref_pos += ref["length"]
        raise ValueError("sequence name %s not found" % ref_name)
    
    @staticmethod
    def parse_fai(fai_filename):
        """
        Parse fasta index file.
        :param fai_filename: filename
        :return: list of references from fasta index file
        """
        data = []
        ref_keys = ["name", "length", "offset", "linebases"]
        with open(fai_filename, "r") as fai_file:
            for line in fai_file:
                raw_record = line.rstrip().split()
                ref_name = raw_record[0]
                ref_data = map(int, raw_record[1:])
                
                # noinspection PyTypeChecker
                ref = [ref_name] + list(ref_data)
                record_map = dict(zip(ref_keys, ref))
                data.append(record_map)
        
        return data
    
    def is_before_snv(self, alignment, snv):
        """
        Check if alignment is mapped before snv.
        :return: True if alignment end is mapped before snv
        """
        # reference_end points to one past the last aligned residue
        ref_end = alignment.reference_end - 1
        alignment_end = self.pos2index(alignment.reference_name, ref_end)
        return alignment_end < snv['index']
    
    def is_after_snv(self, alignment, snv):
        """
        Check if alignment is mapped after snv.
        :return: True if alignment start is mapped after snv
        """
        alignment_start = self.pos2index(alignment.reference_name, alignment.reference_start)
        return alignment_start > snv['index']
    
    @classmethod
    def count_bases(cls, base_pileup):
        """
        Count allele count in base pileup column.
        :param base_pileup: list of base occurences
        :return: list of DNA base frequencies
        """
        alt_ac = [0] * 4
        for base in base_pileup:
            # skip unknown base
            if base != cls.UNKNOWN_BASE:
                try:
                    alt_ac[cls.BASES.index(base)] += 1
                except KeyError:
                    raise ValueError("Illegal DNA base %s" % base)
        
        return alt_ac
    
    @classmethod
    def create_mut_map(cls, alt_ac, ref_ac, rnd):
        """
        Creates mutation mapping for base pileup column.
        pileup base -> mutated base
        :param alt_ac: list of DNA bases frequencies
        :param ref_ac:
        :param rnd: random number generator
        :return: dict which is specific mutation mapping
        """
        ref_bases = cls.BASES[:]
        alt_bases = cls.BASES[:]
        
        # add random value to distinguish tied values
        alt_ac = [ac + rnd.random() for ac in alt_ac]
        
        # init mutation mapping
        mut_map = dict.fromkeys(cls.BASES, None)
        # unknown base is always mapped to itself
        mut_map[cls.UNKNOWN_BASE] = cls.UNKNOWN_BASE
        # map bases but skip last unmapped base
        for i in range(len(cls.BASES) - 1):
            # draw ref base with multinomial probability
            ref_base_id = cls.multi_random(ref_ac, rnd)
            # draw most abundant base from alt alleles
            alt_base_id = np.argmax(alt_ac)
            # add mapping
            mut_map[alt_bases[alt_base_id]] = ref_bases[ref_base_id]
            
            # delete processed items
            del ref_bases[ref_base_id]
            del ref_ac[ref_base_id]
            del alt_bases[alt_base_id]
            del alt_ac[alt_base_id]
        
        # last base mapping is obvious
        mut_map[alt_bases[0]] = ref_bases[0]
        
        return mut_map
    
    @staticmethod
    def next_item(iterable):
        try:
            return next(iter(iterable))
        except StopIteration:
            return None
    
    def read_vac_snv(self, vac_file):
        byte_string = vac_file.read(12)
        if byte_string == "":
            return None
        data_list = list(struct.unpack('<IHHHH', byte_string))
        ref_name, ref_pos = self.index2pos(data_list[0])
        return {'index': data_list[0], 'ref_name': ref_name, 'ref_pos': ref_pos, 'ac': data_list[1:]}
    
    @staticmethod
    def ref_pos2seq_pos(alignment, ref_pos):
        """
        Retrieve base position in sequence string with
        :param alignment: pysam.AlignedSegment
        :param ref_pos: reference position of base
        :return: 0-based position of base with specified reference position in sequence string
        """
        for seq_pos, current_ref_pos in alignment.get_aligned_pairs(matches_only=True, with_seq=False):
            # search for base in snv position
            if current_ref_pos == ref_pos:
                return seq_pos
        return None
    
    def __write_alignment(self, out_file, alignment):
        out_file.write(alignment)
        self.alignment_counter += 1
        if self.verbose and self.alignment_counter % 10000 == 0:
            print("%d alignments processed" % self.alignment_counter)
    
    @staticmethod
    def get_base_pileup(snv_alignments):
        pileup_col = []
        for snv_alignment in snv_alignments:
            if snv_alignment.snv_pos is not None:
                # alignment is mapped at snv position
                # not sure if query_sequence or query_alignment_sequence
                base = snv_alignment.alignment.query_alignment_sequence[snv_alignment.snv_pos]
                pileup_col.append(base)
        
        return pileup_col
    
    def mutate(self, vac_filename, bam_filename, out_filename):
        """
        Mutate BAM file SNVs.
        :param vac_filename: variant allele count filename
        :param bam_filename: input bam filename
        :param out_filename: output bam filename
        """
        # init counters
        self.alignment_counter = 0  # counts written alignments
        self.mut_counter = 0
        self.ns_mut_counter = 0  # non synonymous mutation counter
        self.snv_counter = 1  # counts read snvs
        
        self.unmapped_counter = 0
        self.overlapping_counter = 0
        self.max_snv_alignments = 0
        
        with pysam.AlignmentFile(bam_filename, "rb") as sam_file, \
                open(vac_filename, "rb") as vac_file, \
                pysam.AlignmentFile(out_filename, "wb", template=sam_file) as out_file:
            
            snv_alignments = []
            alignment = self.next_item(sam_file)
            snv = self.read_vac_snv(vac_file)
            
            while True:
                if snv is None or alignment is None:
                    if self.verbose:
                        if snv is None:
                            print("EOF VAC")
                        else:
                            print("EOF BAM")
                    # write remaining snv alignments
                    for snv_alignment in snv_alignments:
                        self.__write_alignment(out_file, snv_alignment.alignment)
                    # write remaining alignments (in EOF VAC case only)
                    while alignment is not None:
                        self.__write_alignment(out_file, alignment)
                        alignment = self.next_item(sam_file)
                    break
                
                elif alignment.is_unmapped:
                    self.__write_alignment(out_file, alignment)
                    self.unmapped_counter += 1
                
                elif self.is_before_snv(alignment, snv):
                    # print("alignment is before snv")
                    # time.sleep(0.1)
                    self.__write_alignment(out_file, alignment)
                    alignment = self.next_item(sam_file)
                
                elif self.is_after_snv(alignment, snv):
                    self.mutate_snv_alignments(snv_alignments, snv)
                    # done with this snv, read next
                    snv = self.read_vac_snv(vac_file)
                    # process snv_alignments with new snv
                    snv_alignments = self.process_snv_alignments(snv_alignments, snv, out_file)
                
                else:  # alignment is overlapping snv
                    # alignment dont have to be mapped to the snv (indel)
                    # find sequence position of snv
                    snv_seq_pos = self.ref_pos2seq_pos(alignment, snv['ref_pos'])
                    snv_alignments.append(SnvAlignment(alignment, snv_seq_pos))
                    alignment = self.next_item(sam_file)
                    
                    if self.verbose:
                        self.max_snv_alignments = max(len(snv_alignments), self.max_snv_alignments)
                        self.overlapping_counter += 1
            
            if self.verbose:
                print("alignments %d" % self.alignment_counter)
                print("unmapped alignments %d" % self.unmapped_counter)
                print("overlapping alignments %d" % self.overlapping_counter)
                print("max snv alignments %d" % self.max_snv_alignments)
                print("snvs %d" % self.snv_counter)
                print("mutations %d" % self.mut_counter)
                print("ns mutations %d" % self.ns_mut_counter)
    
    def mutate_snv_alignments(self, snv_alignments, snv):
        # get pileup of bases from alignments at snv mapping position
        base_pileup = self.get_base_pileup(snv_alignments)
        alt_ac = self.count_bases(base_pileup)
        mut_map = self.create_mut_map(alt_ac=alt_ac, ref_ac=snv['ac'], rnd=self.rnd)
        
        for snv_alignment in snv_alignments:
            if snv_alignment.snv_pos is not None:
                # alignment has snv to mutate
                self.mutate_alignment(snv_alignment.alignment, snv_alignment.snv_pos, mut_map)
                # mutation can be done only once
                snv_alignment.snv_pos = None
    
    def mutate_alignment(self, alignment, snv_pos, mut_map):
        # there is base to be mutated
        self.mut_counter += 1
        
        # alignment is mapped at snv position
        # not sure if query_sequence or query_alignment_sequence
        base = alignment.query_alignment_sequence[snv_pos]
        mut_base = mut_map[base]
        
        if base != mut_base:
            # base has been mutated to another base
            self.ns_mut_counter += 1
            self.mutate_alignment(alignment, snv_pos, mut_base)
            self.check_cigar_str(alignment)
    
    @staticmethod
    def mutate_alignment(alignment, mut_pos, mut_base):
        """
        Replace base at snv position with mutated base
        :param alignment: pysam.AlignedSegment
        :param mut_pos: 0-based position to mutate
        :param mut_base: mutated base
        :return: mutated sequence string
        """
        mut_pos = mut_pos + alignment.query_alignment_start
        mut_seq = str(alignment.query_sequence[:mut_pos])
        mut_seq += str(mut_base)
        mut_seq += str(alignment.query_sequence[mut_pos + 1:])
        alignment.query_sequence = mut_seq
    
    @classmethod
    def check_cigar_str(cls, snv_alignment):
        """
        Validate cigar operation at snv position.
        TODO: Update cigar and quality strings.
        :param snv_alignment:
        """
        # start at -1 to work with 0-based position
        seq_pos = -1
        for cig_op_id, cig_op_len in snv_alignment.alignment.cigartuples:
            if cig_op_id == cls.CIGAR_MATCH or cig_op_id == cls.CIGAR_INS:
                # match and insert are present in the query_alignment_sequence
                curr_seq_pos = seq_pos + cig_op_len
            elif cig_op_id == cls.CIGAR_DEL:
                # deletion is not present in the query_alignment_sequence
                continue
            else:
                raise ValueError("Unsupported CIGAR operation %d" % cig_op_id)
            
            if curr_seq_pos < snv_alignment.seq_pos:
                # curr_seq_pos has not reached snv_seq_pos
                seq_pos = curr_seq_pos
            else:
                # cigar segment covers snv position, it must be a match
                if cig_op_id != cls.CIGAR_MATCH:
                    # snv is not matched
                    # TODO update mapping quality
                    # TODO update cigar string
                    
                    print("snv_seq_pos %d" % snv_alignment.seq_pos)
                    print("seq_pos %d" % seq_pos)
                    print(snv_alignment.cigarstring)
                    print(snv_alignment.get_aligned_pairs(matches_only=False, with_seq=False))
                    
                    raise ValueError("Unexpected CIGAR operation %d" % cig_op_id)
                break
    
    def process_snv_alignments(self, snv_alignments, new_snv, out_file):
        """
        Update ovelapping alignments with new snv. Save alignments before new snv to file
        :param snv_alignments: snv alignments overlapping previous snvs
        :param new_snv:
        :param out_file:
        :return: list of remaining alignments
        """
        new_snv_alignments = snv_alignments[:]
        for i in range(len(snv_alignments)):
            snv_alignment = snv_alignments[i]
            # new_snv alignment is either before or overlapping next new_snv
            if self.is_before_snv(snv_alignment.alignment, new_snv):
                self.__write_alignment(out_file, snv_alignment.alignment)
                # remove written alignment
                del new_snv_alignments[i]
            else:  # is overlapping new_snv
                # find sequence position of new_snv
                snv_seq_pos = self.ref_pos2seq_pos(snv_alignment.alignment, new_snv['ref_pos'])
                if snv_seq_pos is not None:
                    # alignment is mapped at another new_snv position
                    snv_alignment.snv_pos = snv_seq_pos
        
        return new_snv_alignments


class SnvAlignment:
    def __init__(self, alignment, snv_pos):
        """
        :param alignment: pysam.AlignedSegment
        :param snv_pos: reference position of snv
        """
        self.alignment = alignment
        self.snv_pos = snv_pos


if __name__ == "__main__":
    BAM_FILENAME = "/data/projects/varlock/mapping/1020b_S23_chr22.bam"
    FAI_FILENAME = "/data/genome/human/hg19/hg19.fa.fai"
    VAC_FILENAME = "/data/projects/varlock/scripts/in/chr22.vac"
    OUT_FILENAME = "/data/projects/varlock/mapping/1020b_S23_chr22_mut.bam"
    
    # python3 /data/projects/varlock/scripts/mutate.py
    
    # written alignments 517362
    # unmapped alignments 0
    # overlapping alignments 495465
    # max snv alignments 667
    # read snvs 1059416, 101 omitted
    # total snvs 1059517
    
    mutator = Mutator(fai_filename=FAI_FILENAME, verbose=True)
    mutator.mutate(vac_filename=VAC_FILENAME, bam_filename=BAM_FILENAME, out_filename=OUT_FILENAME)
