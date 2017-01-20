import pysam

BASES = ("A", "T", "G", "C")
UNKNOWN_BASE = "N"


class SnvAlignment:
    def __init__(self, alignment, snv_pos):
        """
        :param alignment: pysam.AlignedSegment
        :param snv_pos: reference position of snv
        """
        self.alignment = alignment
        self.snv_pos = snv_pos


class Snv:
    def __init__(self, index, ref_name, ref_pos, ac):
        self.index = index
        self.ref_name = ref_name
        self.ref_pos = ref_pos
        self.ac = ac


class FaiReference:
    def __init__(self, ref_id, name, start, length):
        """
        :param ref_id: reference id
        :param name: reference name
        :param start: 0-based first position
        :param length: reference length - number of bases
        """
        self.ref_id = ref_id
        self.name = name
        self.start = start
        self.length = length


def parse_fai(fai_filename):
    """
    Parse fasta index file.
    :param fai_filename: filename
    :return: list of references from fasta index file
    """
    data = []
    with open(fai_filename, "r") as fai_file:
        counter = 0
        start = 0
        for line in fai_file:
            raw_record = line.rstrip().split()
            name = strip_chr(raw_record[0])
            length = int(raw_record[1])
            
            reference = FaiReference(ref_id=counter, name=name, start=start, length=length)
            data.append(reference)
            counter += 1
            start += reference.length
    
    return data


def fai_list2dict(fai_list):
    """
    :param fai_list: parsed FAI as list
    :return:
    """
    return dict((reference.name, reference) for reference in fai_list)


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
        # each outcome has equal probability
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


def index2pos(index, fai_list):
    """
    Convert absolute position (index) on genome to reference position.
    :param index: 0-based position on genome
    :param fai_list: parsed fai as list
    :return: reference position as tuple (<reference name>, <0-based position>)
    """
    ref_pos = index
    for ref in fai_list:
        new_ref_pos = ref_pos - ref.length
        if new_ref_pos < 0:
            return ref.name, ref_pos
        else:
            ref_pos = new_ref_pos
    raise ValueError("reference position for index %d not found" % index)


def pos2index(ref_name, ref_pos, fai_dict):
    """
    Convert reference position to absolute position (index) on genome.
    :param ref_name: reference name
    :param ref_pos: 0-based reference position
    :param fai_dict: parsed fai as dict
    :return: 0-based position on genome
    """
    return fai_dict[strip_chr(ref_name)].start + ref_pos


def strip_chr(value):
    """
    Strip 'chr' prefix from string.
    :param value: string ot strip
    :return: stripped string
    """
    if value[:3] == 'chr':
        return value[3:]
    else:
        return value


def count_bases(base_pileup):
    """
    Count allele count in base pileup column.
    :param base_pileup: list of base occurences
    :return: list of DNA base frequencies
    """
    alt_ac = [0] * 4
    for base in base_pileup:
        # skip unknown base
        if base != UNKNOWN_BASE:
            try:
                alt_ac[BASES.index(base)] += 1
            except KeyError:
                raise ValueError("Illegal DNA base %s" % base)
    
    return alt_ac


def sam2bam(sam_filename, bam_filename):
    with pysam.AlignmentFile(sam_filename, "r") as sam_file, \
            pysam.AlignmentFile(bam_filename, "wb", template=sam_file) as bam_file:
        for alignment in sam_file:
            bam_file.write(alignment)


def bam2sam(bam_filename, sam_filename):
    with pysam.AlignmentFile(bam_filename, "rb") as bam_file, \
            pysam.AlignmentFile(sam_filename, "w", template=bam_file) as sam_file:
        for alignment in bam_file:
            sam_file.write(alignment)
