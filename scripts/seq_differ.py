from bitarray import bitarray
from random import SystemRandom


class SeqDiffer:
    BASE_2_BIN = {'A': '00', 'T': '01', 'G': '10', 'C': '11'}
    BIN_2_BASE = {'00': 'A', '01': 'T', '10': 'G', '11': 'C'}
    
    BASES = ['A', 'T', 'G', 'C']
    
    # variant bit size
    VAR_BIT_SIZE = 16
    VAR_POS_BIT_SIZE = 14
    # MAX_VARIANT_POS = 2**VAR_POS_BIT_SIZE-1
    
    def __init__(self):
        self.rnd = SystemRandom()
    
    @staticmethod
    def diff(seq, mut_seq):
        """
        Compact sequence diff.
        :param seq: Original sequence.
        :param mut_seq: Mutated sequence.
        :return: Diff as bitarray, preserving values of original sequence.
        """
        diff = bitarray()
        abs_pos = 0
        if len(seq) != len(mut_seq):
            raise Exception("sequences must have same length")
        for i in range(len(seq)):
            if seq[i] != mut_seq[i]:
                # variant found
                bin_base = SeqDiffer.BASE_2_BIN[seq[i]]
                diff += SeqDiffer.var2bytes(i - abs_pos, bin_base)
                abs_pos = i
        
        return diff
    
    @staticmethod
    def var2bytes(pos, base, max_pos_bits=VAR_POS_BIT_SIZE):
        """
        :param pos: Relative position.
        :param base: Base as binary string.
        :param max_pos_bits: Max position bin length
        :return: Encoded variant in 2 bytes (16 bits).
        """
        variant = bitarray()
        # use position relative to last difference
        bin_pos = bitarray(bin(pos)[2:])
        if bin_pos.length() > max_pos_bits:
            raise Exception("position is too big %d" % pos)
        
        padding_size = (max_pos_bits - bin_pos.length())
        bin_padding = bitarray('0' * padding_size)
        
        # add variant
        variant += bin_padding + bin_pos
        variant += base
        
        return variant
    
    def mutate(self, seq, mut_p=0, mut_map_p={}):
        """
        Randomly mutate sequence with defined probability.
        :param seq:
        :param mut_p:
        :param mut_map_p:
        :return: Mutated sequence.
        """
        mut_seq = ""
        for i in range(len(seq)):
            base = seq[i]
            p = self.rnd.random()
            if i in mut_map_p:
                # specific probability distribution
                mut_seq += SeqDiffer.BASES[self.multi_random(mut_map_p[i])]
            elif mut_p > p:
                mut_seq += self.rnd.choice(SeqDiffer.BASES)
            else:
                mut_seq += base
        
        return mut_seq
    
    @staticmethod
    def restore(mut_seq, diff):
        """
        Restore original sequence from mutated sequence.
        :param mut_seq: Mutated sequence.
        :param diff: Bitarray representing the diff, containing values of the original sequence.
        :return: original sequence
        """
        seq = mut_seq
        diff_count = diff.length() / SeqDiffer.VAR_BIT_SIZE
        absolute_pos = 0
        for i in range(0, diff_count):
            start = i * SeqDiffer.VAR_BIT_SIZE
            pos_bits = diff[start:start + SeqDiffer.VAR_POS_BIT_SIZE]
            base_bits = diff[start + SeqDiffer.VAR_POS_BIT_SIZE:start + SeqDiffer.VAR_BIT_SIZE]
            
            absolute_pos += int(pos_bits.to01(), 2)
            seq = seq[:absolute_pos] + SeqDiffer.BIN_2_BASE[base_bits.to01()] + seq[absolute_pos + 1:]
        
        return seq

    def multi_random(self, p_dist):
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
            return self.rnd.integer(0, len(p_dist) - 1)
    
        p_level = 0
        rnd_value = self.rnd.random()
        p_value = rnd_value * sum(p_dist)
    
        #     print(p_dist)
        #     print("random %f" % rnd_value)
        #     print("sum(p_dist) %f" % sum(p_dist))
        #     print("p_value %f" % p_value)
        #     print("------")
    
        for i in range(len(p_dist)):
            p_level += p_dist[i]
            if p_value < p_level:
                return i
    
        raise ValueError(
            "sum of probability distribution %d must be greater then the probability value %d" % (sum(p_dist), p_value))
