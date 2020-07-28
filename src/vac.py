import struct
import typing

import pyfaidx

from src.common import *
from src.fasta_index import FastaIndex
from src.po import FastaSequence


# TODO refactor like file object interface
class Vac:
    # TODO add header with VCF's checksum ?
    """
    Class for handling VAC file. VAC stands for Variant Allele Count.
    VAC is binary file, consisting of two blocks.
    The first block contains SNV records, the second one INDEL records.
    VAC uses two temporary files when writing. These files are deleted after writing is finished.
    
    VAC header:
    number of SNVs, 4B
    number of INDELs, 4B
    
    VAC SNV record:
    index 4B, 0-based absolute genomic position
    reference index 1B, index of reference allele in the A,T,G,C list
    A allele count 4B
    T allele count 4B
    C allele count 4B
    G allele count 4B
    
    VAC INDEL record:
    index 4B, 0-based absolute genomic position
    length 1B, number of alternatives
    list:, of tuples in format (allele_count, allele_sequence)
        reference allele is the first tuple, alternative alleles are ordered as in VCF
        
        INDEL allele count 4B
        INDEL base length 2B, length of INDEL sequence
        list:
            INDEL base sequence, 1B
    """
    VCF_CHROM_ID = 0
    VCF_POS_ID = 1
    VCF_REF_ID = 3
    VCF_ALT_ID = 4
    VCF_INFO_ID = 7

    VCF_COL_SEP = "\t"
    VCF_LIST_SEP = ","
    VCF_INFO_SEP = ";"

    VAC_COL_SEP = "\t"
    VAC_KEYVAL_SEP = ":"
    VAC_ITEM_SEP = ","

    SHORT_SIZE = 2
    INT_SIZE = 4
    # LONG_SIZE = 8

    HEADER_SIZE = INT_SIZE * 2

    MAX_SHORT_INT = 2 ** 16 - 1
    """
    Max value of 2 bytes
    """

    SNV_FORMAT = "<IBIIII"  # int, byte, int, int, int, int
    SNV_RECORD_SIZE = 21  # bytes

    SNV_TEMP_EXT = '.snv.temp'
    INDEL_TEMP_EXT = '.indel.temp'

    def __init__(self, fai: FastaIndex, verbose=False):
        """
        :param fai: FastaIndex
        :param verbose:
        """
        self.fai = fai
        self.current_pos = None
        self.current_chrom = None
        self.verbose = verbose

    @staticmethod
    def is_snv(allele_list: list):
        """
        :param allele_list: variant reference and alternatives
        :return: True if allele_list represents SNP
        """
        # must contain at least one reference and one allele
        assert len(allele_list) >= 2

        return all(allele in BASES for allele in allele_list)

    @staticmethod
    def is_indel(allele_list: list):
        """
        :param allele_list: variant reference and alternatives
        :return: True if allele_list represents INDEL
        """
        # must contain at least one reference and one allele
        assert len(allele_list) >= 2

        seq_found = False
        for allele in allele_list:
            if len(allele) >= 2:
                # at least one allele must have multiple bases
                seq_found = True
            if any([base not in BASES for base in allele]):
                return False

        return seq_found

    @classmethod
    def parse_ac(cls, info_list: list):
        """
        Parse allele count from VCF INFO column.
        :param info_list: list of VCF INFO key=value strings
        :return: list of allele counts
        """
        for item in info_list:
            # find allele count item
            if item[:3] == 'AC=':
                # AC=value1,value2
                return tuple(map(int, item[3:].split(cls.VCF_LIST_SEP)))

        raise ValueError("VCF record is missing INFO AC attribute")

    @staticmethod
    def parse_an(info_list: list):
        """
        Parse allele number from VCF INFO column
        :param info_list: list of VCF INFO key=value strings
        :return: allele number
        """
        for item in info_list:
            # find total allele number
            if item[:3] == 'AN=':
                # AN=value
                return int(item[3:])

        raise ValueError("VCF record is missing INFO AN attribute")

    @staticmethod
    def snv_count_map(allele_list: tuple, count_list: tuple) -> dict:
        """
        :return: SNV count map
        """
        assert len(allele_list) == len(count_list)

        count_map = dict(zip(allele_list, count_list))

        for base in BASES:
            if base not in count_map:
                # add missing base
                count_map[base] = 0

        return count_map

    @staticmethod
    def _indel_count_map(allele_list: tuple, count_list: tuple) -> tuple:
        """
        :return: INDEL count map as list of tuples
        """
        return tuple(zip(count_list, allele_list))

    @classmethod
    def read_header(cls, vac_file):
        return struct.unpack('<II', vac_file.read(cls.HEADER_SIZE))

    @classmethod
    def write_header(cls, vac_file, snv_count: int, indel_count: int):
        vac_file.write(struct.pack('<II', snv_count, indel_count))

    @classmethod
    def read_snv_record(cls, snv_file):
        """
        :param snv_file:
        :return: index, ref_id, allele count tuple
        """
        byte_str = snv_file.read(cls.SNV_RECORD_SIZE)
        if len(byte_str) == 0:
            raise EOFError()
        index, ref_id, a_ac, t_ac, c_ac, g_ac = struct.unpack(cls.SNV_FORMAT, byte_str)

        return index, ref_id, [a_ac, t_ac, c_ac, g_ac]

    @classmethod
    def _write_snv_record(cls, snv_file, index: int, ref_id: int, ac_tuple: tuple):
        """
        :param snv_file:
        :param index: genomic index
        :param ref_id: reference base id in ac_tuple
        :param ac_tuple: counts of bases in format (A,T,G,C)
        :return:
        """
        # print('ws', index, ac_tuple)
        # assert 0 <= ref_id < len(BASES)
        snv_file.write(struct.pack(cls.SNV_FORMAT, index, ref_id, *ac_tuple))

    @classmethod
    def read_indel_record(cls, indel_file):
        """
        :param indel_file:
        :return: genomic index, allele frequencies, allele sequences
        Reference allele is represented by the first item of returned lists.
        """
        byte_str = indel_file.read(cls.INT_SIZE + 1)
        if len(byte_str) == 0:
            raise EOFError()

        # index = struct.unpack('<I', byte_str)[0]
        # length = int.from_bytes(indel_file.read(1), byteorder='little')
        index, length = struct.unpack('<IB', byte_str)
        freqs = [0] * length
        seqs = [''] * length
        for i in range(length):
            freq, base_length = struct.unpack('<IH', indel_file.read(cls.INT_SIZE + cls.SHORT_SIZE))
            seq_byte_size = math.ceil(base_length / 4)
            seq = bytes2seq(indel_file.read(seq_byte_size), base_length)
            freqs[i] = freq
            seqs[i] = seq

        # print('ri', index, freqs, seqs)
        return index, freqs, seqs

    # TODO parameters: index, counts, seqs
    @classmethod
    def _write_indel_record(cls, indel_file, index: int, indel_map: tuple):
        """
        :param indel_file:
        :param index: genomic index
        :param indel_map: descending sorted list of tuples [(SEQ_A,COUNT_A), (SEQ_B:COUNT_B), ...]
        """
        if len(indel_map) < 2 or len(indel_map) > 255:
            raise IndelError('number of INDEL sequences %d out of range <2, 255>' % len(indel_map))

        record = struct.pack('<IB', index, len(indel_map))
        for allele_count, sequence in indel_map:

            if len(sequence) > cls.MAX_SHORT_INT:
                raise IndelError('length of INDEL sequence %d too large (> %d)' % (len(sequence), cls.MAX_SHORT_INT))

            record += struct.pack('<IH', allele_count, len(sequence))
            record += seq2bytes(sequence)

        indel_file.write(record)

    @staticmethod
    def find_reference_allele(reference: pyfaidx.Sequence, allele_list: list, start_pos: int) \
            -> typing.Union[int, None]:
        """
        :param reference: reference sequence starting on the position of variant
        :param allele_list: list of alleles
        :param start_pos: starting position of variant in chromosome
        :return: index of reference allele, None if not found
        """
        for i, allele in enumerate(allele_list):
            try:
                if str(reference[start_pos:start_pos + len(allele)]) == allele:
                    return i
            except ValueError:
                return None

        return None

    def vcf2vac(
            self,
            vcf_file: pysam.VariantFile,
            vac_file: typing.BinaryIO,
            ref_fasta: dict = None,
            skip_indels: bool = False,
            ac_tag: str = 'AC',
            an_tag: str = 'AN'
    ):
        # TODO enable gzipped VCF file as input
        """
        Converts VCF to binary VAC file.
        Can preserve compatibility between different BAM and VCF reference genomes
        by stripping chr prefix via FastaIndex.
        :param vcf_file: input VCF file
        :param vac_file: output VAC file
        :param ref_fasta: fasta file
        :param skip_indels: whether to skip indels and keep only SNPs
        """
        variant_cnt = 0
        snv_count = 0
        indel_count = 0
        incorrect_count = 0

        snv_filename = vac_file.name + self.SNV_TEMP_EXT
        indel_filename = vac_file.name + self.INDEL_TEMP_EXT

        last_index = -1
        same_pos_records = 0
        not_found_in_ref = 0

        with open(snv_filename, 'wb') as snv_file, \
                open(indel_filename, 'wb') as indel_file:

            for sequence in self.fai:  # type: FastaSequence
                if sequence.name not in vcf_file.header.contigs:
                    # TODO warning
                    continue

                for variant in vcf_file.fetch(sequence.name):  # type: pysam.VariantRecord

                    if self.verbose and variant_cnt % 100000 == 0:
                        print(f'{variant_cnt} variants processed; current variant {variant.chrom}:{variant.pos}')

                    variant_cnt += 1

                    if len(variant.alleles) < 2:
                        if self.verbose:
                            print(f'WARNING: incomplete variant at {variant.chrom}:{variant.pos}')
                        continue

                    is_snv = self.is_snv(variant.alleles)
                    is_indel = self.is_indel(variant.alleles)

                    if not is_snv and not is_indel:
                        if self.verbose:
                            print(f'WARNING: unknown variant type at {variant.chrom}:{variant.pos} '
                                  f'with alleles: {variant.alleles}')
                        continue

                    if is_indel:
                        # check consistency of alleles
                        if skip_indels:
                            continue

                        incorrect = False
                        for alt_allele in variant.alts:
                            if variant.ref[0] != alt_allele[0]:
                                incorrect = True
                                break

                        if incorrect:
                            # skip INDEL if first base of an alternative allele differs
                            if self.verbose:
                                print(f'WARNING: first base mismatch at {variant.chrom}:{variant.pos} '
                                      f'for alleles {variant.alleles}')
                            incorrect_count += 1
                            continue

                    # TODO consider unknown bases
                    # TODO refactor this stupid shit
                    index = self.fai.pos2index(variant.chrom, variant.start)

                    # check if the reference is correct (and switch if not)
                    if ref_fasta is not None:
                        reference = ref_fasta[variant.chrom]
                        ref_allele_id = self.find_reference_allele(reference, variant.alleles, variant.start)
                        if ref_allele_id is None:
                            not_found_in_ref += 1
                            if self.verbose:
                                print(f'WARNING: reference allele not found at {variant.chrom}:{variant.pos} '
                                      f'in available alleles {variant.alleles}')
                        elif ref_allele_id != 0:
                            # swap them if it is not first
                            if self.verbose:
                                print(f'WARNING: reference allele mismatch at {variant.chrom}:{variant.pos} '
                                      f'for allele {variant.alleles[ref_allele_id]}')

                    # skip same position variants:
                    if last_index == index:
                        # TODO merge them ?
                        if self.verbose:
                            print(f'WARNING: skipping duplicated position at {variant.chrom}:{variant.pos}')
                        same_pos_records += 1
                        continue

                    # TODO check if AN and AC exist

                    # an_tag = 'AN_nfe'  # 'DP'
                    # ac_tag = 'AC_nfe'  # 'VD'

                    total_count = variant.info[an_tag]
                    alt_counts = variant.info[ac_tag]

                    if isinstance(total_count, tuple):
                        assert len(total_count) == 1
                        total_count = total_count[0]

                    assert isinstance(total_count, int)

                    if isinstance(variant.info[ac_tag], int):
                        alt_counts = (variant.info[ac_tag],)
                    else:
                        assert isinstance(alt_counts, tuple)

                    if sum(alt_counts) == 0:
                        # nothing to do here
                        # print(f'WARNING: skipping zero alternative allele count at {variant.chrom}:{variant.pos}')
                        continue

                    try:
                        count_list = self.allele_counts(alt_counts, total_count)
                    except NegativeRefCountError:
                        print(f'WARNING: skipping negative reference allele count at {variant.chrom}:{variant.pos}')
                        continue

                    if count_list[0] == 0:
                        if self.verbose:
                            print(f'WARNING: zero reference allele count at {variant.chrom}:{variant.pos} '
                                  f'for alleles {variant.alleles} with counts {count_list}')

                    if is_snv:
                        # allele count map
                        count_map = self.snv_count_map(variant.alleles, count_list)
                        count_tuple = (count_map['A'], count_map['T'], count_map['G'], count_map['C'])
                        self._write_snv_record(
                            snv_file=snv_file,
                            index=index,
                            ref_id=BASES.index(variant.ref),
                            ac_tuple=count_tuple
                        )
                        last_index = index
                        snv_count += 1

                    elif is_indel:
                        indel_map = self._indel_count_map(variant.alleles, count_list)
                        try:
                            self._write_indel_record(
                                indel_file=indel_file,
                                index=index,
                                indel_map=indel_map
                            )
                            last_index = index
                            indel_count += 1
                        except IndelError as e:
                            print(f'ERROR: INDEL {variant.chrom}:{variant.pos} cause: {e}')

        self.write_header(vac_file, snv_count, indel_count)

        if self.verbose:
            print("variants read %d" % variant_cnt)
            print("SNVs written %d" % snv_count)
            print("INDELs written %d" % indel_count)
            print("total same_pos_records %d" % same_pos_records)
            print("total not found in reference %d" % not_found_in_ref)
            print("total incorrect starting point %d" % incorrect_count)

        self.__final_merge(vac_file, snv_filename, indel_filename)

    def __final_merge(self, vac_file, snv_filename: str, indel_filename: str):
        if self.verbose:
            print("merging temporary files")

        with open(snv_filename, 'rb') as snv_file:
            # for chunk in iter(lambda: snv_file.read(4096), b""):
            #     vac_file.write(chunk)

            vac_file.write(snv_file.read())

        with open(indel_filename, 'rb') as indel_file:
            # for chunk in iter(lambda: indel_file.read(4096), b""):
            #     vac_file.write(chunk)

            vac_file.write(indel_file.read())

        os.remove(snv_filename)
        os.remove(indel_filename)

    # deprecta
    def parse_allele_counts(self, info_value: str) -> tuple:
        """
        :param info_value: VCF INFO value
        :return: list of allele counts
        [ref_count, allele_1_count, allele_2_count, ...]
        """
        info_list = info_value.split(self.VCF_INFO_SEP)
        # alternative allele count list
        info_ac = self.parse_ac(info_list)
        # total allele number
        info_an = self.parse_an(info_list)

        return self.allele_counts(info_ac, info_an)

    @staticmethod
    def allele_counts(info_ac: tuple, info_an: int) -> tuple:
        ref_count = info_an - sum(info_ac)

        # sanity check
        if ref_count < 0:
            raise NegativeRefCountError()

        return (ref_count,) + info_ac

    @classmethod
    def text2vac(cls, text_filepath, vac_filepath):
        """
        Text genomic indices are 1-based. After conversion to binary they change to 1-based.
        :param text_filepath:
        :param vac_filepath:
        :return:
        """
        with open(vac_filepath, 'wb') as vac_file, \
                open(text_filepath, 'rt') as text_file:

            snv_count = int(text_file.readline().rstrip())
            vac_file.write(struct.pack('<I', snv_count))

            indel_count = int(text_file.readline().rstrip())
            vac_file.write(struct.pack('<I', indel_count))

            # write SNVs
            for i in range(snv_count):
                line = text_file.readline().rstrip()
                if not line:
                    break
                index_str, ref_id_str, ac_str = line.split(cls.VAC_COL_SEP)
                try:
                    ac_tuple = tuple(map(int, ac_str.split(Vac.VAC_ITEM_SEP)))
                except ValueError:
                    message = "SNV records has invalid format - "
                    message += "probably INDEL record at SNV position"
                    raise ValueError(message)

                cls._write_snv_record(
                    snv_file=vac_file,
                    index=int(index_str) - 1,
                    ref_id=int(ref_id_str),
                    ac_tuple=ac_tuple
                )

            # write INDELs
            for i in range(indel_count):
                line = text_file.readline().rstrip()
                if not line:
                    break
                index_str, indel_str = line.split('\t')
                indel_map = []
                for item in indel_str.split(Vac.VAC_ITEM_SEP):
                    allele_count_str, sequence = item.split(Vac.VAC_KEYVAL_SEP)
                    indel_map.append((int(allele_count_str), sequence))

                cls._write_indel_record(
                    indel_file=vac_file,
                    index=int(index_str) - 1,
                    indel_map=indel_map
                )

    @classmethod
    def vac2text(cls, vac_filepath, text_filepath):
        """
        Indices are always 0-based.
        :param vac_filepath:
        :param text_filepath:
        :return:
        """
        with open(vac_filepath, 'rb') as vac_file, \
                open(text_filepath, 'wt') as text_file:
            snv_count, indel_count = cls.read_header(vac_file)

            text_file.write('%d\n' % snv_count)
            text_file.write('%d\n' % indel_count)
            for i in range(snv_count):
                index, ref_id, count_list = cls.read_snv_record(vac_file)
                ac_string = ','.join(map(str, count_list))
                text_file.write('%d\t%d\t%s\n' % (index, ref_id, ac_string))

            for i in range(indel_count):
                index, counts, seqs = cls.read_indel_record(vac_file)

                text_file.write('%d\t' % index)
                indels = []
                for count, seq in zip(counts, seqs):
                    indels.append('%d:%s' % (count, seq))

                text_file.write(','.join(indels))
                text_file.write('\n')


class IndelError(ValueError):
    pass


class NegativeRefCountError(ValueError):
    pass
