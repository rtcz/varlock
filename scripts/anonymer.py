import pyfaidx
import subprocess
import os
import re


class Anonymer:
    FLAG_UNMAPPED = 4
    FLAG_REVERSE_COMPLEMENT = 16
    FLAG_FIRST = 64
    FLAG_LAST = 128
    
    SEQ_VALUE_LINE_ID = 1  # second line
    
    SAM_SEQ_NAME_ID = 0
    SAM_FLAG_ID = 1
    SAM_CHR_ID = 2
    SAM_POS_ID = 3
    SAM_MAPQ_ID = 4
    SAM_SEQ_ID = 9
    SAM_QUAL_ID = 10
    
    SUBST_1_SUFFIX = ".subst1"
    SUBST_2_SUFFIX = ".subst2"
    OUT_FASTQ_SUFFIX = ".anon.fastq"
    OUT_SAM_SUFFIX = ".sam"
    
    FASTQ_SUFFIX_PATTERN = "[-_]?(1,2)?[Rr]?[12]?\.fastq(.gz)?$"
    
    def __init__(
            self,
            r1,
            r2,
            out_dir,
            gen_file,
            gen_idx,
            verbose=False,
            threads=1,
            keep_temp=False):
        """
        :param r1: r1 fastq file
        :param r2: r2 fastq file
        :param out_dir: output directory
        :param gen_file: genome file .fa
        :param gen_idx: genome index file
        :param verbose: explain what is being done
        :param threads: number of working threads
        :param keep_temp: keep temporary substitution files
        """
        self.genome = pyfaidx.Fasta(gen_file, sequence_always_upper=True)
        self.r1 = r1
        self.r2 = r2
        self.gen_idx = gen_idx
        self.out_dir = out_dir
        self.threads = threads
        self.verbose = verbose
        self.keep_temp = keep_temp
    
    def substitute(self):
        """
        substitutes uniquely mapped sequences with reference sequences
        """
        # sam_file = os.path.join(self.out_dir, 'sam.sam')
        
        sample_name = self.__sample_name(self.r1)
        
        if self.verbose:
            print "> first substitution"
        self.__substitute_fastq_pair(
            r1_in=self.r1,
            r2_in=self.r2,
            r1_out=self.__out_file_path(self.r1, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX),
            r2_out=self.__out_file_path(self.r2, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX),
            sam_out=self.__out_file_path(sample_name, self.SUBST_1_SUFFIX + self.OUT_SAM_SUFFIX)
        )
        
        if self.verbose:
            print "> second substitution"
        self.__substitute_fastq_pair(
            r1_in=self.__out_file_path(self.r1, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX),
            r2_in=self.__out_file_path(self.r2, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX),
            r1_out=self.__out_file_path(self.r1, self.SUBST_2_SUFFIX + self.OUT_FASTQ_SUFFIX),
            r2_out=self.__out_file_path(self.r2, self.SUBST_2_SUFFIX + self.OUT_FASTQ_SUFFIX),
            sam_out=self.__out_file_path(sample_name, self.SUBST_2_SUFFIX + self.OUT_SAM_SUFFIX)
        )
        
        if self.verbose:
            print "> saving anonymized fastq pair with unique mappings"
        self.__save_substitution(
            r1_in=self.r1,
            r2_in=self.r2,
            r1_subst1=self.__out_file_path(self.r1, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX),
            r2_subst1=self.__out_file_path(self.r2, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX),
            sam_subst1=self.__out_file_path(sample_name, self.SUBST_1_SUFFIX + self.OUT_SAM_SUFFIX),
            r1_subst2=self.__out_file_path(self.r1, self.SUBST_2_SUFFIX + self.OUT_FASTQ_SUFFIX),
            r2_subst2=self.__out_file_path(self.r2, self.SUBST_2_SUFFIX + self.OUT_FASTQ_SUFFIX),
            sam_subst2=self.__out_file_path(sample_name, self.SUBST_2_SUFFIX + self.OUT_SAM_SUFFIX),
            r1_out=self.__out_file_path(self.r1, self.OUT_FASTQ_SUFFIX),
            r2_out=self.__out_file_path(self.r2, self.OUT_FASTQ_SUFFIX)
        )
        
        if not self.keep_temp:
            if self.verbose:
                print "> deleting temporary files"
            os.remove(self.__out_file_path(self.r1, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX))
            os.remove(self.__out_file_path(self.r2, self.SUBST_1_SUFFIX + self.OUT_FASTQ_SUFFIX))
            os.remove(self.__out_file_path(self.r1, self.SUBST_2_SUFFIX + self.OUT_FASTQ_SUFFIX))
            os.remove(self.__out_file_path(self.r2, self.SUBST_2_SUFFIX + self.OUT_FASTQ_SUFFIX))
            os.remove(self.__out_file_path(sample_name, self.SUBST_1_SUFFIX + self.OUT_SAM_SUFFIX))
            os.remove(self.__out_file_path(sample_name, self.SUBST_2_SUFFIX + self.OUT_SAM_SUFFIX))
        
        if self.verbose:
            print "> anonymization finished"
    
    def __out_file_path(self, file_name, suffix):
        return self.out_dir + "/" + \
               self.first_basename_segment(file_name) + \
               suffix
    
    def __save_substitution(
            self,
            r1_in,
            r2_in,
            r1_subst1,
            r2_subst1,
            sam_subst1,
            r1_subst2,
            r2_subst2,
            sam_subst2,
            r1_out,
            r2_out):
        with \
                open(r1_in, 'r') as r1_in_file, \
                open(r2_in, 'r') as r2_in_file, \
                open(r1_subst1, 'r') as r1_subst1_file, \
                open(r2_subst1, 'r') as r2_subst1_file, \
                open(sam_subst1, 'r') as sam_subst1_file, \
                open(r1_subst2, 'r') as r1_subst2_file, \
                open(r2_subst2, 'r') as r2_subst2_file, \
                open(sam_subst2, 'r') as sam_subst2_file, \
                open(r1_out, 'w') as r1_out_file, \
                open(r2_out, 'w') as r2_out_file:
            counter = 0
            for r1_in_line in r1_in_file:
                r2_in_line = r2_in_file.readline()
                
                r1_subst1_line = r1_subst1_file.readline()
                r2_subst1_line = r2_subst1_file.readline()
                
                # TODO use input fastq and bams only
            
                r1_subst2_line = r1_subst2_file.readline()
                r2_subst2_line = r2_subst2_file.readline()
                
                seq_line_id = counter % 4
                if seq_line_id == self.SEQ_VALUE_LINE_ID:
                    # line with sequence letters
                    
                    # if self.is_sam_flag(sam_line_1[self.SAM_FLAG_ID], self.FLAG_FIRST):
                    #     r1_sam_line = sam_line_1
                    #     r2_sam_line = sam_line_2
                    # else:
                    #     r1_sam_line = sam_line_2
                    #     r2_sam_line = sam_line_1

                    sam_subst1_r1, sam_subst1_r2 = self.__sam_pair_2_arr(
                        sam_subst1_file.readline(),
                        sam_subst1_file.readline()
                    )
                    sam_subst2_r1, sam_subst2_r2 = self.__sam_pair_2_arr(
                        sam_subst2_file.readline(),
                        sam_subst2_file.readline()
                    )
                    
                    # if sam_1_line[self.SAM_MAPQ_ID] != sam_2_line[self.SAM_MAPQ_ID]:
                    #     raise ValueError("Unequal mapq between sam files at line " + str(counter / 4 + 1))

                    # sequence is uniquely mapped, (mapping has not changed between subsitutions)
                    # is_uniquely_mapped = r1_subst1_line == r1_subst2_line and r2_subst1_line == r2_subst2_line

                    is_same_chr_r1 = sam_subst1_r1[self.SAM_CHR_ID] == sam_subst2_r1[self.SAM_CHR_ID]
                    is_same_pos_r1 = sam_subst1_r1[self.SAM_POS_ID] == sam_subst2_r1[self.SAM_POS_ID]

                    is_same_chr_r2 = sam_subst1_r2[self.SAM_CHR_ID] == sam_subst2_r2[self.SAM_CHR_ID]
                    is_same_pos_r2 = sam_subst1_r2[self.SAM_POS_ID] == sam_subst2_r2[self.SAM_POS_ID]
                    
                    is_same_chr = is_same_chr_r1 and is_same_chr_r2
                    is_same_pos = is_same_pos_r1 and is_same_pos_r2
                    
                    # sufficient mapping quality
                    # r1_quality = int(sam_subst1_r1[self.SAM_MAPQ_ID])
                    # r2_quality = int(sam_subst1_r2[self.SAM_MAPQ_ID])
                    # has_sufficient_mapq = r1_quality > 1 and r2_quality > 1
                    
                    if is_same_chr and is_same_pos: # and is_uniquely_mapped: # and has_sufficient_mapq:
                        # sequence is uniquely mapped, (mapping has not changed between subsitutions) use it
                        r1_out_file.write(r1_subst1_line)
                        r2_out_file.write(r2_subst1_line)
                    else:
                        # sequence is not uniquely mapped, use original sequence
                        r1_out_file.write(r1_in_line)
                        r2_out_file.write(r2_in_line)
                else:
                    # other than sequence letter lines are same as in the original fastq
                    r1_out_file.write(r1_in_line)
                    r2_out_file.write(r2_in_line)
                counter += 1
    
    def __sam_pair_2_arr(self, first_line, second_line):
        first_arr = self.__sam_line_2_arr(first_line)
        second_arr = self.__sam_line_2_arr(second_line)
    
        if self.is_sam_flag(first_arr[self.SAM_FLAG_ID], self.FLAG_FIRST):
            r1_arr = first_arr
            r2_arr = second_arr
        else:
            r1_arr = second_arr
            r2_arr = first_arr
            
        return r1_arr, r2_arr
    
    def __save_subsitution(self, in_fastq, subst1_fastq, subst2_fastq, out_fastq):
        with open(in_fastq, 'r') as in_file, \
                open(subst1_fastq, 'r') as subst1_file, \
                open(subst2_fastq, 'r') as subst2_file, \
                open(out_fastq, 'w') as out_file:
            counter = 0
            for in_line in in_file:
                seq_line_id = counter % 4
                
                subst1_line = subst1_file.readline()
                subst2_line = subst2_file.readline()
                
                if seq_line_id == self.SEQ_VALUE_LINE_ID:
                    # line with sequence letters
                    if subst1_line == subst2_line:
                        # sequence is uniquely mapped, (it has not changed between subsitutions) use it
                        out_file.write(subst1_line)
                    else:
                        # sequence is not uniquely mapped, use original sequence
                        out_file.write(in_line)
                else:
                    # other than sequence letter lines are same as in the original fastq
                    out_file.write(in_line)
                counter += 1
    
    def __substitute_fastq_pair(self, r1_in, r2_in, r1_out, r2_out, sam_out):
        with open(r1_out, 'w') as r1_out_file, open(r2_out, 'w') as r2_out_file, open(sam_out, 'w') as sam_file:
            self.counter = 0
            command = self.__build_bowtie_command(r1_in, r2_in)
            proc = subprocess.Popen(command, stdout=subprocess.PIPE)
            if self.verbose:
                print "bowtie2 process started"
            while True:
                line = proc.stdout.readline()
                sam_file.write(line)
                if line != '':
                    self.__substitute_fastq_line(line, r1_out_file, r2_out_file)
                else:
                    break
            r1_out_file.close()
            r2_out_file.close()
            # sam_file.close()
    
    @staticmethod
    def __sam_line_2_arr(line):
        return line.rstrip().split('\t')
    
    def __substitute_fastq_line(self, line, r1_file, r2_file):
        alignment = self.__sam_line_2_arr(line)
        seq_name = alignment[self.SAM_SEQ_NAME_ID]
        bitwise_flag = alignment[self.SAM_FLAG_ID]
        if self.is_sam_flag(bitwise_flag, self.FLAG_FIRST):
            # first member of a fastq pair
            self.__substitute_fastq_line_to_file(
                bitwise_flag=bitwise_flag,
                chr_id=alignment[self.SAM_CHR_ID],
                chr_pos=alignment[self.SAM_POS_ID],
                mapq=alignment[self.SAM_MAPQ_ID],
                seq=alignment[self.SAM_SEQ_ID],
                seq_qual=alignment[self.SAM_QUAL_ID],
                out_file=r1_file,
                seq_name=seq_name + "/1"
            )
        else:
            # last member of a fastq pair
            self.__substitute_fastq_line_to_file(
                bitwise_flag=bitwise_flag,
                chr_id=alignment[self.SAM_CHR_ID],
                chr_pos=alignment[self.SAM_POS_ID],
                mapq=alignment[self.SAM_MAPQ_ID],
                seq=alignment[self.SAM_SEQ_ID],
                seq_qual=alignment[self.SAM_QUAL_ID],
                out_file=r2_file,
                seq_name=seq_name + "/2"
            )
        self.counter += 1
        if self.verbose and self.counter % 2000 == 0:
            print "%d read pairs processed" % (self.counter / 2)
    
    def __build_bowtie_command(self, r1, r2):
        """
        builds command as array for subprocess.Popen
        :param r1: r1 fastq file
        :param r2: r2 fastq file
        :return: array
        """
        command = []
        command.append('bowtie2')
        command.append('-x ' + self.gen_idx)
        command.append('-1 ' + r1)
        command.append('-2 ' + r2)
        command.append('--very-sensitive')
        command.append('--reorder')
        command.append('--no-head')
        command.append('--quiet')
        command.append('--threads ' + str(self.threads))
        return command
    
    def __substitute_fastq_line_to_file(self, bitwise_flag, seq, seq_qual, mapq, chr_id, chr_pos, out_file, seq_name):
        if not self.is_sam_flag(bitwise_flag, self.FLAG_UNMAPPED):
            # is mapped sequence, substitute it with reference sequence
            chr_pos_start = int(chr_pos) - 1
            chr_pos_end = chr_pos_start + len(seq)
            seq = self.genome[chr_id][chr_pos_start:chr_pos_end].seq
        
        # sequence identifier
        out_file.write("@" + seq_name + "\n")
        if self.is_sam_flag(bitwise_flag, self.FLAG_REVERSE_COMPLEMENT):
            # reversed complemented sequence
            out_file.write(self.reverse_complement(seq) + "\n")
            out_file.write("+" + mapq + "\n")
            out_file.write(seq_qual[::-1] + "\n")
        else:
            # straight sequence
            out_file.write(seq + "\n")
            out_file.write("+\n")
            out_file.write(seq_qual + "\n")
    
    def __sample_name(self, fastq):
        """
        extracts name from fastq file, keeps file path
        example: path/file_R1.fasq -> path/file
        """
        return re.sub(self.FASTQ_SUFFIX_PATTERN, "", fastq)
    
    @staticmethod
    def is_sam_flag(bitwise_flag, bit_value):
        return int(bitwise_flag) & int(bit_value) == bit_value
    
    @staticmethod
    def reverse_complement(seq):
        complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return "".join([complements[base] for base in reversed(seq)])
    
    @staticmethod
    def first_basename_segment(path):
        """
        example: "/a/b/file.fastq.gz" -> "file"
        :param path: file name / path to file
        :return: first segment of file
        """
        return os.path.basename(path).split(".")[0]
