import pyfaidx
import subprocess
import os


class Anonymer:
    FLAG_UNMAPPED = 4
    FLAG_REVERSE_COMPLEMENT = 16
    FLAG_FIRST = 64
    FLAG_LAST = 128
    
    SEQ_VALUE_LINE_ID = 1  # second line
    
    TEMP_SUFFIX = ".subst"
    OUT_SUFFIX = ".anon.fastq"
    
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
        
        temp_first_suffix = self.TEMP_SUFFIX + "1"
        temp_second_suffix = self.TEMP_SUFFIX + "2"
        
        if self.verbose:
            print "> first substitution"
        self.__substitute_fastq_pair(
            r1_in=self.r1,
            r2_in=self.r2,
            r1_out=self.r1 + temp_first_suffix,
            r2_out=self.r2 + temp_first_suffix
        )
        
        if self.verbose:
            print "> second substitution"
        self.__substitute_fastq_pair(
            r1_in=self.r1 + temp_first_suffix,
            r2_in=self.r2 + temp_first_suffix,
            r1_out=self.r1 + temp_second_suffix,
            r2_out=self.r2 + temp_second_suffix
        )
        
        if self.verbose:
            print "> saving anonymized fastq pair with unique mappings"
        self.__save_substitution(
            r1_in=self.r1,
            r2_in=self.r2,
            r1_subst1=self.r1 + temp_first_suffix,
            r2_subst1=self.r2 + temp_first_suffix,
            r1_subst2=self.r1 + temp_second_suffix,
            r2_subst2=self.r2 + temp_second_suffix,
            r1_out=self.__out_path(self.r1),
            r2_out=self.__out_path(self.r2)
        )
        
        if not self.keep_temp:
            if self.verbose:
                print "> deleting temporary substitution files"
            os.remove(self.r1 + temp_first_suffix)
            os.remove(self.r2 + temp_first_suffix)
            os.remove(self.r1 + temp_second_suffix)
            os.remove(self.r2 + temp_second_suffix)
        
        if self.verbose:
            print "> anonymization finished"
    
    def __out_path(self, in_file):
        return self.out_dir + "/" + \
               self.first_basename_segment(in_file) + \
               self.OUT_SUFFIX
    
    def __save_substitution(
            self,
            r1_in,
            r2_in,
            r1_subst1,
            r2_subst1,
            r1_subst2,
            r2_subst2,
            r1_out,
            r2_out):
        with \
                open(r1_in, 'r') as r1_in_file, \
                open(r2_in, 'r') as r2_in_file, \
                open(r1_subst1, 'r') as r1_subst1_file, \
                open(r2_subst1, 'r') as r2_subst1_file, \
                open(r1_subst2, 'r') as r1_subst2_file, \
                open(r2_subst2, 'r') as r2_subst2_file, \
                open(r1_out, 'w') as r1_out_file, \
                open(r2_out, 'w') as r2_out_file:
            counter = 0
            for r1_in_line in r1_in_file:
                r2_in_line = r2_in_file.readline()
                r1_subst1_line = r1_subst1_file.readline()
                r2_subst1_line = r2_subst1_file.readline()
                r1_subst2_line = r1_subst2_file.readline()
                r2_subst2_line = r2_subst2_file.readline()

                seq_line_id = counter % 4
                if seq_line_id == self.SEQ_VALUE_LINE_ID:
                    # line with sequence letters
                    if r1_subst1_line == r1_subst2_line and r2_subst1_line == r2_subst2_line:
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
    
    def __substitute_fastq_pair(self, r1_in, r2_in, r1_out, r2_out):
        with open(r1_out, 'w') as r1_out_file, open(r2_out, 'w') as r2_out_file:  # open(sam_file, 'w') as sam_file:
            self.counter = 0
            command = self.__build_bowtie_command(r1_in, r2_in)
            proc = subprocess.Popen(command, stdout=subprocess.PIPE)
            if self.verbose:
                print "bowtie2 process started"
            while True:
                line = proc.stdout.readline()
                if line != '':
                    self.__substitute_fastq_line(line, r1_out_file, r2_out_file)
                else:
                    break
            r1_out_file.close()
            r2_out_file.close()
            # sam_file.close()
    
    def __substitute_fastq_line(self, line, r1_file, r2_file):
        alignment = line.strip().split('\t')
        # sam_file.write(line)
        seq_name = alignment[0]
        bitwise_flag = alignment[1]
        if self.is_sam_flag(bitwise_flag, self.FLAG_FIRST):
            # first member of a fastq pair
            self.__substitute_fastq_line_to_file(
                bitwise_flag=bitwise_flag,
                chr_id=alignment[2],
                chr_pos=alignment[3],
                seq=alignment[9],
                seq_qual=alignment[10],
                out_file=r1_file,
                seq_name=seq_name + "/1"
            )
        else:
            # last member of a fastq pair
            self.__substitute_fastq_line_to_file(
                bitwise_flag=bitwise_flag,
                chr_id=alignment[2],
                chr_pos=alignment[3],
                seq=alignment[9],
                seq_qual=alignment[10],
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
    
    def __substitute_fastq_line_to_file(self, bitwise_flag, seq, seq_qual, chr_id, chr_pos, out_file, seq_name):
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
            out_file.write("+\n")
            out_file.write(seq_qual[::-1] + "\n")
        else:
            # straight sequence
            out_file.write(seq + "\n")
            out_file.write("+\n")
            out_file.write(seq_qual + "\n")
    
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
