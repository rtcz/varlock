import pyfaidx
import subprocess
import os


class Anonymer:
	FLAG_UNMAPPED = 4
	FLAG_REVERSE_COMPLEMENT = 16
	FLAG_FIRST = 64
	FLAG_LAST = 128
	
	TEMP_SUFFIX = ".temp"
	FINAL_SUFFIX = ".anon.fastq"
	
	def __init__(self,
	             r1,
	             r2,
	             out_dir,
	             gen_file,
	             gen_idx,
	             verbose=False,
	             threads=1):
		"""
		:param r1: r1 fastq file
		:param r2: r2 fastq file
		:param out_dir: output directory
		:param gen_file: genome file .fa
		:param gen_idx: genome index file
		:param verbose: print what is being done
		:param threads: number of working threads
		"""
		self.genome = pyfaidx.Fasta(gen_file, sequence_always_upper=True)
		self.r1 = r1
		self.r2 = r2
		self.gen_idx = gen_idx
		self.out_dir = out_dir
		self.threads = threads
		self.verbose = verbose
	
	def substitute(self):
		"""
		substitutes uniquely mapped sequences with reference sequences
		"""
		# sam_file = os.path.join(self.out_dir, 'sam.sam')
		
		TEMP_FIRST_SUFFIX = self.TEMP_SUFIX + "1"
		TEMP_SECOND_SUFFIX = self.TEMP_SUFFIX + "2"
		
		if self.verbose:
			print "first substitution"
		self.__substitute_fastq_pair(in_suffix="", out_suffix=TEMP_FIRST_SUFFIX)
		
		if self.verbose:
			print "second substitution"
		self.__substitute_fastq_pair(in_suffix=TEMP_FIRST_SUFFIX, out_suffix=TEMP_SECOND_SUFFIX)
		
		if self.verbose:
			print "reverting non unique mappings"
		
		self.r1 + TEMP_FIRST_SUFFIX
		self.r1 + TEMP_SECOND_SUFFIX
		
		self.r2 + TEMP_FIRST_SUFFIX
		self.r2 + TEMP_SECOND_SUFFIX
		self.__save_subsitution(TEMP_FIRST_SUFFIX)
		
		if self.verbose:
			print "anononymization finished"
	
	def __in2out_fastq(self, fastq_file):
		return os.path.join(self.out_dir, self.basic_basename(fastq_file) + self.OUT_POSTFIX)
	
	def __save_subsitution(self, orig_fastq, subst1_fastq, subst2_fastq):
		out_fastq = orig_fastq
		with open(orig_fastq, 'r') as orig_file, \
				open(subst1_fastq, 'r') as subst1_file, \
				open(subst2_fastq, 'r') as subst2_file, \
				open(out_fastq) as out_file:
			for orig_line in orig_file:
				subst1_line = subst1_file.readline()
				subst2_line = subst2_file.readline()
				if subst1_line != subst2_line:
					# sequence is not uniquely mapped, use original sequence
					out_file.write(orig_line)
				else:
					# sequence is uniquely mapped, use it
					out_file.write(subst1_line)
				break
	
	def __substitute_fastq_pair(self, in_suffix, out_suffix):
		r1_in = self.r1 + in_suffix
		r2_in = self.r2 + in_suffix
		r1_out = self.r1 + out_suffix
		r2_out = self.r2 + out_suffix
		
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
			self.__substitute_fastq_line_to_file(bitwise_flag=bitwise_flag,
			                                     chr_id=alignment[2],
			                                     chr_pos=alignment[3],
			                                     seq=alignment[9],
			                                     seq_qual=alignment[10],
			                                     out_file=r1_file,
			                                     seq_name=seq_name + "/1")
		else:
			# last member of a fastq pair
			self.__substitute_fastq_line_to_file(bitwise_flag=bitwise_flag,
			                                     chr_id=alignment[2],
			                                     chr_pos=alignment[3],
			                                     seq=alignment[9],
			                                     seq_qual=alignment[10],
			                                     out_file=r2_file,
			                                     seq_name=seq_name + "/2")
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
		command.append('-1 ' + self.r1)
		command.append('-2 ' + self.r2)
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
			chr_pos_end = chr_pos_start + len(seq) - 1
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
	def basic_basename(path):
		"""
		example: "/a/b/file.fastq.gz" -> "file"
		:param path: file name / path to file
		:return: first segment of file
		"""
		return os.path.basename(path).split(".")[0]
