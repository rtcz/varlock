#!/usr/bin/python

# /data/projects/dante/scripts/arguments/input.py
# /data/projects/dante/scripts/parser/readfile.py

# python anonymizer.py -1 "../reads/original/TEST_10000_R1.fastq" -2 "../reads/original/TEST_10000_R2.fastq" -o "../reads/original/" -v

import re
import os
import argparse
import subprocess

import pyfaidx

# TODO
GENOME_FILE = "/data/genome/human/hg38/hg38.fa"
GENOME_INDEX = "/data/genome/human/hg38/bowtie2_index/hg38"
THREAD_COUNT = 8
# reverse complement bitwise flag value
FLAG_REVERSE_COMPLEMENT = 16
FLAG_FIRST = 64
FLAG_LAST = 128
OUT_POSTFIX = ".anon.fastq"


def main():
	args = parse_args()
	r1_filename = os.path.join(args.out_dir, basic_basename(args.r1) + OUT_POSTFIX)
	r2_filename = os.path.join(args.out_dir, basic_basename(args.r2) + OUT_POSTFIX)
	
	# fasta = pyfaidx.Fasta(GENOME_FILE)
	sam_file = os.path.join(args.out_dir, 'sam.sam')
	
	with open(r1_filename, 'w') as r1_file, open(r2_filename, 'w') as r2_file, open(sam_file, 'w') as sam_file:
		counter = 0
		proc = subprocess.Popen(build_bowtie_command(args), stdout=subprocess.PIPE)
		if args.verbose:
			print "bowtie2 process started"
		while True:
			line = proc.stdout.readline()
			alignment = line.strip().split('\t')
			if line != '':
				sam_file.write(line)
				if is_sam_flag(alignment[1], FLAG_FIRST):
					# first member of a fastq pair
					alignment2fastq_file(alignment, r1_file, 1)
				else:
					# last member of a fastq pair
					alignment2fastq_file(alignment, r2_file, 2)
				counter += 1
				if args.verbose and counter % 2000 == 0:
					print "%d read pairs processed" % (counter / 2)
			else:
				break
		r1_file.close()
		r2_file.close()
		sam_file.close()


def build_bowtie_command(args):
	command = []
	command.append('bowtie2')
	command.append('-x ' + GENOME_INDEX)
	command.append('-1 ' + args.r1)
	command.append('-2 ' + args.r2)
	command.append('--very-sensitive')
	command.append('--reorder')
	# command.append('--no-sq')
	command.append('--no-head')
	command.append('--quiet')
	command.append('--threads ' + str(THREAD_COUNT))
	
	# TODO out to file (stderr vs stdout ?)
	# command.append('2>>2')
	
	return command


def is_sam_flag(bitwise_flag, bit_value):
	return int(bitwise_flag) & int(bit_value) == bit_value


def alignment2fastq_file(alignment, out_file, file_id):
	# sequence identifier
	out_file.write("@" + alignment[0] + "/" + str(file_id) + "\n")
	if is_sam_flag(alignment[1], FLAG_REVERSE_COMPLEMENT):
		# reversed complemented sequence
		out_file.write(reverse_complement(alignment[9]) + "\n")
		out_file.write("+\n")
		out_file.write(alignment[10][::-1] + "\n")
	else:
		# straight sequence
		out_file.write(alignment[9] + "\n")
		out_file.write("+\n")
		out_file.write(alignment[10] + "\n")


def reverse_complement(seq):
	complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
	return "".join([complements[base] for base in reversed(seq)])


def parse_args():
	parser = argparse.ArgumentParser()
	required = parser.add_argument_group("Required")
	required.add_argument('-1', '--r1', type=is_fastq_file, help="read 1 fastq file", required=True)
	required.add_argument('-2', '--r2', type=is_fastq_file, help="read 2 fastq file", required=True)
	
	optional = parser.add_argument_group("Optional")
	optional.add_argument('-o', '--out-dir', type=is_dir, help="output directory", default=os.getcwd())
	optional.add_argument('-t', '--threads', type=is_pos_int, help="number of working threads", default=THREAD_COUNT)
	optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
	# optional.add_argument('-s', '--sam', help="sam file name for saving sam file")
	return parser.parse_args()


def is_dir(value):
	if os.path.isdir(value):
		return os.path.realpath(value)
	else:
		raise argparse.ArgumentTypeError("Value %s is not a directory." % value)


def basic_basename(path):
	"""
	example: "/a/b/file.fastq.gz" -> "file"
	:param path: file name / path to file
	:return: first segment of file
	"""
	return os.path.basename(path).split(".")[0]


def is_pos_int(value):
	pos_int_pattern = re.compile("^[1-9]\d*$")
	if pos_int_pattern.match(value):
		return int(value)
	else:
		raise argparse.ArgumentTypeError("Value %s is not a non-zero positive integer." % value)


def is_fastq_file(value):
	fastq_pattern = re.compile(".*\.(fastq|fastq\.gz)?$")
	if os.path.isfile(value) and fastq_pattern.match(value):
		return value
	else:
		raise argparse.ArgumentTypeError("Value %s is not a fastq file." % value)


if __name__ == "__main__":
	main()
