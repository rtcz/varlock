#!/usr/bin/python

# /data/projects/dante/scripts/arguments/input.py
# /data/projects/dante/scripts/parser/readfile.py

# python anonymize.py -1 "../reads/original/TEST_10000_R1.fastq" -2 "../reads/original/TEST_10000_R2.fastq" -o "../reads/original/" -v -k

# bowtie2 -x /data/genome/human/hg38/bowtie2_index/hg38 -1 ../reads/original/TEST_10000_R1.fastq -2 ../reads/original/TEST_10000_R2.fastq -S ../reads/original/TEST_10000.sam --very-sensitive --reorder --no-head --threads 8

import argparse
import os
import re
import sys
from datetime import datetime

from src.anonymer import Anonymer

GENOME_FILE = "/data/genome/human/hg38/hg38.fa"
GENOME_INDEX = "/data/genome/human/hg38/bowtie2_index/hg38"
PROCESS_COUNT = 8

def main():
    # print start time
    start_time = datetime.now()
    print('''Anonymizer - genome data anonymization tool \nAnonymizer Starting: {start: %Y-%m-%d %H:%M:%S}'''.format(start=start_time))

    # read arguments
    args = parse_args()
    if args.verbose:
        print(args)

    # run anonymizer
    anonym = Anonymer(r1=args.r1,
                      r2=args.r2,
                      out_dir=args.out_dir,
                      gen_file=args.genome,
                      gen_idx=args.bowtie_index,
                      threads=PROCESS_COUNT,
                      keep_temp=args.keep_temp,
                      verbose=args.verbose)
    anonym.substitute()

    # print the time of the end:
    end_time = datetime.now()
    print('Anonymizer Stopping: {finish:%Y-%m-%d %H:%M:%S}'.format(finish=end_time))
    print('Total time of run  : {duration}'.format(duration=end_time - start_time))

def parse_args():
    """
    Parse command arguments.
    :return: argparse arguments
    """
    parser = argparse.ArgumentParser()
    
    required = parser.add_argument_group("Required")
    required.add_argument('-1', '--r1', type=is_fastq_file, help="read 1 fastq file", required=True)
    required.add_argument('-2', '--r2', type=is_fastq_file, help="read 2 fastq file", required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-i', '--fai', type=is_file, help="fai file")
    optional.add_argument('-o', '--out-dir', type=is_dir, help="output directory", default=os.getcwd())
    optional.add_argument('-p', '--processes', type=is_pos_int, help="number of working processes", default=PROCESS_COUNT)
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    optional.add_argument('-k', '--keep-temp', action='store_true', help="do not delete temporary files")
    optional.add_argument('--genome', type=is_file, help="reference genome, default=%s" % GENOME_FILE, default=GENOME_FILE)
    optional.add_argument('--bowtie-index', type=str, help="bowtie index, default=%s" % GENOME_INDEX, default=GENOME_INDEX)
    return parser.parse_args()


def is_dir(value):
    if os.path.isdir(value):
        return os.path.abspath(value)
    else:
        try:
            os.makedirs(value)
        except OSError:
            raise argparse.ArgumentTypeError("Value %s is not a directory and cannot be created." % value)
        return value


def is_file(value):
    if os.path.isfile(value):
        return os.path.abspath(value)
    else:
        raise argparse.ArgumentTypeError("Value %s is not a file." % value)


def is_pos_int(value):
    pos_int_pattern = re.compile("^[1-9]\d*$")
    if pos_int_pattern.match(value):
        return int(value)
    else:
        raise argparse.ArgumentTypeError("Value %s is not a non-zero positive integer." % value)


def is_fastq_file(value):
    fastq_pattern = re.compile(".*\.(fastq|fastq\.gz)?$")
    if os.path.isfile(value) and fastq_pattern.match(value):
        return os.path.realpath(value)
    else:
        raise argparse.ArgumentTypeError("Value %s is not a fastq file." % value)


if __name__ == "__main__":
    status = main()
    sys.exit(status)
