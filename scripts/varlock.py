import argparse
import os
import sys

from Crypto.PublicKey import RSA

from varlocker import Varlocker


def main():
    try:
        locker = Varlocker()
        command, args = parse_command()
        if command == 'encrypt':
            # python3 varlock.py encrypt --pub_key resources/jozko.pub --bam examples/resources/sample.bam --vac examples/resources/sample.vac --out_bam resources/out.mut.bam --out_diff resources/out.diff.enc
            parsed_args = parse_encrypt_args(args)
            with open(parsed_args.pub_key, 'rb') as pub_key_file:
                rsa_key = RSA.importKey(pub_key_file.read())
                locker.encrypt(
                    bam_filename=parsed_args.bam,
                    vac_filename=parsed_args.vac,
                    rsa_key=rsa_key,
                    out_bam_filename=parsed_args.out_bam,
                    out_enc_diff_filename=parsed_args.out_diff,
                    verbose=parsed_args.verbose
                )
        
        elif command == 'decrypt':
            parsed_args = parse_decrypt_args(args)
        elif command == 'reencrypt':
            parsed_args = parse_reencrypt_args(args)
        elif command == 'vac':
            parsed_args = parse_vac_args(args)
        else:
            print("unrecognized command '%s'" % command)
    except InvalidCommandError:
        print_usage()


def parse_command():
    if len(sys.argv) < 2:
        raise InvalidCommandError("Command is missing.")
    else:
        command = sys.argv[1]
        args = sys.argv[2:]
        return command, args


def print_usage():
    print('Usage:\t\tvarlock <command> [options]')
    print()
    print('Command:\tencrypt\t\tcreate mutated BAM along with encrypted DIFF and it\'s key')
    print('\t\tdecrypt\t\trevert mutated BAM to original')
    print('\t\treencrypt\t\tcreate encrypted slice of DIFF')
    print('\t\tvac\t\tconvert VCF to VAC')


def parse_encrypt_args(args):
    parser = argparse.ArgumentParser(prog='varlock encrypt')
    
    required = parser.add_argument_group("Required")
    required.add_argument('-k', '--pub_key', type=is_file, help='public key', required=True)
    required.add_argument('-b', '--bam', type=is_file, help='BAM file', required=True)
    required.add_argument('-c', '--vac', type=is_file, help='VAC file', required=True)
    required.add_argument('-m', '--out_bam', type=str, help='output mutated BAM file', required=True)
    required.add_argument('-d', '--out_diff', type=str, help='output encrypted DIFF file', required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    return parser.parse_args(args)


def parse_decrypt_args(args):
    parser = argparse.ArgumentParser(prog='varlock decrypt')
    
    required = parser.add_argument_group("Required")
    required.add_argument('-k', '--key', type=is_file, help='private key', required=True)
    required.add_argument('-m', '--bam', type=is_file, help='mutated BAM file', required=True)
    required.add_argument('-d', '--diff', type=is_file, help='encrypted DIFF file', required=True)
    required.add_argument('-b', '--out_bam', type=str, help='output restored BAM file', required=True)
    
    optional = parser.add_argument_group("Optional")
    range_help = "range in one of following formats: 'chr1', 'chr1:chr2', 'chr1:10000:20000', 'chr1:10000:chr2:20000'"
    optional.add_argument('-r', '--range', type=is_sam_range, help=range_help)
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    return parser.parse_args(args)


def parse_reencrypt_args(args):
    parser = argparse.ArgumentParser(prog='varlock reencrypt')
    return parser.parse_args(args)


def parse_vac_args(args):
    parser = argparse.ArgumentParser(prog='varlock vac')
    return parser.parse_args(args)


def is_sam_range(value):
    """
    A genomic region, stated relative to a reference sequence.
    A region consists of reference name (‘chr1’), start (10000), and end (20000).
    Start and end can be omitted for regions spanning a whole chromosome.
    If end is missing, the region will span from start to the end of the chromosome.
    Within pysam, coordinates are 0-based, half-open intervals, i.e.,
    the position 10,000 is part of the interval, but 20,000 is not.
    An exception are samtools compatible region strings such as ‘chr1:10000:20000’,
    which are closed, i.e., both positions 10,000 and 20,000 are part of the interval.
    """
    range_args = value.split(':')
    # TODO
    if len(range_args) == 1:
        pass
    elif len(range_args) == 2:
        pass
    elif len(range_args) == 3:
        pass
    elif len(range_args) == 4:
        pass
    else:
        raise argparse.ArgumentTypeError("Value %s is not a sam range." % value)


def is_file(value):
    if os.path.isfile(value):
        return os.path.realpath(value)
    else:
        raise argparse.ArgumentTypeError("Value %s is not a file." % value)


class InvalidCommandError(Exception):
    pass


if __name__ == "__main__":
    main()
