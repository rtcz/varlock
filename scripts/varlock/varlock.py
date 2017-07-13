import argparse
import os
import sys

from Crypto.PublicKey import RSA

from varlock import Varlocker


def main():
    try:
        locker = Varlocker()
        command, args = parse_command()
        if command == 'encrypt':
            # TODO add option for ommiting optional fields
            # python3 varlock.py encrypt --pub_key resources/jozko.pub --bam resources/input.bam --vac resources/input.vac --out_bam resources/out.mut.bam --out_diff resources/out.diff.enc
            parsed_args = parse_encrypt_args(args)
            with open(parsed_args.pub_key, 'r') as pub_key_file:
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
            # python3 varlock.py decrypt --key resources/jozko --bam resources/out.mut.bam --diff resources/out.diff.enc --out_bam out.bam
            parsed_args = parse_decrypt_args(args)
            with open(parsed_args.key, 'r') as key_file:
                rsa_key = RSA.importKey(key_file.read(), passphrase=parsed_args.password)
                
                rsa_ver_key = None
                if parsed_args.ver_key is not None:
                    with open(parsed_args.ver_key, 'r') as ver_key_file:
                        rsa_ver_key = RSA.importKey(ver_key_file.read())
                
                locker.decrypt(
                    rsa_key=rsa_key,
                    bam_filename=parsed_args.bam,
                    enc_diff_filename=parsed_args.diff,
                    out_bam_filename=parsed_args.out_bam,
                    start_ref_name=parsed_args.range[0],
                    start_ref_pos=parsed_args.range[1],
                    end_ref_name=parsed_args.range[2],
                    end_ref_pos=parsed_args.range[3],
                    include_unmapped=parsed_args.include_unmapped,
                    unmapped_only=parsed_args.unmapped_only,
                    rsa_ver_key=rsa_ver_key,
                    verbose=parsed_args.verbose
                )
        elif command == 'reencrypt':
            # python3 varlock.py reencrypt -d resources/jozko -e resources/jozko.pub -b resources/out.mut.bam -s resources/input.diff -o resources/output.diff -v -p password
            parsed_args = parse_reencrypt_args(args)
            with open(parsed_args.key, 'r') as key_file, \
                    open(parsed_args.pub_key, 'r') as pub_key_file:
                rsa_key = RSA.importKey(key_file.read(), passphrase=parsed_args.password)
                rsa_enc_key = RSA.importKey(pub_key_file.read())
                
                rsa_ver_key = None
                if parsed_args.ver_key is not None:
                    with open(parsed_args.ver_key, 'r') as ver_key_file:
                        rsa_ver_key = RSA.importKey(ver_key_file.read())
                
                locker.reencrypt(
                    rsa_key=rsa_key,
                    rsa_enc_key=rsa_enc_key,
                    bam_filename=parsed_args.bam,
                    enc_diff_filename=parsed_args.diff,
                    out_enc_diff_filename=parsed_args.out_diff,
                    start_ref_name=parsed_args.range[0],
                    start_ref_pos=parsed_args.range[1],
                    end_ref_name=parsed_args.range[2],
                    end_ref_pos=parsed_args.range[3],
                    include_unmapped=parsed_args.include_unmapped,
                    unmapped_only=parsed_args.unmapped_only,
                    rsa_ver_key=rsa_ver_key,
                    verbose=True
                )
        
        elif command == 'vac':
            # python3 varlock.py vac --bam examples/resources/sample.bam --vcf examples/resources/sample.vcf.gz --vac examples/resources/sample.vac
            parsed_args = parse_vac_args(args)
            locker.create_vac(
                bam_filename=parsed_args.bam,
                vcf_filename=parsed_args.vcf,
                out_vac_filename=parsed_args.vac,
                verbose=parsed_args.verbose
            )
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
    required.add_argument('-k', '--pub_key', type=is_file, help='public key for encryption', required=True)
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
    required.add_argument('-k', '--key', type=is_file, help='private key for decryption', required=True)
    required.add_argument('-m', '--bam', type=is_file, help='mutated BAM file', required=True)
    required.add_argument('-d', '--diff', type=is_file, help='encrypted DIFF file', required=True)
    required.add_argument('-b', '--out_bam', type=str, help='output restored BAM file', required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-p', '--password', type=str, help='private key password')
    region_formats = "'chr1', 'chr1:chr2', 'chr1:10000:20000', 'chr1:10000:chr2:20000'"
    range_help = "range in one of following formats: " + region_formats + "; positions are 1-based"
    optional.add_argument('-r', '--range', type=is_sam_range, help=range_help, default=(None, None, None, None))
    optional.add_argument('-i', '--include_unmapped', action='store_true', help="include all unplaced unmapped reads")
    optional.add_argument('-u', '--unmapped_only', action='store_true', help="only unmapped reads")
    optional.add_argument('-s', '--ver_key', help="public key for verification", default=None)
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    return parser.parse_args(args)


def parse_reencrypt_args(args):
    parser = argparse.ArgumentParser(prog='varlock reencrypt')
    required = parser.add_argument_group("Required")
    required.add_argument('-d', '--key', type=is_file, help='private key for decryption', required=True)
    required.add_argument('-e', '--pub_key', type=is_file, help='public key for encryption', required=True)
    required.add_argument('-b', '--bam', type=is_file, help='mutated BAM file', required=True)
    required.add_argument('-s', '--diff', type=is_file, help='source DIFF', required=True)
    required.add_argument('-o', '--out_diff', type=str, help='output DIFF', required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-p', '--password', type=str, help='private key password')
    region_formats = "'chr1', 'chr1:chr2', 'chr1:10000:20000', 'chr1:10000:chr2:20000'"
    range_help = "range in one of following formats: " + region_formats + "; positions are 1-based"
    optional.add_argument('-r', '--range', type=is_sam_range, help=range_help, default=(None, None, None, None))
    optional.add_argument('-i', '--include_unmapped', action='store_true', help="include all unplaced unmapped reads")
    optional.add_argument('-u', '--unmapped_only', action='store_true', help="only unmapped reads")
    optional.add_argument('-s', '--ver_key', help="public key for verification", default=None)
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    
    return parser.parse_args(args)


def parse_vac_args(args):
    parser = argparse.ArgumentParser(prog='varlock vac')
    
    required = parser.add_argument_group("Required")
    required.add_argument('-b', '--bam', type=is_file, help='BAM file', required=True)
    required.add_argument('-f', '--vcf', type=is_file, help='VCF file', required=True)
    required.add_argument('-c', '--vac', type=is_file, help='output VAC file', required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    
    return parser.parse_args(args)


def parse_sam_range(value: str):
    """
    1-based positions are converted to 0-based
    :param value: range string
    :return: range tuple
    """
    start_ref_pos = None
    end_ref_pos = None
    
    range_args = value.split(':')
    if len(range_args) == 1:
        # chr1
        start_ref_name = range_args[0]
        end_ref_name = range_args[0]
    elif len(range_args) == 2:
        # chr1:chr2
        start_ref_name = range_args[0]
        end_ref_name = range_args[1]
    elif len(range_args) == 3:
        # chr1:10000:20000
        start_ref_name = range_args[0]
        start_ref_pos = int(range_args[1]) - 1
        end_ref_name = range_args[0]
        end_ref_pos = int(range_args[2]) - 1
    elif len(range_args) == 4:
        # chr1:10000:chr2:20000
        start_ref_name = range_args[0]
        start_ref_pos = int(range_args[1]) - 1
        end_ref_name = range_args[2]
        end_ref_pos = int(range_args[3]) - 1
    else:
        raise ValueError("Invalid range format")
    
    return start_ref_name, start_ref_pos, end_ref_name, end_ref_pos


def is_sam_range(value):
    try:
        return parse_sam_range(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Value %s has not required format." % value)


def is_file(value):
    if os.path.isfile(value):
        return os.path.realpath(value)
    else:
        raise argparse.ArgumentTypeError("Value %s is not a file." % value)


class InvalidCommandError(Exception):
    pass


if __name__ == "__main__":
    main()
