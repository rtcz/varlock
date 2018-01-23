import argparse
import os
import sys
from datetime import datetime

from Crypto.PublicKey import RSA

from varlock_src.varlocker import Varlocker


def main():
    # print start time
    start_time = datetime.now()
    print('VarLock = "Variant Locker" genome data anonymization tool \nVarLock Starting : {start: %Y-%m-%d %H:%M:%S}'.format(start=start_time))
    
    try:
        command, args = parse_command()
        parse_args_function = {'encrypt': parse_encrypt_args, 'decrypt': parse_decrypt_args,
                               'reencrypt': parse_reencrypt_args, 'vac': parse_vac_args}
        if command in parse_args_function.keys():
            parsed_args = parse_args_function[command](args)
        else:
            raise InvalidCommandError("unknown command %s (use: %s)" % (command, ' '.join(parse_args_function.keys())))
        locker = Varlocker(verbose=parsed_args.verbose)
        
        if command == 'encrypt':
            with open(parsed_args.key, 'r') as key_file, \
                    open(parsed_args.pub_key, 'r') as pub_key_file:
                rsa_key = RSA.importKey(key_file.read(), passphrase=parsed_args.password)
                rsa_pub_key = RSA.importKey(pub_key_file.read())
                locker.encrypt(
                    rsa_sign_key=rsa_key,
                    rsa_enc_key=rsa_pub_key,
                    bam_filename=parsed_args.bam,
                    vac_filename=parsed_args.vac,
                    out_bam_filename=parsed_args.out_bam,
                    out_enc_diff_filename=parsed_args.out_diff,
                    mut_p=parsed_args.mut_p,
                    seed=parsed_args.seed
                )
        
        elif command == 'decrypt':
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
                    rsa_ver_key=rsa_ver_key
                )
        elif command == 'reencrypt':
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
                    rsa_ver_key=rsa_ver_key
                )
        
        elif command == 'vac':
            parsed_args = parse_vac_args(args)
            locker.create_vac(
                bam_filename=parsed_args.bam,
                vcf_filename=parsed_args.vcf,
                out_vac_filename=parsed_args.vac,
                ref_fasta_filename=parsed_args.ref_fasta,
                skip_indels=parsed_args.no_indels
            )
        else:
            print("unrecognized command '%s'" % command)
    except InvalidCommandError:
        print_usage()
    
    # print the time of the end:
    end_time = datetime.now()
    print('VarLock Stopping : {finish:%Y-%m-%d %H:%M:%S}'.format(finish=end_time))
    print('Total time of run: {duration}'.format(duration=end_time - start_time))


def parse_command():
    if len(sys.argv) < 2:
        raise InvalidCommandError("command is missing")
    else:
        command = sys.argv[1]
        args = sys.argv[2:]
        return command, args


def print_usage():
    print('Usage:\t\tvarlock <command> [options]')
    print()
    print('Command:\tencrypt\t\tcreate mutated BAM along with encrypted DIFF and it\'s key')
    print('\t\tdecrypt\t\trevert mutated BAM to original')
    print('\t\treencrypt\tcreate encrypted slice of DIFF')
    print('\t\tvac\t\tconvert VCF to VAC')


def parse_encrypt_args(args):
    parser = argparse.ArgumentParser(prog='varlock encrypt')
    
    required = parser.add_argument_group("Required")
    required.add_argument('-k', '--key', type=is_file, help='private key for signing', required=True)
    required.add_argument('-e', '--pub_key', type=is_file, help='public key for encryption', required=True)
    required.add_argument('-b', '--bam', type=is_file, help='BAM file', required=True)
    required.add_argument('-c', '--vac', type=is_file, help='VAC file', required=True)
    required.add_argument('-m', '--out_bam', type=str, help='output mutated BAM file', required=True)
    required.add_argument('-d', '--out_diff', type=str, help='output encrypted DIFF file', required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-p', '--password', type=str, help='private key password')
    mut_p_help = 'random variant probability per genome base'
    optional.add_argument('-t', '--mut_p', type=is_mut_p, help=mut_p_help, default=0.0)
    optional.add_argument('-s', '--seed', type=str, help='random number generator seed')
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
    optional.add_argument('-f', '--ver_key', help="public key for verification", default=None)
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    return parser.parse_args(args)


def parse_reencrypt_args(args):
    parser = argparse.ArgumentParser(prog='varlock reencrypt')
    required = parser.add_argument_group("Required")
    required.add_argument('-k', '--key', type=is_file, help='private key for decryption', required=True)
    required.add_argument('-e', '--pub_key', type=is_file, help='public key for encryption', required=True)
    required.add_argument('-b', '--bam', type=is_file, help='mutated BAM file', required=True)
    required.add_argument('-c', '--diff', type=is_file, help='encrypted DIFF input', required=True)
    required.add_argument('-o', '--out_diff', type=str, help='encrypted DIFF output', required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-p', '--password', type=str, help='private key password')
    region_formats = "'chr1', 'chr1:chr2', 'chr1:10000:20000', 'chr1:10000:chr2:20000'"
    range_help = "range in one of following formats: " + region_formats + "; positions are 1-based"
    optional.add_argument('-r', '--range', type=is_sam_range, help=range_help, default=(None, None, None, None))
    optional.add_argument('-i', '--include_unmapped', action='store_true', help="include all unplaced unmapped reads")
    optional.add_argument('-u', '--unmapped_only', action='store_true', help="only unmapped reads")
    optional.add_argument('-f', '--ver_key', help="public key for verification", default=None)
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    
    return parser.parse_args(args)


def parse_vac_args(args):
    parser = argparse.ArgumentParser(prog='varlock vac')
    
    required = parser.add_argument_group("Required")
    required.add_argument('-b', '--bam', type=is_file, help='BAM file', required=True)
    required.add_argument('-f', '--vcf', type=is_file, help='VCF file', required=True)
    required.add_argument('-c', '--vac', type=str, help='output VAC file', required=True)
    
    optional = parser.add_argument_group("Optional")
    optional.add_argument('-r', '--ref-fasta', type=is_file, help='Reference FASTA file', default=None)
    optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
    optional.add_argument('--no-indels', action='store_true', help="Skip INDELs from VCF file.")
    
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


def is_mut_p(value):
    try:
        float_val = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Value %s is not a float." % value)
    
    if 0 <= float_val < 0.1:
        return float_val
    else:
        raise argparse.ArgumentTypeError("Value %s is not from interval <0, 0.1)." % value)


class InvalidCommandError(Exception):
    pass


if __name__ == "__main__":
    main()
