import os

RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')


# def main():
#     args = parse_args()
#     anonym = Anonymer(r1=args.r1,
#                       r2=args.r2,
#                       out_dir=args.out_dir,
#                       gen_file=GENOME_FILE,
#                       gen_idx=GENOME_INDEX,
#                       threads=THREAD_COUNT,
#                       keep_temp=args.keep_temp,
#                       verbose=args.verbose)
#     anonym.substitute()
#
#
# def parse_args():
#     parser = argparse.ArgumentParser()
#
#     required = parser.add_argument_group("Required")
#     required.add_argument('-1', '--r1', type=is_fastq_file, help="read 1 fastq file", required=True)
#     required.add_argument('-2', '--r2', type=is_fastq_file, help="read 2 fastq file", required=True)
#
#     optional = parser.add_argument_group("Optional")
#     optional.add_argument('-i', '--fai', type=is_file, help="fai file")
#     optional.add_argument('-o', '--out-dir', type=is_dir, help="output directory", default=os.getcwd())
#     optional.add_argument('-t', '--threads', type=is_pos_int, help="number of working threads", default=THREAD_COUNT)
#     optional.add_argument('-v', '--verbose', action='store_true', help="explain what is being done")
#     optional.add_argument('-k', '--keep-temp', action='store_true', help="do not delete temporary files")
#     return parser.parse_args()

class Varlock:
    AES_KEY_LENGTH = 32
    
    def __init__(self, dir):
        self.dir = dir
    
    def add_user(self, username, password):
        """
        Create asymetric key pair. Private key is encrypted by provided password.
        :return:
        """
        pass
    
    def remove_user(self, username):
        """
        Remove asymetric key pair.
        :return:
        """
        pass
    
    def update_password(self, username, old_password, new_password):
        """
        Decrypt asymetric key pair by old password and encrypt it with new password.
        :param username:
        :param old_password:
        :param new_password:
        :return:
        """
        pass
    
    def dasdas(self):
        pass
    
    @classmethod
    def gen_password(cls):
        return os.urandom(cls.AES_KEY_LENGTH)


if __name__ == "__main__":
    # Varlock.main()
    pass
