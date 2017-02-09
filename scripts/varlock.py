import io
import os

import varlock


class Varlock:
    RSA_KEY_LENGTH = 4096
    AES_KEY_LENGTH = 32
    
    @staticmethod
    def create_vac(bam_file, vcf_file, out_vac_file):
        """
        :param bam_file:
        :param vcf_file:
        :param out_vac_file:
        :return:
        """
        vac = varlock.Vac(bam_file)
        vac.vcf2vac(vcf_file, out_vac_file)
    
    @classmethod
    def encrypt(
            cls,
            bam_file,
            vac_file,
            public_key,
            mut_bam_file,
            diff_enc_file,
            verbose=False
    ):
        """
        Mutate BAM file and store it along with DIFF file encrypted by AES key which is encrypted by RSA public key.
        Output formats:
        .mut.bam
        .diff.enc
        .aes.enc
        :param bam_file: pysam.AlignmentFile
        :param vac_file: VAC file
        :param public_key: An RSA key object (`RsaKey`) containing the public key
        :param mut_bam_file: output pysam.AlignmentFile
        :param diff_enc_file: output binary file
        :param verbose:
        :return: AES key (to diff file) encrypted with public key
        """
        aes_key = os.urandom(cls.AES_KEY_LENGTH)
        out_diff = io.BytesIO()
        
        mut = varlock.Mutator(verbose=verbose)
        mut.mutate(
            in_vac_file=vac_file,
            in_bam_file=bam_file,
            out_bam_file=mut_bam_file,
            out_diff_file=out_diff
        )
        
        # encrypt diff with new AES key
        aes = varlock.FileAES(aes_key)
        aes.encrypt(out_diff, diff_enc_file)
        
        # encrypt AES key with public key
        return public_key.encrypt(aes_key)
    
    def decrypt(
            self,
            mut_bam_file,
            private_key,
            aes_enc_key,
            diff_enc_file,
            bam_file,
            ref_name=None,
            start_pos=None,
            end_pos=None,
            verbose=False
    ):
        """
        Reverse operation for mutate. Restore original bam file.
        Output formats:
        .bam
        :param mut_bam_file: input BAM
        :param private_key: An RSA key object (`RsaKey`) containing the private key
        :param aes_enc_key:
        :param diff_enc_file:
        :param bam_file: output BAM
        # :param ref_name:
        # :param start_pos:
        # :param end_pos:
        :param verbose:
        """
        if not self.has_privilege(
                private_key,
                aes_enc_key,
                diff_enc_file,
                mut_bam_file,
                # ref_name,
                # start_pos,
                # end_pos
        ):
            raise PermissionError("SHIT")
        
        # decrypt aes key
        aes_key = private_key.decrypt(aes_enc_key)
        
        # decrypt diff
        diff_file = io.BytesIO()
        aes = varlock.FileAES(aes_key)
        aes.decrypt(diff_enc_file, diff_file)
        
        # unmutate
        mut = varlock.Mutator(verbose)
        mut.unmutate(
            mut_bam_file=mut_bam_file,
            diff_file=diff_file,
            out_bam_file=bam_file,
            ref_name=ref_name,
            start_ref_pos=start_pos,
            end_ref_pos=end_pos,
        )
    
    @classmethod
    def add_privilege(
            cls,
            public_key,
            private_key,
            aes_enc_key,
            mut_bam_file,
            diff_enc_file,
            out_diff_enc_file,
            ref_name=None,
            start_pos=None,
            end_pos=None
    ):
        """
        Output formats:
        .diff.enc
        .aes.enc
        :param public_key:
        :param private_key:
        :param aes_enc_key:
        :param diff_enc_file:
        :param out_diff_enc_file:
        :param mut_bam_file:
        :param ref_name:
        :param start_pos:
        :param end_pos:
        :return: AES key (to diff file) encrypted with public key
        """
        if not cls.has_privilege(
                private_key,
                aes_enc_key,
                diff_enc_file,
                mut_bam_file,
                ref_name,
                start_pos,
                end_pos
        ):
            raise PermissionError("SHIT")
        
        # decrypt AES key
        aes_key = private_key.decrypt(aes_enc_key)
        aes = varlock.FileAES(aes_key)
        
        # decrypt diff
        diff_file = io.BytesIO()
        aes.decrypt(diff_enc_file, diff_file)
        
        fai_list = varlock.get_fai_list(mut_bam_file)
        fai_dict = varlock.fai_list2dict(fai_list)
        
        # slice diff
        out_diff_file = cls.__slice(diff_file, ref_name, start_pos, end_pos, fai_dict)
        
        # encrypt diff with aes key
        out_aes_key = os.urandom(cls.AES_KEY_LENGTH)
        out_aes = varlock.FileAES(out_aes_key)
        out_aes.encrypt(out_diff_file, out_diff_enc_file)
        
        # encrypt aes key
        return public_key.encrypt(out_aes_key)
    
    @classmethod
    def __start_index(cls, ref_name, start_pos, fai_dict):
        if start_pos is None:
            # access is given from the start of the reference
            start_index = fai_dict[ref_name].start
        else:
            start_index = varlock.pos2index(ref_name, start_pos, fai_dict)
        
        return start_index
    
    @classmethod
    def __end_index(cls, ref_name, end_pos, fai_dict):
        if end_pos is None:
            # access is given from the start of the reference
            end_index = fai_dict[ref_name].start + fai_dict[ref_name].length - 1
        else:
            end_index = varlock.pos2index(ref_name, end_pos, fai_dict)
        
        return end_index
    
    @classmethod
    def __slice(cls, diff_file, ref_name, start_pos, end_pos, fai_dict):
        if ref_name is not None:
            if ref_name not in fai_dict:
                raise ValueError("Reference name %s not found." % ref_name)
        
            return varlock.Diff.slice(
                diff_file=diff_file,
                start_index=cls.__start_index(ref_name, start_pos, fai_dict),
                end_index=cls.__end_index(ref_name, end_pos, fai_dict)
            )
    
    @classmethod
    def has_privilege(
            cls,
            private_key,
            aes_enc_key,
            diff_enc_file,
            mut_bam_file,
            ref_name=None,
            start_pos=None,
            end_pos=None
    ):
        """
        :param private_key:
        :param aes_enc_key:
        :param diff_enc_file:
        :param mut_bam_file:
        :param ref_name:
        :param start_pos:
        :param end_pos:
        :return:
        """
        if ref_name is not None and start_pos is not None and end_pos is not None:
            # range is specified
            if start_pos >= end_pos:
                # wrong range
                raise ValueError("End position must be greater than start position.")
        
        # TODO 1 try to decrypt aes
        aes_key = private_key.decrypt(aes_enc_key)
        
        # TODO 2 try to decrypt diff, only header part
        
        # TODO 3 check the range specified in diff's header
        
        # TODO 4 compare checksum of BAM to checksum from diff's header
        
        # TODO 5 check if range is present in BAM
        
        return True
