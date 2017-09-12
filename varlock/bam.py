from pysam.calignmentfile import VALID_HEADER_TYPES, KNOWN_HEADER_FIELDS, AlignmentFile

MUT_TAG = 'mt'
MUT_BAM_TAG = 'bm'
MUT_VAC_TAG = 'vc'

# add mut tag to pysam's header format
VALID_HEADER_TYPES[MUT_TAG] = dict
KNOWN_HEADER_FIELDS[MUT_TAG] = {MUT_BAM_TAG: str, MUT_VAC_TAG: str}


def mut_header(header: dict, bam_checksum: str, vac_checksum: str):
    if MUT_TAG in header:
        raise ValueError("File appears to be already mutated.")
    else:
        result_header = header.copy()
        result_header[MUT_TAG] = {MUT_BAM_TAG: bam_checksum, MUT_VAC_TAG: vac_checksum}
        return result_header


def unmut_header(header_map):
    if MUT_TAG in header_map:
        del header_map[MUT_TAG]
        return header_map
    else:
        raise ValueError('File does not appear to be mutated.')


def open_bam(*args, **kwargs):
    return AlignmentFile(*args, **kwargs)
