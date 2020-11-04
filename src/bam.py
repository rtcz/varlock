import pysam
from pysam.libcalignmentfile import VALID_HEADER_TYPES, KNOWN_HEADER_FIELDS, AlignmentFile, AlignmentHeader

MUT_TAG = 'mt'
MUT_BAM_TAG = 'bm'
MUT_VAC_TAG = 'vc'

# add mut tag to pysam's header format
VALID_HEADER_TYPES[MUT_TAG] = dict
KNOWN_HEADER_FIELDS[MUT_TAG] = {MUT_BAM_TAG: str, MUT_VAC_TAG: str}


def mut_header(header: AlignmentHeader, bam_checksum: str, vac_checksum: str) -> AlignmentHeader:
    if MUT_TAG in header:
        raise ValueError("File appears to be already mutated.")
    else:
        bm = MUT_BAM_TAG + ':' + bam_checksum
        vc = MUT_VAC_TAG + ':' + vac_checksum
        mut_line = '@' + MUT_TAG + '\t' + bm + '\t' + vc + '\n'
        text_header = str(header) + mut_line
        return AlignmentHeader.from_text(text_header)


def unmut_header(header: AlignmentHeader) -> AlignmentHeader:
    # drop empty line
    header_lines = str(header).split('\n')[:-1]
    last_line = header_lines[-1]
    segments = last_line.split('\t')

    if segments and '@' + MUT_TAG == segments[0]:
        # drop mut line from the header
        text_header = '\n'.join(header_lines[:-1]) + '\n'
        return AlignmentHeader.from_text(text_header)
    else:
        # TODO warning
        print('File does not appear to be mutated.')
        # return the original header
        return header


def open_bam(*args, **kwargs):
    # https://github.com/pysam-developers/pysam/issues/939
    pysam.set_verbosity(0)
    return AlignmentFile(*args, **kwargs)
