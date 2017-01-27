import os

RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')

FAI_FILENAME = os.path.join(RESOURCES_DIR, "hg19.fa.fai")
FAI2_FILENAME = os.path.join(RESOURCES_DIR, "hs37d5.fa.fai")

TEXT_VAC_FILENAMES = [
    os.path.join(RESOURCES_DIR, "input_01.vac.txt"),
    os.path.join(RESOURCES_DIR, "input_02.vac.txt")
]
VAC_FILENAMES = [
    os.path.join(RESOURCES_DIR, "input_01.vac"),
    os.path.join(RESOURCES_DIR, "input_02.vac")
]

SAM_FILENAME = os.path.join(RESOURCES_DIR, "input.sam")
BAM_FILENAME = os.path.join(RESOURCES_DIR, "input.bam")

DESIRED_MUT_FILENAMES = [
    os.path.join(RESOURCES_DIR, "desired_mut_01.sam"),
    os.path.join(RESOURCES_DIR, "desired_mut_02.sam")
]
MUT_BAM_FILENAMES = [
    os.path.join(RESOURCES_DIR, "mut_01.bam"),
    os.path.join(RESOURCES_DIR, "mut_02.bam")
]
MUT_SAM_FILENAMES = [
    os.path.join(RESOURCES_DIR, "mut_01.sam"),
    os.path.join(RESOURCES_DIR, "mut_02.sam")
]

DIFF_FILENAMES = [
    os.path.join(RESOURCES_DIR, "mut_01.diff"),
    os.path.join(RESOURCES_DIR, "mut_02.diff")
]

TEXT_DIFF_FILENAMES = [
    os.path.join(RESOURCES_DIR, "mut_01.diff.txt"),
    os.path.join(RESOURCES_DIR, "mut_02.diff.txt")
]

DESIRED_DIFF_FILENAMES = [
    os.path.join(RESOURCES_DIR, "desired_mut_01.diff.txt"),
    os.path.join(RESOURCES_DIR, "desired_mut_02.diff.txt")
]

UNMUT_BAM_FILENAMES = [
    os.path.join(RESOURCES_DIR, "unmut_01.bam"),
    os.path.join(RESOURCES_DIR, "unmut_02.bam")
]

UNMUT_SAM_FILENAMES = [
    os.path.join(RESOURCES_DIR, "unmut_01.sam"),
    os.path.join(RESOURCES_DIR, "unmut_02.sam")
]

DESIRED_UNMUT_FILENAMES = [
    os.path.join(RESOURCES_DIR, "desired_unmut_01.sam"),
    os.path.join(RESOURCES_DIR, "desired_unmut_02.sam")
]
