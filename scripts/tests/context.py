import os

RESOURCES_DIR = os.path.join(os.path.dirname(__file__), 'resources')

FAI_FILENAME = os.path.join(RESOURCES_DIR, "hg19.fa.fai")
FAI2_FILENAME = os.path.join(RESOURCES_DIR, "hs37d5.fa.fai")

TEXT_VAC_FILENAMES = [
    os.path.join(RESOURCES_DIR, "input_01.vac.txt"),
    os.path.join(RESOURCES_DIR, "input_02.vac.txt")
]
IN_VAC_FILENAMES = [
    os.path.join(RESOURCES_DIR, "input_01.vac"),
    os.path.join(RESOURCES_DIR, "input_02.vac")
]

IN_SAM_FILENAME = os.path.join(RESOURCES_DIR, "input.sam")
IN_BAM_FILENAME = os.path.join(RESOURCES_DIR, "input.bam")

DESIRED_BAM_FILENAMES = [
    os.path.join(RESOURCES_DIR, "desired_01.sam"),
    os.path.join(RESOURCES_DIR, "desired_02.sam")
]
OUT_BAM_FILENAMES = [
    os.path.join(RESOURCES_DIR, "output_01.bam"),
    os.path.join(RESOURCES_DIR, "output_02.bam")
]
OUT_SAM_FILENAMES = [
    os.path.join(RESOURCES_DIR, "output_01.sam"),
    os.path.join(RESOURCES_DIR, "output_02.sam")
]

OUT_DIFF_FILENAMES = [
    os.path.join(RESOURCES_DIR, "output_01.diff"),
    os.path.join(RESOURCES_DIR, "output_02.diff")
]

DIFF_FILENAME = 'output.diff'
