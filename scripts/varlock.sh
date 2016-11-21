#!/usr/bin/env bash

set -u # treat unset variables as an error
set -e # exit immediatelly after return of non-zero status

# "constants"
SCRIPT_NAME="$0"
SCRIPT_DIR=`dirname $SCRIPT_NAME`
R1_SUFFIX="_R1.fastq"
R2_SUFFIX="_R2.fastq"
DATE_FORMAT="+%Y-%m-%d_%H:%M:%S"

source "$SCRIPT_DIR/varlock_f.sh"

function usage {
	echo "Usage:"
	echo "$SCRIPT_NAME [options]* -d <WORK_DIR> -s <SAMPLE>"
	echo "  <WORK_DIR>  working directory, root of the project tree"
	echo "  <SAMPLE>    name of paired reads sample"
	echo "Options:"
	echo "  -h      print this usage message"
	echo "  -f      force option, overwrites all output"
	echo "  -l      log to file <WORKDIR>/log/varlock.log"
	exit 1
}

# initialise variables
WORK_DIR=""
SAMPLE=""
IS_FORCED=false
LOG2FILE=false

# parse options
while getopts hd:s:f:l opt
do
	case $opt in
		h) usage;;
		d) WORK_DIR=`readlink -f ${OPTARG}`;;
		s) SAMPLE=`basename ${OPTARG}`;;
		f) IS_FORCED=true;;
		l) LOG2FILE=true;;
		*) usage;;
	esac
done

# check mandatory parameters
if [[ -z "$WORK_DIR" || -z "$SAMPLE" ]]
then
    usage
fi

MAPPING_DIR="$WORK_DIR/mapping"
READS_ORIG_DIR="$WORK_DIR/reads/original"
READS_TRIM_DIR="$WORK_DIR/reads/trimmed"
LOG_DIR="$WORK_DIR/log"

R1_FILENAME="$SAMPLE$R1_SUFFIX"
R2_FILENAME="$SAMPLE$R2_SUFFIX"
BAM_FILE="$MAPPING_DIR/$SAMPLE.bam"

# check sources
if [ ! -d "$READS_ORIG_DIR" ]
then
	echo "Source directory $READS_ORIG_DIR does not exists, exiting."
	exit 1
elif [ ! -f "$READS_ORIG_DIR/$R1_FILENAME" ]
then
	echo "Source file $READS_ORIG_DIR/$R1_FILENAME does not exists, exiting."
	exit 1
elif [ ! -f "$READS_ORIG_DIR/$R2_FILENAME" ]
then
	echo "Source file $READS_ORIG_DIR/$R2_FILENAME does not exists, exiting."
	exit 1
fi

# create missing working directories if they not exist
mkdir -p "$MAPPING_DIR"
mkdir -p "$READS_TRIM_DIR"
mkdir -p "$LOG_DIR"

# log configuration
if [ $LOG2FILE == true ]
then
	exec 1>>"$LOG_DIR/varlock.log"
	exec 2>>"$LOG_DIR/varlock.err"
fi

echo "SAMPLE $SAMPLE BEGIN"

R1_BASENAME=`get_true_basename "$R1_FILENAME"`
R2_BASENAME=`get_true_basename "$R2_FILENAME"`

# evaluating original reads
echo `date $DATE_FORMAT` "- ORIGINAL FASTQC"
run_fastqc "$READS_ORIG_DIR/$R1_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/orig_R1.fastqc.log"
run_fastqc "$READS_ORIG_DIR/$R2_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/orig_R2.fastqc.log"

# trimming
echo `date $DATE_FORMAT` "- TRIMMING"
trim_fastq "$SAMPLE" \
	"$READS_ORIG_DIR" \
	"$READS_TRIM_DIR" \
	$IS_FORCED \
	"$LOG_DIR/$SAMPLE/trimming.log"

# evaluating trimmed reads
echo `date $DATE_FORMAT` "- TRIMMED FASTQC at "
run_fastqc "$READS_TRIM_DIR/$R1_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/trim_R1.fastqc.log"
run_fastqc "$READS_TRIM_DIR/$R2_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/trim_R2.fastqc.log"

# mapping fastq files to sam files
echo `date $DATE_FORMAT` "- MAPPING"
fastq2bam "$SAMPLE" "$READS_ORIG_DIR" "$BAM_FILE" $IS_FORCED "$LOG_DIR/$SAMPLE/mapping.log"

# evaluating bam
echo `date $DATE_FORMAT` "- QUALIMAP"
run_qualimap "$BAM_FILE" $IS_FORCED "$LOG_DIR/$SAMPLE/qualimap.log"

echo "SAMPLE $SAMPLE END"


