#!/usr/bin/env bash

set -u # treat unset variables as an error
set -e # exit immediatelly after return of non-zero status


# "constants"
SCRIPT_NAME="$0"
SCRIPT_DIR=`dirname $SCRIPT_NAME`
QSUB_THREAD_COUNT=8
READ_FORMAT="*.fastq"

source "$SCRIPT_DIR/varlock_functions.sh"

function usage {
	echo "Usage:"
	echo "$SCRIPT_NAME [options]* -d <WORK_DIR> [-s <SAMPLE>]"
	echo "  <WORK_DIR>  working directory, root of the project tree"
	echo "  <SAMPLE>    name of paired reads sample"
	echo "Options:"
	echo "  -h      print this usage message"
	echo "  -f      force option, overwrites all output"
	echo "  -l      log to file <WORKDIR>/log/varlock.log"
	exit 1
}

# initialise parameters
WORK_DIR=""
SAMPLE=""
IS_FORCED=false
LOG2FILE=false

# parse options
while getopts hd:s:fl opt
do
	case $opt in
		h) usage;;
		d) WORK_DIR=`readlink -f "$OPTARG"`;;
		s) SAMPLE=`get_true_basename "$OPTARG"`;;
		f) IS_FORCED=true;;
		l) LOG2FILE=true;;
		*) usage;;
	esac
done

# check mandatory parameters
if [[ -z "$WORK_DIR" ]]
then
    usage
fi

MAPPING_DIR="$WORK_DIR/mapping"
READS_ORIG_DIR="$WORK_DIR/reads/original"
READS_TRIM_DIR="$WORK_DIR/reads/trimmed"
LOG_DIR="$WORK_DIR/log"

# check source dir
if [ ! -d "$READS_ORIG_DIR" ]
then
	echo "Source directory $READS_ORIG_DIR does not exists, exiting."
	exit 1
fi

# create missing working directories if they not exist
mkdir -p "$MAPPING_DIR"
mkdir -p "$READS_TRIM_DIR"
mkdir -p "$LOG_DIR"

if [[ -z "$SAMPLE" ]]
then
	# undefined sample
	SAMPLE_COUNT=`ls "$READS_ORIG_DIR/"$READ_FORMAT | wc -l`
	read -p "No sample specified, run all $SAMPLE_COUNT samples under $WORK_DIR via qsub (y/n)?" choice
	case "$choice" in
		y|Y )
		# run all samples
		VISITED=()
		for FILE in "$READS_ORIG_DIR/"$READ_FORMAT
		do
			SAMPLE=`get_sample_name "$FILE"`
			if [[ `array_contains VISITED "$SAMPLE"` == false ]]
			then
				QSUB_LOG_FILE="$LOG_DIR/$SAMPLE/qsub.log"
				QSUB_ERR_FILE="$LOG_DIR/$SAMPLE/qsub.err"

				mkdir -p `dirname "$QSUB_LOG_FILE"`
				mkdir -p `dirname "$QSUB_ERR_FILE"`

				echo "$SCRIPT_DIR/varlock.sh -d "$WORK_DIR" -s "$SAMPLE" -l" \
					| qsub \
					-l thr=$QSUB_THREAD_COUNT \
					-cwd \
					-o $QSUB_LOG_FILE \
					-e $QSUB_ERR_FILE \
					-N _$SAMPLE \
					-p 100
					visited+=$SAMPLE
			fi
			VISITED+=("$SAMPLE")
		done
			;;
		*) exit 0 ;;
	esac
else
	if [ $LOG2FILE == true ]
	then
		exec 1>>"$LOG_DIR/$SAMPLE/varlock.log"
		exec 2>>"$LOG_DIR/$SAMPLE/varlock.err"
	fi

	# remove old logs
	rm "$LOG_DIR/$SAMPLE/*"

	# run specific sample
	run_sample "$SAMPLE" \
		"$READS_ORIG_DIR" \
		"$READS_TRIM_DIR" \
		"$MAPPING_DIR" \
		"$LOG_DIR" \
		$IS_FORCED
fi

