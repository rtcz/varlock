#!/usr/bin/env bash

# reference: /data/projects/exome/scripts/process_sample.sh

# "constants"
THREAD_COUNT=8
MIN_MAPQ=2
GENOME_FILE="/data/genome/human/hg38/hg38.fa"
GENOME_INDEX="/data/genome/human/hg38/bowtie2_index/hg38"
R1_SUFFIX="_R1.fastq"
R2_SUFFIX="_R2.fastq"
R_SUFFIX_FORMAT="_R[1,2]*"

# get basename without extension
# example:
# "path/file.ext1.ext2" returns "file"
#
# string    FILE
function get_true_basename {
	local FILE=$1
	local FILENAME=`basename "$FILE"`;
	echo "${FILENAME%%.*}"
}

function get_sample_name {
	local FILE=$1
	local FILENAME=`basename "$FILE"`;
	echo "${FILENAME%$R_SUFFIX_FORMAT}"
}

# check if array contains value
#
# array
# string
function array_contains {
	local ARRAY="$1[@]"
	local VALUE="$2"
	local RESULT=false
	set +u
	for ELEMENT in "${!ARRAY}"
	do
		if [[ $ELEMENT == $VALUE ]]
		then
			echo true
			break
		fi
	done
	set -u
	echo $RESULT
}

# get log file, optionaly defined
#
# string file
function get_log_file {
	local LOG_FILE="/dev/stdout"
	if [[ ${1+defined} == defined ]]
	then
		mkdir -p `dirname "$1"`
		LOG_FILE="$1"
	fi
	echo "$LOG_FILE"
}

# get value
#
# any   VALUE   optional value
# any   DEFAULT default value
function get_optional {
	local VALUE=$2
	if [[ ${1+defined} == defined ]]
	then
		VALUE=$1
	fi
	echo $VALUE
}

# create fastqc quality report
#
# string    FASTQ_FILE  file to evaluate
# boolean   IS_FORCED   (optional) regenerates existing files
# string    LOG_FILE    (optional)
function generate_fastqc {
	local FASTQC_DIRNAME="stats"
	local FASTQC_FILE_SUFFIX="_fastqc.html"
	local IS_FORCED=`get_optional "$2" false`

	# local FASTQ_FILENAME=`basename $FASTQ_FILE`;
	local FASTQ_FILE="$1";
	local FASTQ_DIR=`dirname "$FASTQ_FILE"`;
	local FASTQ_BASENAME=`get_true_basename "$FASTQ_FILE"`

	local FASTQC_FILENAME="$FASTQ_BASENAME$FASTQC_FILE_SUFFIX"
	local FASTQC_DIR="$FASTQ_DIR/$FASTQC_DIRNAME";
	local FASTQC_FILE="$FASTQC_DIR/$FASTQC_FILENAME";

	mkdir -p "$FASTQC_DIR"

	local LOG_FILE=`get_log_file "$3"`
	local ERR_FILE=`get_log_file "$3"`

	if [[ -f "$FASTQC_FILE" && "$IS_FORCED" == false ]]
	then
		echo "-- file $FASTQC_FILE already exists, skipping"
	else
		echo "++ generating $FASTQC_FILE"
		fastqc \
			-t "$THREAD_COUNT" \
			-o "$FASTQC_DIR" \
			"$FASTQ_FILE" \
			1>>"$LOG_FILE" \
			2>>"$ERR_FILE"
	fi
}

# bam qualitimap
#
# string    BAM_FILE
# boolean   IS_FORCED (optional) regenerates existing files
# string    LOG_FILE (optional)
function generate_qualimap {
	local QUALIMAP_DIRNAME="stats"
	local QUALIMAP_FILENAME="qualimapReport.html"
	local IS_FORCED=`get_optional "$2" false`

	local BAM_FILE="$1"
	local BAM_DIR=`dirname $BAM_FILE`
	# local BAM_FILENAME=`basename $BAM_FILE` # ???
	local SAMPLE=`get_true_basename "$BAM_FILE"`

	local QUALIMAP_DIR="$BAM_DIR/$QUALIMAP_DIRNAME/$SAMPLE"
	local QUALIMAP_FILE="$QUALIMAP_DIR/$QUALIMAP_FILENAME"

	mkdir -p "$QUALIMAP_DIR"

	local LOG_FILE=`get_log_file "$3"`
	local ERR_FILE=`get_log_file "$3"`

	if [[ -f "$QUALIMAP_FILE" && "$IS_FORCED" == false ]]
	then
		echo "-- file $QUALIMAP_FILE already exists, skipping"
	else
		echo "++ generating $QUALIMAP_FILE"
		qualimap bamqc \
			-bam $BAM_FILE \
			-nt $THREAD_COUNT \
			-outdir $QUALIMAP_DIR \
			--outside-stats \
			1>>$LOG_FILE \
			2>>$ERR_FILE
	fi
}

# create trimmed fastq files from original files
#
# string    SAMPLE      sample name
# string    SOURCE_DIR  directory with fastq pair files
# string    OUT_DIR     output directory for trimmed fastq pair files
# boolean   IS_FORCED   regenerates existing files
# string    LOG_FILE    (optional)
function trim_fastq {
	local SAMPLE="$1";
	local SOURCE_DIR="$2";
	local OUT_DIR=`readlink -f "$3"`
	local IS_FORCED=`get_optional "$4" false`

	local ORIG_R1_FILE="$SOURCE_DIR/$SAMPLE$R1_SUFFIX"
	local ORIG_R2_FILE="$SOURCE_DIR/$SAMPLE$R2_SUFFIX"

	local TRIM_R1_FILE="$OUT_DIR/$SAMPLE$R1_SUFFIX"
	local TRIM_R2_FILE="$OUT_DIR/$SAMPLE$R2_SUFFIX"

	local LOG_FILE=`get_log_file "$5"`
	local ERR_FILE=`get_log_file "$5"`

	if [ ! -d "$SOURCE_DIR" ]
	then
		echo "Source directory $SOURCE_DIR does not exists, skipping sample trimming."
		return 1
	fi

	if [ ! -d "$OUT_DIR" ]
	then
		echo "Output directory $OUT_DIR does not exists, skipping sample trimming"
		return 1
	fi

	if [[ -f $TRIM_R1_FILE && -f $TRIM_R2_FILE && "$IS_FORCED" == false ]]
	then
		echo "-- trimmed reads already exist, skipping"
	else
		echo "++ trimming fastq pair $ORIG_R1_FILE, $ORIG_R2_FILE to $TRIM_R1_FILE, $TRIM_R2_FILE"
		java -jar /usr/local/tools/trimmomatic-0.32/trimmomatic-0.32.jar PE \
			-threads $THREAD_COUNT \
			-phred33 \
			$ORIG_R1_FILE \
			$ORIG_R2_FILE \
			$TRIM_R1_FILE \
			/dev/null \
			$TRIM_R2_FILE \
			/dev/null \
			SLIDINGWINDOW:5:25 \
			MINLEN:30 \
			1>>"$LOG_FILE" \
			2>>"$ERR_FILE"
			# MAXINFO:35:0.8 \
	fi
}

# mapping fastq to sam file
#
# string    SAMPLE      sample name
# string    SOURCE_DIR  directory with fastq pair files
# string    BAM_FILE    output bam file
# boolean   IS_FORCED   (optional)
# string    LOG_FILE    (optional)
function fastq2bam {
	local SAMPLE="$1"
	local SOURCE_DIR="$2"
	local BAM_FILE="$3"
	local IS_FORCED=`get_optional "$4" false`
	local ERR_FILE=`get_log_file "$5"`

	local FASTQ_R1_FILE="$SOURCE_DIR/$SAMPLE$R1_SUFFIX"
	local FASTQ_R2_FILE="$SOURCE_DIR/$SAMPLE$R2_SUFFIX"
	# local BAM_DIR=`dirname "$BAM_FILE"`

	local BAM_NAME=`get_true_basename "$BAM_FILE"`
	local BAM_DIR=`dirname "$BAM_FILE"`

	if [ ! -d "$SOURCE_DIR" ]
	then
		echo "Source directory $SOURCE_DIR does not exists, skipping sample mapping."
		return 1
	fi

	if [ ! -d "$BAM_DIR" ]
	then
		echo "Output directory $BAM_DIR does not exists, skipping sample mapping"
		return 1
	fi

	if [[ -f "$BAM_FILE" && "$IS_FORCED" == false ]]
	then
		echo "-- reads already mapped to $BAM_FILE, skipping"
	else
		echo "++ mapping files $FASTQ_R1_FILE, $FASTQ_R2_FILE to $BAM_FILE"
		bowtie2 \
			-x $GENOME_INDEX \
			-1 $FASTQ_R1_FILE -2 $FASTQ_R2_FILE \
			--very-sensitive \
			-p $THREAD_COUNT \
			2>>"$ERR_FILE" \
		| samtools view \
			-bS \
			-T $GENOME_FILE \
			-q $MIN_MAPQ \
			- \
			2>>"$ERR_FILE" \
		| bamtools filter \
			-isMapped true \
			-isPaired true \
			-isProperPair true \
			2>>"$ERR_FILE" \
		| samtools sort \
			-@ $THREAD_COUNT \
			- \
			"$BAM_DIR"/"$BAM_NAME" \
			2>>"$ERR_FILE"

		echo "++ indexing $BAM_FILE"
		samtools index $BAM_FILE
	fi
}

# run single sampe through pipeline
#
# string    SAMPLE
# string    ORIG_DIR
# string    TRIM_DIR
# string    MAPPIG_DIR
# string    LOG_DIR
# boolean   IS_FORCED
function run_sample {
	local SAMPLE="$1"
	local READS_ORIG_DIR="$2"
	local READS_TRIM_DIR="$3"
	local MAPPING_DIR="$4"
	local LOG_DIR="$5"
	local IS_FORCED="$6"

	local R1_SUFFIX="_R1.fastq"
	local R2_SUFFIX="_R2.fastq"
	local DATE_FORMAT="+%Y-%m-%d_%H:%M:%S"

	local R1_FILENAME="$SAMPLE$R1_SUFFIX"
	local R2_FILENAME="$SAMPLE$R2_SUFFIX"
	local BAM_FILE="$MAPPING_DIR/$SAMPLE.bam"

	# check source files
	if [ ! -f "$READS_ORIG_DIR/$R1_FILENAME" ]
	then
		echo "Source file $READS_ORIG_DIR/$R1_FILENAME does not exists, skipping sample."
		return 1
	elif [ ! -f "$READS_ORIG_DIR/$R2_FILENAME" ]
	then
		echo "Source file $READS_ORIG_DIR/$R2_FILENAME does not exists, skipping sample."
		return 1
	fi

	echo "SAMPLE $SAMPLE BEGIN"

	R1_BASENAME=`get_true_basename "$R1_FILENAME"`
	R2_BASENAME=`get_true_basename "$R2_FILENAME"`

	# evaluating original reads
	echo `date $DATE_FORMAT` "- ORIGINAL FASTQC"
	generate_fastqc "$READS_ORIG_DIR/$R1_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/orig_R1.fastqc.log"
	generate_fastqc "$READS_ORIG_DIR/$R2_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/orig_R2.fastqc.log"

	# trimming
	echo `date $DATE_FORMAT` "- TRIMMING"
	trim_fastq "$SAMPLE" \
		"$READS_ORIG_DIR" \
		"$READS_TRIM_DIR" \
		$IS_FORCED \
		"$LOG_DIR/$SAMPLE/trimming.log"

	# evaluating trimmed reads
	echo `date $DATE_FORMAT` "- TRIMMED FASTQC at "
	generate_fastqc "$READS_TRIM_DIR/$R1_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/trim_R1.fastqc.log"
	generate_fastqc "$READS_TRIM_DIR/$R2_FILENAME" $IS_FORCED "$LOG_DIR/$SAMPLE/trim_R2.fastqc.log"

	# mapping fastq files to sam files
	echo `date $DATE_FORMAT` "- MAPPING"
	fastq2bam "$SAMPLE" "$READS_TRIM_DIR" "$BAM_FILE" $IS_FORCED "$LOG_DIR/$SAMPLE/mapping.log"

	# evaluating bam
	echo `date $DATE_FORMAT` "- QUALIMAP"
	generate_qualimap "$BAM_FILE" $IS_FORCED "$LOG_DIR/$SAMPLE/qualimap.log"

	echo "SAMPLE $SAMPLE END"
}

