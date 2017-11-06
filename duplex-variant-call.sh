#!/bin/bash
# DS bash script
# Version 2.1
# 

# Step 1: Setup variables for run:
clear

# Set up error checking
# Stop on any error
set -e
# Stop on an error inside a pipeline
set -o pipefail
# Throw an error on calling an unassigned variable
#set -u
#decide which python version to use
#. $HOME/pyEnv/py2/bin/activate

##########
# Config #
##########
export COLOR_RED="\e[31;40m"; export  COLOR_GREEN="\e[32;40m"; export  COLOR_YELLOW="\e[33;40m"; export  COLOR_BLUE="\e[34;40m"; export  COLOR_MAGENTA="\e[35;40m"; export  COLOR_CYAN="\e[36;40m"; export COLOR_RED_BOLD="\e[31;1m"; export COLOR_RED_LIGHT="\e[91m"; export COLOR_GREEN_BOLD="\e[32;1m"; export COLOR_GREEN_LIGHT="\e[92m"; export  COLOR_YELLOW_BOLD="\e[33;1m"; export  COLOR_BLUE_BOLD="\e[34;1m"; export  COLOR_MAGENTA_BOLD="\e[35;1m"; export  COLOR_CYAN_BOLD="\e[36;1m"; export COLOR_END="\e[0m"; 


export PATH=/home/sgeadmin/git/vygr_bioinformatics/:/home/sgeadmin/git/vygr_bioinformatics/utils/:$PATH

PIPELINE_DIRECTORY=/home/sgeadmin/git/vygr_bioinformatics/
###############
#functions  ###
###############
. $PIPELINE_DIRECTORY/utils/user_def.functions.sh
. $PIPELINE_DIRECTORY/utils/common.functions.sh

#########
# USAGE #
#########

usage () {
cat << EOF
    duplex-variant-call.sh
	a pipeline for Duplex Sequencing at Voyager; adapted from github: https://github.com/loeblab/Duplex-Sequencing
        -h usage && exit 0 ;;

        required:
        -i input.r1.fq[.gz]
        -r input.r2.fq[.gz]
        
        optional:
        -d Duplex sequencing code path
        -F prefix of output file
        -s insert size for pair-end reads
        -m minimal reads number to form a tag family
        -M maxmal reads number to form a tag family
        -o cutoff [default: 0.7 ]
        -N nCutoff [default: 1 ]
        -c cpu [8]
        -L read length [default: 101]
        -B barcode length [default: 12]
        -S spacer length [5]
        -f filter set [default:'os']
        -t read types [default: 'dpm']
        -p reaFile 
        -k print out progress every K number of reads
        -O output directory [default: current work dir]
EOF

}

PIPELINE_VERSION=1.0
#############################
# ARGS reading and checking #
############################
while getopts "hi:r:d:F:s:m:M:C:N:L:B:S:f:t:p:k:c:O:D" OPTION; do
    case $OPTION in
        h)  usage && exit 0 ;;
        i)  read1in=`readlink -f "${OPTARG}"` ;;
        r)  read2in=`readlink -f "${OPTARG}"` ;;
        d)  DSpath=$OPTARG ;;
        F)  filePrefix=$OPTARG ;;
        s)  iSize=$OPTARG ;;
        m)  minMem=$OPTARG ;;
        M)  maxMem=$OPTARG ;;
        C)  cutOff=$OPTARG ;;
        N)  nCutOff=$OPTARG ;;
        L)  readLength=$OPTARG ;;
        B)  barcodeLength=$OPTARG ;;
        S)  spacerLength=$OPTARG ;;
        f)	filtersSet=$OPTARG ;;
        t)	readTypes=$OPTARG ;;
        p)	repFilt=$OPTARG ;;
        k)  readOut=$OPTARG ;;
        c)  CPU=$OPTARG ;;
        O)  OUTDIR=$OPTARG;;
        D)  CLEAN=1;;
        *)  usage && exit 1 ;;
    esac
done


# if INPUT_FASTQ is undefined, print out usage and exit
[[ -z "${read1in}" ]] && usage && echo2 "Missing option -i for input R1 fastq file  or file does not exist" "error"
[[ -z "${read2in}" ]] && usage && echo2 "Missing option -r for input R2 fastq file or file does not exist" "error"


##check the existence of fastq files
[ ! -f "${read1in}" ] && echo2 "Cannot find input file ${read1in}" "error"
[ ! -f "${read2in}" ] && echo2 "Cannot find input file ${read2in}" "error"
file=`basename "${read1in}"` && export filePrefix=${file%%.*}

##default values for parameters
[ ! -z "$DSpath" ] || DSpath='/home/sgeadmin/git/Duplex-Sequencing/Nat_Protocols_Version'
[ ! -z "$alignRef" ] || alignRef='/home/sgeadmin/bowtieIndexes/mm10/mm10.fa'
[ ! -z "$filePrefix" ] || filePrefix='SRR1613972_bash_test'
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z "$iSize" ] || iSize=-1
[ ! -z "$minMem" ] || minMem=3
[ ! -z "$maxMem" ] || maxMem=1000
[ ! -z "$cutOff" ] || cutOff=0.7
[ ! -z "$nCutOff" ] || nCutOff=1

[ ! -z "$readLength" ] || readLength=101
[ ! -z "$barcodeLength" ] || barcodeLength=12
[ ! -z "$spacerLength" ] || spacerLength=5
#FINAL_READ_LENGTH
readLength=$((readLength-barcodeLength-spacerLength))


[ ! -z "$filtersSet" ] || filtersSet='os'
[ ! -z "$readTypes" ] || readTypes='dpm'
[ ! -z "$repFilt" ] || repFilt=9
[ ! -z "$readOut" ] || readOut=9
echo $OUTDIR
##check the output directory and write permission
[ ! -z "$OUTDIR" ] || OUTDIR=$PWD
[ "OUTDIR" != `readlink -f "${PWD}"` ] && (mkdir -p "${OUTDIR}" || echo2 "Cannot create directory ${OUTDIR}" "error")
echo $OUTDIR
cd ${OUTDIR} || (echo2 "Cannot access directory ${OUTDIR}... Exiting..." "error")
echo $PWD
touch .writing_permission && rm -rf .writing_permission || (echo2 "Cannot write to directory ${OUTDIR}... Exiting..." "error")



#LOG_FILE_NAME
logFile=${filePrefix}.log.txt

#Output folder
output_folder=${filePrefix}


# Load required software into path using the Environment Modules Project (http://modules.sourceforge.net)
#module load Python
#module load BWA
#module load SAMtools

# Print out options used to log file
touch $logFile
echo "read 1:" $read1in | tee -a ${logFile}
echo "read 2:" $read2in | tee -a ${logFile}
echo "Run identifier: " $filePrefix | tee -a ${logFile}
echo "workign directory: " $OUTDIR | tee -a ${logFile}
echo "Program path: " $DSpath | tee -a ${logFile}
echo "Reference genome: " $alignRef | tee -a ${logFile}
echo "Barcode length: " $barcodeLength | tee -a ${logFile}
echo "Spacer length: " $spacerLength | tee -a ${logFile}
echo "insert size: " $iSize| tee -a ${logFile}
echo "Post-tag_to_header read length: " $readLength | tee -a ${logFile}
echo "Repetitive tag filter length: " $repFilt | tee -a ${logFile}
echo "Minimum family size: " $minMem | tee -a ${logFile}
echo "Maximum family size: " $maxMem | tee -a ${logFile}
echo "Consensus cutoff: " $cutOff | tee -a ${logFile}
echo "Consensus N cutoff: " $nCutOff | tee -a ${logFile}
echo "Read types: " $readTypes | tee -a ${logFile}
echo "Filters: " $filtersSet | tee -a ${logFile}
echo "" | tee -a ${logFile}

#  Step 2: Run tag_to_header.py on imput files
STEP=1
JOBUID=`echo "$read1in" | md5sum | cut -d" " -f1`
echo "Starting Run" | tee -a ${logFile}
echo "tag_to_header starting"  | tee -a ${logFile}
date | tee -a ${logFile}
echo "" | tee -a ${logFile}
[ ! -f .${JOBUID}.status.${STEP}.tag2head ] && \
	python ${DSpath}/tag_to_header.py --infile1 $read1in --infile2 $read2in --outprefix ${filePrefix} --tagstats --spacerlen ${spacerLength} --taglen ${barcodeLength} && \
touch .${JOBUID}.status.${STEP}.tag2head
[ ! -f .${JOBUID}.status.${STEP}.tag2head ] && echo2 "run tag to head failed" "error"
STEP=$((STEP+1))

# Step 3 and 4: Align sequences and sort

echo "Aligning with BWA and Sorting aligned sequences" | tee -a ${logFile}
date | tee -a ${logFile}
if [[ $read1in =~ \.t?gz$ ]]
then
	xumi_seq1=${filePrefix}.seq1.smi.fq.gz
	xumi_seq2=${filePrefix}.seq2.smi.fq.gz
else
	xumi_seq1=${filePrefix}.seq1.smi.fq
	xumi_seq1=${filePrefix}.seq1.smi.fq
fi
[ ! -f .${JOBUID}.status.${STEP}.bwaAlignment ] && \
	bwa aln -t ${CPU} $alignRef $xumi_seq1 > ${filePrefix}.seq1.aln && \
	bwa aln -t ${CPU} $alignRef $xumi_seq2 > ${filePrefix}.seq2.aln && \
	bwa sampe -s $alignRef ${filePrefix}.seq1.aln ${filePrefix}.seq2.aln $xumi_seq1 $xumi_seq2 > ${filePrefix}.pe.sam && \
	samtools view -Sbu ${filePrefix}.pe.sam | samtools sort - -o ${filePrefix}.pe.sort && \
touch .${JOBUID}.status.${STEP}.bwaAlignment
[ ! -f .${JOBUID}.status.${STEP}.bwaAlignment ] && echo2 "run 1st bwa alignment failed" "error"
STEP=$((STEP+1))

# Step 5 and 6: Run Consensus Maker and Sort SSCSs
echo "Starting Consensus Maker" | tee -a ${logFile}
date | tee -a ${logFile}
[ ! -f .${JOBUID}.status.${STEP}.consensusMaker ] && \
	python ${DSpath}/ConsensusMaker.py --infile ${filePrefix}.pe.sort.bam --tag_file ${filePrefix}.pe.tagcounts --tag_stats ${filePrefix}.pe.tagstats  --outfile ${filePrefix}.sscs.bam --minmem $minMem --maxmem $maxMem --read_length $readLength --cut_off $cutOff --Ncut_off $nCutOff --read_type $readTypes --filt $filtersSet --isize $iSize && \
	samtools view -bu ${filePrefix}.sscs.bam | samtools sort - ${filePrefix}.sscs.sort && \
touch .${JOBUID}.status.${STEP}.consensusMaker ]
[ ! -f .${JOBUID}.status.${STEP}.consensusMaker ] && echo2 "consensus making failed" "error"
STEP=$((STEP+1))

# Step 7: Run Duplex Maker
echo "Starting Duplex Maker" | tee -a ${logFile}
date  | tee -a ${logFile}
[ ! -f .${JOBUID}.status.${STEP}.duplexMaker ] && \
	python ${DSpath}/DuplexMaker.py --infile ${filePrefix}.sscs.sort.bam --outfile ${filePrefix}.dcs.bam --Ncutoff $nCutOff --readlength $readLength && \
touch .${JOBUID}.status.${STEP}.duplexMaker ]
[ ! -f .${JOBUID}.status.${STEP}.duplexMaker ] && echo2 "duplex consensus making failed" "error"
STEP=$((STEP+1))

# Step 8, 9 and 10: Align DCSs and sort aligned DCSs; index bam files
echo "Aligning DCSs and sort and index" | tee -a ${logFile}
date | tee -a ${logFile}
[ ! -f .${JOBUID}.status.${STEP}.duplexbwaAlignment ] && \
	bwa aln -t ${CPU}  $alignRef ${filePrefix}.dcs.r1.fq > ${filePrefix}.dcs.r1.aln && \
	bwa aln -t ${CPU} $alignRef ${filePrefix}.dcs.r2.fq > ${filePrefix}.dcs.r2.aln && \
	bwa sampe -s $alignRef ${filePrefix}.dcs.r1.aln ${filePrefix}.dcs.r2.aln ${filePrefix}.dcs.r1.fq ${filePrefix}.dcs.r2.fq > ${filePrefix}.dcs.sam && \
	samtools view -Sbu ${filePrefix}.dcs.sam | samtools sort - -o ${filePrefix}.dcs.aln.sort && \
	samtools index ${filePrefix}.dcs.aln.sort.bam && \
touch .${JOBUID}.status.${STEP}.duplexbwaAlignment ]
[ ! -f .${JOBUID}.status.${STEP}.duplexbwaAlignment ] && echo2 "duplex consensus alignment failed" "error"
STEP=$((STEP+1))

# Step 11: Clean up
echo "Finishing with run.. " $filePrefix | tee -a ${logFile}
echo "Cleaning.." | tee -a ${logFile}
date | tee -a ${logFile} 
python ${DSpath}/clean.py --scripts_folder $(pwd) --output_folder ${output_folder} 
