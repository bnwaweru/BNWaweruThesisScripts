#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 100G
#SBATCH -t 10-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis/scripts/batchLogs/fastqcTrim_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis/scripts/batchLogs/fastqcTrim_\%A_\%a.err # STDERR
#SBATCH --job-name=fastqcTrim
#SBATCH --array=0-7

i=$SLURM_ARRAY_TASK_ID

# initiate log reporting
function log_line() {
	echo $(date) "$1" >&2
}


##==== Objective =====

#-- check quality and correct RNA Seq reads before mapping to genome with STAR

## ==== load tools =============
#===========================================

source package a138abc4-ed60-4774-8db8-80b4770b1710  # FastQC 0.11.9
source package a8a18f99-1c90-4175-8f58-330b0ad61cad  # MultiQC 1.19
source package 50fcf79b-73a3-4f94-9553-5ed917823423  # trimmomatic 0.39


## =======  CONFIG ==========

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis"
cd ${workDir}

readsDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2026-04-30-RNASeqData/origData/X204SC26043180-Z01-F002/01.RawData"

# two time points, 4 biological replicates each, hence eight samples

declare -a queries=("W42_T2_17_04" "W42_T2_17_03" "W42_T2_17_02" "W42_T2_17_01" "W42_T1_11_04" "W42_T1_11_03" "W42_T1_11_02" "W42_T1_11_01")
query=$(($i%8))
query=${queries[$query]}

rawReadsDir=${readsDir}/${query}

# we store all the fastc reports in one folder, so that we can then create one mltic report
fastQC=${workDir}/fastQC
#fastqcOut=${fastQC}/fastqcRaw
#multiqcOut=${fastQC}/multiqcRaw

# if we need to trim the reads they will go here, a folder per sample
trimdReadsDir=${workDir}/trimmedReads/${query}

mkdir -p ${fastQC} ${fastqcOut} ${multiqcOut} ${trimdReadsDir}

##------- Monitor variables -------------
#------------------------------------------
log_line "current working dir: $PWD"
log_line "slum job ID: ${i}"
log_line "current query file: ${query}"
log_line "reading raw reads from: ${rawReadsDir}"
log_line "output directory: ${fastqcOut}"
log_line "fastQC run starting ....."

#================================================================
#============== Run fastqc to check quality of raw reads =========
#================================================================

#for file in $(ls ${rawReadsDir}/*.fq.gz) ;\
#do log_line "now working on $file" ;\
#fastqc -t 24 ${file} -o ${fastqcOut}/;\
#done

log_line "fastQC run complete!"

#================================================================
#============== use multiqc to aggregate fastqc reports =========
#================================================================

#log_line "Running multiQC in file in folder ${fastqcOut}"

#multiqc ${fastqcOut}/ -o ${multiqcOut}

#log_line "Multiqc run complete !!"



#================================================================
#============== Run trimmomatic mainly to trim adaptaters =========
#================================================================

log_line "Trimmomatic run starting ....."

#for line in $(ls ${rawReadsDir}/*1.fq.gz) ;\
#do 
#outName=$(echo $line |cut -f 12 -d "/" | sed 's/_1.fq.gz//g') ;\
#log_line "Sample out file name ${outName}" ;\

#R2=${rawReadsDir}/${outName}_2.fq.gz
#log_line "Forward read ${line}" ;\
#log_line "Reverse read ${R2}" ;\

#trimmomatic PE -threads 32 -trimlog ${trimdReadsDir}/${query}_trimmomaticLog.txt \
# ${line} ${R2} -baseout ${trimdReadsDir}/${outName}_trmd.fq \
#ILLUMINACLIP:/adapters/TruSeq3-PE-2.fa:2:30:10 MINLEN:50;\
#cd ${trimdReadsDir};\
#rm *U.fq ;\
#rm ${trimdReadsDir}/${query}_trimmomaticLog.txt ;\
#done

log_line "Trimmomatic complete !!"

#================================================================
#============== Re-run fastqc to check quality of trimmed reads =========
#================================================================


fastqcOut=${fastQC}/fastqcTrmd
multiqcOut=${fastQC}/multiqcTrmd

mkdir -p ${fastqcOut} ${multiqcOut}

for file in $(ls ${trimdReadsDir}/*P.fq) ;\
do echo $file ;\
fastqc -t 32 ${file} -o ${fastqcOut}/;\
done

log_line "fastQC run complete, running multiqc ....."

#================================================================
#============== use multiqc to aggregate fastqc reports =========
#================================================================

multiqc ${fastqcOut}/ -o ${multiqcOut}

log_line "Multiqc run complete !!"

