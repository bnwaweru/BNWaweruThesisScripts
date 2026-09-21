#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 380G
#SBATCH -t 10-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis/scripts/batchLogs/stringTie_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis/scripts/batchLogs/stringTie_\%A_\%a.err # STDERR
#SBATCH --job-name=stringTie
#SBATCH --array=0-7


# initiate log reporting
function log_line() {
	echo $(date) "$1" >&2
}


##==== Objective =====

#-- Count abundance of transcripts from mapped RNA Seq reads

## ==== load tools =============
#===========================================

source package 187d32ff-ff9f-4e1f-a711-3df9284d38b8 #  stringtie 2.1.1
source package 9117dc0d-5bb0-48ce-83c3-bc0db08e0cd3  # gffread 0.12.8

## =======  CONFIG ==========

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis"
cd ${workDir}

declare -a queries=( \
"W42_T1_11_01" \
"W42_T1_11_02" \
"W42_T1_11_03" \
"W42_T1_11_04" \
"W42_T2_17_01" \
"W42_T2_17_02" \
"W42_T2_17_03" \
"W42_T2_17_04"  
)

query="${queries[$SLURM_ARRAY_TASK_ID]}"

## === main =======
helixergffs=/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/helixerOut
starOut=${workDir}/starOut
stringTieOut=${workDir}/stringTieOut

#== query specific ==========
genesGFF=${helixergffs}/WATDE0042_WG_helixer.gff3
querySortedBam=${starOut}/${query}/${query}_Aligned.sortedByCoord.out.bam
queryOutStringtie=${stringTieOut}/${query}

#=== output files =============
outGtf=${stringTieOut}/WATDE0042_WG_helixer.gtf
ballGownDir=${stringTieOut}/${query}/ballGown
trnscrptsOut=${queryOutStringtie}/${query}_transcripts.gtf
abdnceFile=${queryOutStringtie}/${query}_WGAbundance.txt

# make query specific out dir 
mkdir -p ${stringTieOut} ${queryOutStringtie} ${ballGownDir}

##------- Monitor variables -------------
#------------------------------------------
log_line "current working dir: $PWD"
log_line "slum job ID: ${SLURM_ARRAY_TASK_ID}"
log_line "current genome sample: ${query}"
log_line "current gff3: ${genesGFF}"
log_line "query sorted bam: ${querySortedBam}"
log_line "stringtie out directory: ${queryOutStringtie}"
log_line "mapped bam sorted file: ${querySortedBam}"


#===== Check main input files ==================
#======================================================================================
if [[ -s "$querySortedBam" ]]; then
    log_line "Input sorted bam file for $query exists and is not empty"
    log_line "StringTie run starting ....."  >&2
else 
    log_line "ERROR: Input sorted bam file for $query missing or empty" >&2
    exit 1
fi



#++++++ Step 1: Convert gff3 to gtf accepted by stringtie

# -T → convert GFF3 → GTF
# -o → output GTF

#gffread ${genesGFF} -T -o ${outGtf}

# check gffread output 

#if [[ -s "${outGtf}" ]]; then
#    log_line "DONE: Gff3 succesfully converted to gft file for ${query} !!"
#else
#    log_line "ERROR: Gff3 convertion to gft file for ${query} failed, no output" >&2
#    exit 1
#fi

#++++++ Step 2: Run stringtie guided by existing transcripts, generating an abudance file and ballgown files

# -l name prefix for output transcripts 
# -A gene abundance estimation output file
#  -b enable output of Ballgown table files but these files will be created under the directory path given as <dir_path>

#time stringtie -p 32 \
#${querySortedBam} \
#-G ${outGtf} \
#-o ${trnscrptsOut} \
#-l "TraesWATDE0042" \
#-A ${abdnceFile} \
#-b ${ballGownDir}

abdnceFile=${queryOutStringtie}/${query}_WGAbundance.txt

 # Check abundance file output
if [[ -s "${abdnceFile}" ]]; then
    log_line "✅ DONE: Stringtie abundance and transcripts generation for ${query} complete !!"
else
    log_line "❌ ERROR: Stringtie abundance and transcripts generation for ${query} failed!!" >&2
    exit 1
fi

sort -k5 ${abdnceFile} |  awk '{if ($3 =="chr5B" && $5 >= 10000000 && $6 <= 13000000) print $0}' > ${queryOutStringtie}/${query}_WGAbundance_chr5B_10-13Mb.txt