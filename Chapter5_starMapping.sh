#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 380G
#SBATCH -t 10-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis/scripts/batchLogs/starMapping_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis/scripts/batchLogs/starMapping_\%A_\%a.err # STDERR
#SBATCH --job-name=STARMap
#SBATCH --array=0-7

taskID=${SLURM_ARRAY_TASK_ID:-0}

# initiate log reporting
function log_line() {
	echo $(date) "$1" >&2
}



##==== Objective =====

#-- Map RNS Seq reads of Watkins assemblies to theier genomes and count number of reads mapping

## ==== load tools =============
#===========================================

source package 266730e5-6b24-4438-aecb-ab95f1940339   # STAR 2.7.10a

source package /tgac/software/testing/bin/subread-1.6.0_uv  # subread

source package 7b4623ee-93ca-4c66-b740-361a5e43c0e4   ## fastp 0.23.1

source package cf33bb3d-b7c0-450b-ac7d-ebedff25d497   ## BBMap 39.01 

source package c92263ec-95e5-43eb-a527-8f1496d56f1a   ## samtools 1.18

## =======  CONFIG ==========

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis"
cd ${workDir}

#--------------------------
#---- variables --------
#--------------------------

# reference and gff3
sampleFastaDir="/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies/2025-07-02_Watkins_assemblies"
helixergffs=/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/helixerOut

# two time points, 4 biological replicates each, hence eight samples
## assign the lumr array id to the queries so it reads all and submits jobs for all


declare -a queries=(\
"W42_T1_11_01" \
"W42_T1_11_02" \
"W42_T1_11_03" \
"W42_T1_11_04" \
"W42_T2_17_01" \
"W42_T2_17_02" \
"W42_T2_17_03" \
"W42_T2_17_04"  
)

#="${[$SLURM_ARRAY_TASK_ID]}"
#query="W42_T1_11_04"
query=${queries[$taskID]}

inputFasta=${sampleFastaDir}/IRN_WATDE0042_v1.fasta
genesGFF=${helixergffs}/WATDE0042_WG_helixer.gff3
trimdReadsDir=${workDir}/trimmedReads

starOut=${workDir}/starOut
queryOutIndex=${starOut}/genomeIndex # we have created this previously, hence we just copied it and used a central rather the re-create 
queryOutMapping=${starOut}/${query}

# make query specific out dir 
mkdir -p ${starOut} ${queryOutIndex} ${queryOutMapping}


##------- Monitor variables -------------
#------------------------------------------
log_line "current working dir: $PWD"
log_line "slum job ID: ${SLURM_ARRAY_TASK_ID}"
log_line "current sample: ${query}"
log_line "current genome file: ${inputFasta}"
log_line "current gff3: ${genesGFF}"
log_line "genome output directory: ${queryOutIndex}"
log_line "mapping output directory: ${queryOutMapping}"


#===== Check main input files ==================
#======================================================================================
if [[ ! -s "$inputFasta" ]]; then
    log_line "ERROR: genome fasta file for $query missing or empty" >&2
    exit 1
fi
if [[ ! -s "$genesGFF" ]]; then
    log_line "ERROR: Genes Annot gff3 file for $query missing or empty" >&2
    exit 1
fi

#========= Run STAR geneome index ========

#cp -r /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/starOut/WATDE0042/genomeIndex ${queryOutIndex}

#STAR --runThreadN 24 \
# --runMode genomeGenerate \
# --genomeDir ${queryOutIndex}\
# --genomeFastaFiles ${inputFasta}\
# --sjdbGTFfile ${genesGFF} \
# --sjdbOverhang 148-1 \
# --genomeSAindexNbases 13 \
# --limitGenomeGenerateRAM 100531998474

 # Check indexing is available by checking the SA file 

if [[ -s "${queryOutIndex}/SA" ]]; then
    log_line "STAR genome indexing for ${query} complete !!"
else
    log_line "STAR genome indexing for ${query} failed" >&2
    exit 1
fi

#log_line "Read correction and validation starting ..."

#========== Run STAR mapping reads to the indexed genome ==========


# Concatanate the trimmed files then gzip them, we will do this for sample 2 in T1 and sample 3 in T2

#cat ${trimdReadsDir}/${query}/*_1P.fq > ${trimdReadsDir}/${query}/${query}_trmdMgd_1P.fq

#cat ${trimdReadsDir}/${query}/*_2P.fq > ${trimdReadsDir}/${query}/${query}_trmdMgd_2P.fq

# for the majority of samples

R1=$(ls ${trimdReadsDir}/${query}/${query}_*_1P.fq)
R2=$(ls ${trimdReadsDir}/${query}/${query}_*_2P.fq)

# for the two samples we had to concantenate

#R1=${trimdReadsDir}/${query}/${query}_trmdMgd_1P.fq
#R2=${trimdReadsDir}/${query}/${query}_trmdMgd_2P.fq

log_line "R1 files: ${R1}"
log_line "R2 files: ${R2}"

ulimit -n 65535 ## increase the number of jobs that can be opened and worked on simulatanously so that STAR does not fail

#STAR \
#--runThreadN 32 \
#--genomeDir ${queryOutIndex} \
#--readFilesIn ${R1} ${R2} \
#--outFileNamePrefix ${queryOutMapping}/${query}_ \
#--outSAMtype BAM SortedByCoordinate \
#--outFilterMultimapNmax 10 \
#--outFilterMismatchNmax 10 \
#--alignIntronMax 50000 \
#--alignMatesGapMax 100000 \
#--alignSJoverhangMin 8 \
#--alignSJDBoverhangMin 1 \
#--outFilterScoreMinOverLread 0.66 \
#--outFilterMatchNminOverLread 0.66 \
#--twopassMode Basic \
#--quantMode TranscriptomeSAM GeneCounts \
#--limitBAMsortRAM 180000000000 \
#--outSAMstrandField intronMotif


outBam=${queryOutMapping}/${query}_Aligned.sortedByCoord.out.bam

# Check bam output
if [[ -s "${outBam}" ]]; then
    log_line "DONE: Trimmed reads successfully mapped to genome for ${query} !!"
else
    log_line "ERROR: STAR mapping of reads to genome for ${query} failed, no output" >&2
    exit 1
fi


##========== run feature counts to count number of reads mapped to each transcript ==========


#featureCounts -T 32 \
# -p \
# -a ${genesGFF} \
# -t gene \
# -g ID \
# -o ${queryOutMapping}/${query}_counts_ID.txt ${outBam}

countsFile=${queryOutMapping}/${query}_counts_ID.txt

## index the sorted BAM file 

samtools index --thread 32 -c ${outBam}

 # Check bam output
if [[ -s "${countsFile}" ]]; then
    log_line "✅ DONE: Counts for reads mapping to genes of ${query} genome complete !!"
else
    log_line "❌ ERROR: NO counts generated for reads mapping to genes of ${query} genome" >&2
    exit 1
fi


###======== subset chr5B bam file

chr5BBamFiles="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/2026-May_RNASeqAnalysis/starOut/chr5BBamFiles"

samtools view -b ${outBam} chr5B > ${chr5BBamFiles}/${query}_Aligned.sortedByCoord_chr5B.out.bam

samtools index --thread 32 -c ${chr5BBamFiles}/${query}_Aligned.sortedByCoord_chr5B.out.bam