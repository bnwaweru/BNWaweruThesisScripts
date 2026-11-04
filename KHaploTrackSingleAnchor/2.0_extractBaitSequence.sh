#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH -t 5-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmcExtrBait_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmcExtrBait_\%A_\%a.err # STDERR
#SBATCH --job-name=kmcExtrBait


#=================================================================================================
# Objective: extract the reference sequence spanning the Rht-B1 on chr4B interval,chr4B:30860000-30910000
# generate its 31-mer database, and remove repetitive whole-genome k-mers.
#=================================================================================================


# Timestamped logging to STDERR
function log_line() {
	echo $(date) "$1" >&2
}

## load tools 

#==== load KMC 3.1.1
source package 88c30fea-a263-44d2-a0b6-07681251a41d
source package c92263ec-95e5-43eb-a527-8f1496d56f1a   ## samtools 1.18

#=================================================================================================
# Define variables and paths
#=================================================================================================

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"
cd ${workDir}


## Reference genome sequence
# Chinese Spring is used as the source of the bait sequence.
genome="chinesespring"
WGS_seq="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation/${genome}.fa"


## Whole-genome repetitive k-mer database
repeats_dir="/jic/scratch/projects/watseq/kmer_agis/00_kmers/pangenomes"
kmer_db_MDR_repeats=${repeats_dir}/${genome}/${genome}_MDR_k31_ci2


## Genomic interval containing the target locus
region=("chr4B:30860000-30910000")

log_line "Working directory: ${workDir}"
log_line "Reference genome: ${genome}"
log_line "Reference FASTA: ${WGS_seq}"
log_line "Target region: ${region}"
log_line "Repeat k-mer database: ${kmer_db_MDR_repeats}"


#=================================================================================================
# Extract the target interval from the indexed reference genome
#=================================================================================================

mkdir -p ${workDir}/baitSequences

log_line "Extracting bait sequence"

# The FASTA must already have a samtools faidx index.
# samtools faidx ${WGS_seq}

samtools faidx ${WGS_seq} ${region} > ${workDir}/baitSequences/${genome}_$(echo $region | sed 's/:/_/g').fa

baitSeq=${workDir}/baitSequences/${genome}_$(echo $region | sed 's/:/_/g').fa

log_line "Bait sequence written to: ${baitSeq}"

#=================================================================================================
# Generate a 31-mer database from the extracted bait sequence
#=================================================================================================

# Store each region-specific k-mer database in a separate directory.
outDir=${workDir}/kmcOut_$(echo $region | sed 's/:/_/g')
mkdir -p ${outDir}
mkdir -p ${outDir}/temp/


# Report expanded variables before running KMC
log_line "Running from: $PWD"
log_line "Region being processed: $(echo $region | sed 's/:/_/g')"
log_line "Bait sequence: ${baitSeq}"
log_line "Output directory: ${outDir}"
log_line "Output database: ${outDir}/${genome}_$(echo $region | sed 's/:/_/g')"
log_line "KMC temporary directory: ${outDir}/temp"
log_line "KMC run started"


# -k31: use 31-mers
# -ci1: retain k-mers occurring at least once
# -t24: use 24 threads
# -m100: use up to 100 GB memory
# -fm: input is FASTA
# Canonical k-mers are used because the -b option is not specified.

time kmc \
-k31 \
-ci1 \
-t24 \
-m100 \
-fm ${baitSeq} ${outDir}/${genome}_$(echo $region | sed 's/:/_/g') ${outDir}/temp/

log_line "KMC database generation of bait sequence completed "


#=================================================================================================
# Dump bait k-mers to a text file
#=================================================================================================

kmc_db_bait=${outDir}/${genome}_$(echo $region | sed 's/:/_/g')

log_line "Dumping bait k-mer database to text"

kmc_tools transform \
	${kmc_db_bait} -ci1 \
	dump -s ${kmc_db_bait}.txt

log_line "KMC dump of bait kmer seqs completed "
log_line "Total bait k-mers: $(wc -l < ${kmc_db_bait}.txt)"


#=================================================================================================
# Remove k-mers present in the whole-genome repetitive k-mer database
#=================================================================================================

# kmers_subtract retains k-mers present in the bait database but absent from
# the repetitive whole-genome database.

log_line "Filtering repetitive k-mers from the bait database"

kmc_tools -t24 simple \
	${kmc_db_bait} -ci1 \
	${kmer_db_MDR_repeats} -ci2 \
	kmers_subtract ${kmc_db_bait}-target-filtered-no-repeats

log_line "Repeat filtering of repeat kmers from bait kmer set completed "


# Dump the remaining non-repetitive bait k-mers to text.
log_line "Dumping non-repetitive bait k-mers to text"

kmc_tools transform \
	${kmc_db_bait}-target-filtered-no-repeats -ci1 \
	dump -s ${kmc_db_bait}-target-filtered-no-repeats.txt

log_line "Total non-repetitive bait k-mers: $(wc -l < ${kmc_db_bait}-target-filtered-no-repeats.txt)"
log_line "Bait k-mer extraction and repeat filtering completed"