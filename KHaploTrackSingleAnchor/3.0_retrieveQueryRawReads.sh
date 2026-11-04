#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 50G
#SBATCH -t 5-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmcFilterReads_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmcFilterReads_\%A_\%a.err # STDERR
#SBATCH --job-name=kmcExtrctReads
#SBATCH --array=0-1


# Select query using the SLURM array index
taskID=${SLURM_ARRAY_TASK_ID:-0}
SLURM_CPUS_PER_TASK=32


# Timestamped logging to STDERR
function log_line() {
	echo $(date) "$1" >&2
}

log_line "Starting k-mer-based read extraction"
log_line "SLURM array task ID: ${taskID}"


#=================================================================================================
# Objective: use repeat-filtered bait k-mers to retrieve reads from the target query genotype
# interval and recover their corresponding read mates.
#=================================================================================================


#=================================================================================================
# Load required packages
#=================================================================================================

source package 88c30fea-a263-44d2-a0b6-07681251a41d   ## KMC 3.1.1
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6  ## SeqKit 0.12.0

log_line "KMC and SeqKit packages loaded"


#=================================================================================================
# Define variables and paths
#=================================================================================================

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"
cd ${workDir}

log_line "Working directory: ${workDir}"


# Reference genome used to generate the bait sequence
genome="chinesespring"

# Target Rht-B1 interval on CS v1.0
region=("chr4B:30860000-30910000")


#=================================================================================================
# Define region-dependent input and output paths
#=================================================================================================

kmcOut=${workDir}/kmcOut_$(echo $region | sed 's/:/_/g')
outDirReads=${kmcOut}/${genome}-filtered_reads
mkdir -p ${outDirReads}

# Repeat-filtered bait k-mer database
kmc_db_repeat_free=${kmcOut}/${genome}_$(echo $region | sed 's/:/_/g')-target-filtered-no-repeats

# Raw-read locations for pangenome references and Watkins accessions
pgRefsRawReads="/jic/scratch/projects/watseq/kmer_agis/raw_reads"
watkinsRawReads="/jic/archive/projects/watseq"


#=================================================================================================
# Accessions from which target-region reads will be extracted
#=================================================================================================

declare -a queries=(
    "jagger/jagger2" \
    "stanley/stanley3"
    )

query=${queries[$taskID]}

# Use the first path component as the accession output name
outFileName=$(echo ${query} | cut -d "/" -f 1)

log_line "Processing query: ${query}"
log_line "Target region: ${region}"


#=================================================================================================
# Determine whether the query is a pangenome reference or Watkins accession
#=================================================================================================

pgRefNames=("arinaLrFor" \
    "chinesespring" \
    "jagger" \
    "julius" \
    "lancer" \
    "landmark" \
	"norin61" \
    "mace" \
    "stanley" \
    "symattis" \
	"cadenza" \
	"claire" \
	"paragon" \
	"weebil"
	)


# Create a text file containing the absolute paths of all input FASTQ files
if [[ "$(declare -p pgRefNames)" =~ "${outFileName}" ]]; then
    log_line "${outFileName} is in the pangenome reference array"

	# Input reads for pangenome references
	log_line "Retrieving read paths from: ${pgRefsRawReads}"
	ls -d ${pgRefsRawReads}/${query}/*.fastq.gz > ${kmcOut}/${outFileName}_reads.txt
else
    log_line "${outFileName} is not in the pangenome reference array"

	# Input reads for Watkins accessions or other non-pangenome samples
	log_line "Retrieving read paths from: ${watkinsRawReads}"
	find "${watkinsRawReads}" -type f -name "*${outFileName}*.gz" -print > "${kmcOut}/${outFileName}_reads.txt"
fi

log_line "Number of input read files: $(wc -l < ${kmcOut}/${outFileName}_reads.txt)"


#=================================================================================================
# Report expanded variables before read filtering
#=================================================================================================

log_line "Running from: $PWD"
log_line "SLURM array task ID: ${taskID}"
log_line "Reads-list file: ${kmcOut}/${outFileName}_reads.txt"
log_line "Sample: ${query}; output name: ${outFileName}"
log_line "Repeat-free bait database: ${kmc_db_repeat_free}"
log_line "Extracted-read output directory: ${outDirReads}"


#=================================================================================================
# Filter reads using the repeat-free bait k-mer database
#=================================================================================================

# Retain reads containing at least 10 matching bait k-mers.
kmer_filt=10

log_line "Starting KMC read filtering with threshold: ${kmer_filt} k-mers"

kmc_tools filter \
${kmc_db_repeat_free} -ci1 \
@${kmcOut}/${outFileName}_reads.txt -ci${kmer_filt} \
${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered.fastq


log_line "KMC read filtering completed"
log_line "Filtered reads saved to: ${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered.fastq"
log_line "Number of initially retrieved reads: $(awk 'END {print NR/4}' ${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered.fastq)"


#=================================================================================================
# Extract the identifiers of reads retrieved by KMC
#=================================================================================================

awk 'NR%4==1 {print substr($1,2)}' \
${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered.fastq \
> ${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered_readIDs.txt

log_line "Read identifiers extracted"
log_line "Number of read identifiers: $(wc -l < ${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered_readIDs.txt)"


#=================================================================================================
# Retrieve both members of each read pair using SeqKit
#=================================================================================================

log_line "Extracting read mates using the retrieved read identifiers"

# Load the raw-read paths into an array
mapfile raw_reads_file_list < ${kmcOut}/${outFileName}_reads.txt

# File containing the IDs of reads identified by KMC
ids_file=${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered_readIDs.txt

# Combined output file containing the forward and reverse reads
RF_Reads_file=${outDirReads}/${outFileName}_${kmer_filt}kmers_non_repeat_target_filtered_sequences_RF.fastq


# Remove an existing output file to avoid appending duplicate reads
if [ -e "$RF_Reads_file" ]; then
    log_line "Existing F/R read file for ${query} found; deleting it"
    rm "$RF_Reads_file"
else
    log_line "No existing F/R read file found; continuing"
fi


# Search each raw-read file and append matching reads to the combined output
for read_file in "${raw_reads_file_list[@]}";
do
	log_line "Searching read file: ${read_file}"

	seqkit grep \
	--pattern-file ${ids_file} \
	${read_file} >> ${RF_Reads_file}
done

log_line "Read-mate extraction completed"
log_line "Combined F/R reads saved to: ${RF_Reads_file}"
log_line "Total reads in combined F/R file: $(awk 'END {print NR/4}' ${RF_Reads_file})"
log_line "K-mer-based read extraction completed for ${query}"
