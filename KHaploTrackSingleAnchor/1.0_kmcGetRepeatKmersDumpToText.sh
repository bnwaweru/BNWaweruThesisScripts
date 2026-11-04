#!/bin/bash
#SBATCH --partition=jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 100G
#SBATCH -t 5-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmcWGSrepeats_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmcWGSrepeats_\%A_\%a.err # STDERR
#SBATCH --job-name=kmc
#SBATCH --array=0


# Timestamped logging to STDERR
function log_line() {
	echo $(date) "$1" >&2
}

#=== declare the array iterator
taskID=${SLURM_ARRAY_TASK_ID:-0}
SLURM_CPUS_PER_TASK=32

log_line "SLURM array task ID: ${taskID}"
log_line "CPUs assigned: ${SLURM_CPUS_PER_TASK}"


#==== load KMC 3.1.1
source package 88c30fea-a263-44d2-a0b6-07681251a41d

#==== Objective  ==========================
# Generate whole-genome 31-mer databases for repeat filtering in IBSpy analyses.
# Repeated k-mers are defined here as k-mers occurring at least twice.


# move into the working directory
workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"
cd ${workDir}

# === Directory with the genome fasta's
#GenomeDir="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation"
pgFastaDir="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation"


#=================================================================================================
# Use whole-genome k-mer databases to identify repetitive k-mers.
# A k-mer size of 31 is used for consistency with the IBSpy analysis.
#=================================================================================================

## Genome reference names

declare -a references=( \
#"arinaLrFor" \
"chinesespring" \
#"jagger" \
#"julius" \
#"lancer" \
#"landmark" \
#"mace" \
#"norin61" \
#"symattis" \
#"stanley"
)

# Select genome using the SLURM array index
genome=${references[$taskID]}


# directory to store genome k-mer databases
outDirKmc=/jic/scratch/projects/watseq/kmer_agis/00_kmers/pangenomes/${genome}
tempDir=${outDirKmc}/temp/
mkdir -p ${outDirKmc} ${tempDir}

# input genome sequence location
WGS_seq=${pgFastaDir}/${genome}.fa

# prefix of the MDR kmer library (ci2)
kmc_db_WGS_MDR=${outDirKmc}/${genome}_MDR_k31_ci2

log_line "Working directory: ${workDir}"
log_line "Genome selected: ${genome}"
log_line "Input genome FASTA: ${WGS_seq}"
log_line "KMC output directory: ${outDirKmc}"
log_line "Repeat k-mer database prefix: ${kmc_db_WGS_MDR}"
log_line "Starting KMC whole-genome k-mer counting"

#=======================================================================================
## Generate a database containing k-mers occurring at least twice
#=======================================================================================
# -k31: k-mer length 31
# -ci2: minimum k-mer count of 2
# -t32: 32 threads
# -m100: 100 GB memory
# -fm: multi-FASTA input

log_line "Starting repeated k-mer database generation"

time kmc \
-k31 \
-ci2 \
-t32 \
-m100 \
-fm ${WGS_seq} ${kmc_db_WGS_MDR} ${tempDir}

log_line "Finished MDR 2 repeated k-mer database generation"


#=======================================================================================
## Generate a database containing all whole-genome k-mers
#=======================================================================================
# -ci1 retains k-mers occurring at least once.

kmc_db_WGS=${outDirKmc}/${genome}_k31_ci1

log_line "All-k-mer database prefix: ${kmc_db_WGS}"
log_line "Starting complete whole-genome k-mer database generation"

time kmc \
-k31 \
-ci1 \
-t32 \
-m100 \
-fm ${WGS_seq} ${kmc_db_WGS} ${tempDir}

log_line "Finished complete whole-genome k-mer database generation"

log_line "KMC processing completed for ${genome}"
