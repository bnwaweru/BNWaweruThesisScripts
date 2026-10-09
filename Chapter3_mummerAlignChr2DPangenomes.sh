#!/bin/bash
#SBATCH --partition=RG-Cristobal-Uauy,nbi-medium,jic-medium,nbi-long,jic-long,RG-Cristobal-Uauy
#SBATCH --nodes=1
#SBATCH --cpus=2
#SBATCH --mem 240G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023-11_mummer_alignments_pangenomes/scripts/batch_logs/mummer_panrefs.%N.%j.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023-11_mummer_alignments_pangenomes/scripts/batch_logs/mummer_panrefs.%N.%j.err # STDERR
#SBATCH --job-name=mummer_2D
#SBATCH --array=1

# Align chromosome 2D sequences from wheat pangenome assemblies against
# Chinese Spring chromosome 2D using NUCmer (MUMmer4).
#
# The computational commands and workflow in this script are unchanged from
# the original version. Chromosome-specific names and paths have been updated
# from 2B to 2D, and explanatory comments have been added.

# S B A T C H --localscratch=ssd:200

# =============================================================================
# Objective
# =============================================================================

# Investigate the Limagrain chromosome 2D target interval corresponding to
# Chinese Spring RefSeq v1.0 chr2D:635-651 Mb.
#
# The synteny analysis recovered comparatively short intervals in some
# pangenome assemblies (approximately 3 Mb on average), whereas the Chinese
# Spring target interval is approximately 15 Mb. Whole-chromosome alignments
# were therefore used to examine whether deletions, inversions, translocations,
# other structural differences, or possible assembly-related differences may
# explain the short syntenic intervals.
#
# This interval is close to the distal region of Chinese Spring chromosome 2D,
# where complex or unexpected alignment patterns may occur.

# =============================================================================
# Load software
# =============================================================================

# Load MUMmer v4.0.0beta2.
source package eab121cb-2eb8-49c1-a9a5-a33754ea9fea

# =============================================================================
# Declare paths and variables
# =============================================================================

# Main working directory containing input sequences and output directories.
work_dir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023-11_mummer_alignments_pangenomes"

# Run all subsequent commands from the main working directory.
cd ${work_dir}

# Use the SLURM array-task identifier to select a pangenome assembly from the
# array declared below.

i=$SLURM_ARRAY_TASK_ID

# Pangenome assemblies to use as queries. 
declare -a references=(\
"arinaLrFor" \
"chinesespring" \
"jagger" \
"julius" \
"lancer" \
"landmark" \
"mace" \
"norin61" \
"stanley" \
"symattis" \
)

# Convert the current SLURM task identifier to the zero-based array index used
# by Bash. With one active query and --array=1, 1 % 1 evaluates to index 0.
# If additional queries are activated, both the #SBATCH --array range and the
# divisor in this expression must be updated to match the number of queries.
reference=$(($i%1))

# Directory containing the complete pangenome assembly FASTA files.
panref_genomes='/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation'

# Complete FASTA file for the pangenome assembly selected by the array task.
query=${panref_genomes}/${references[$reference]}.fa

# Directory containing chromosome 2D FASTA sequences extracted from each
# pangenome assembly.

pangenome_2D_seqs=${work_dir}/pangenome_2D_seqs

mkdir -p ${work_dir}/pangenome_2D_seqs

# Short name of the pangenome assembly selected for the current array task.
ref=${references[$reference]}

# output-directory commands are defined below.
ref_dir=${work_dir}/${ref}_mummer_out
mkdir -p ${work_dir}/${ref}_mummer_out

# =============================================================================
# Logging and run tracking
# =============================================================================

# Print timestamped messages to STDERR so that run information is captured in
# the SLURM error log.
function log_line() {
	echo $(date) "$1" >&2
}

# Record the run date for tracking.
date=$(date +'%Y-%m-%d')

log_line "Running from: ${work_dir}"
log_line "current pangenome ref: ${ref}"
log_line "current pangenome query: ${query}"
log_line "output dir for pangenome_2D_seqs: ${pangenome_2D_seqs}"
log_line "job array number: ${i}"
log_line "output file for chr2D seq: ${pangenome_2D_seqs}/${references[$reference]}_chr2D.fa"
log_line "date run started: $date"

# =============================================================================
# Extract chromosome 2D from the selected pangenome assembly
# =============================================================================

awk 'BEGIN {RS=">"} /chr2D/ {print ">"$0}' ${query} > ${pangenome_2D_seqs}/${references[$reference]}_chr2D.fa

# =============================================================================
# Run the chromosome 2D alignment
# =============================================================================

# Main directory for chromosome 2D MUMmer alignment outputs.
mummer_out=${work_dir}/mummer_out_pangenomes_2D

mkdir -p ${work_dir}/mummer_out_pangenomes_2D

# Create a separate output directory for the query assembly.
mkdir -p ${work_dir}/mummer_out_pangenomes_2D/${ref}_mummer_out

ref_dir=${work_dir}/mummer_out_pangenomes_2D/${ref}_mummer_out

# Chinese Spring chromosome 2D is used as the reference, and chromosome 2D
# from the selected pangenome assembly is used as the query. Including Chinese
# Spring in the query array would allow a self-alignment to be used as a
# control.


# Run NUCmer from the assembly-specific output directory. The -p option defines
# the prefix applied to the resulting MUMmer output files.
cd ${ref_dir}

time nucmer -p CS_2D_VS_${ref}_2D ${pangenome_2D_seqs}/chinesespring_chr2D.fa ${pangenome_2D_seqs}/${references[$reference]}_chr2D.fa

now=$(date)
log_line "nucmer completed at: $now"
