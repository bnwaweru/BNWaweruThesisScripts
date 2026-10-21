#!/bin/bash

#===============================================================================
# SLURM job configuration
#
# Submit this script as a SLURM array job. Each task uses one compute node,
# 32 CPU cores, and 50 GB of memory. Standard output and standard error are
# written to separate log files containing the SLURM job and array-task IDs.
#===============================================================================

#SBATCH --partition=nbi-medium,jic-medium,nbi-long,jic-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 50G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batchLogs/blastPGs_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batchLogs/blastPGs_\%A_\%a.err # STDERR
#SBATCH --job-name=blastPanrefs
#SBATCH --array=1-2

# Optional SLURM directive for requesting 200 GB of local SSD scratch space.
# This line is inactive because it is written as a standard shell comment.
# S B A T C H --localscratch=ssd:200


#==========================================
#   Objective
#==========================================

# Identify the genomic locations of Chinese Spring anchor genes associated
# with the chromosome 5B QTL in the selected pangenome and Watkins assemblies.
# The whole-genome assemblies are used as BLAST nucleotide queries against a
# database containing the chromosome 5B QTL anchor-gene sequences.

#find location of CS anchor genes for the 5B QTL on releant pangenomes and


#============ load tools ===========================
#===================================================

# Load BLAST+ version 2.12.0 from the available software environment.
source blast+-2.12.0


#============= declare variables ======================
#======================================================

# Define the principal working directory for the chromosome 5B comparative
# analysis and change to this directory before running subsequent commands.
workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis"

cd ${workDir}

# Define the FASTA file containing the four Chinese Spring anchor-gene
# sequences associated with the chromosome 5B QTL. This file is used as the
# nucleotide BLAST database.

# ===== the database is the four anchor genes
database="${workDir}/anchorGeneSeqs/5B_qtl_anchor_genes.fa"


#===============================================================================
# Query-genome selection
#
# Use the SLURM array-task identifier to select one genome assembly from the
# references array. Each selected whole-genome assembly is used as the BLAST
# query against the anchor-gene sequence database.
#===============================================================================

# to submit an array job with the references as the queries,

#========== slurm id used to create the array jobs

# Store the SLURM array-task identifier in a local variable.
i=$SLURM_ARRAY_TASK_ID

# Define the genome assemblies available for analysis. Assemblies that are not
# included in the current run are retained as commented entries to document
# the broader set of pangenome references considered. In the present analysis,
# WATDE0042 and WATDE0100 are the active query genomes.
declare -a references=(\
#"arinaLrFor" \
#"chinesespring" \
#"jagger" \
#"julius" \
#"lancer" \
#"landmark" \
#"mace" \
#"norin61" \
#"stanley" \
#"symattis" \
#"paragon_v3_renamed" \
#"cadenza_v2_renamed" \
#"renan_v2.1_renamed" \
#"kariega_v1_renamed" \
#"fielder_v1_renamed" \
"WATDE0042" \
"WATDE0100" \
)

# Convert the SLURM array-task identifier to an array index using modulo 2,
# corresponding to the two active Watkins genome assemblies.
reference=$(($i%2))


#===============================================================================
# Alternative query-genome paths
#
# The commented blocks below document file locations and naming conventions
# previously used for the ten-plus wheat pangenome assemblies and additional
# downloaded assemblies. These blocks are retained for provenance but are not
# executed in the current analysis.
#===============================================================================

# location of the pan reference genome file for the 10+ genomes
#panref_genomes='/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation'
#query=${panref_genomes}/${references[$reference]}.fa
#out_file=$(echo ${query} | cut -f 7 -d "/" | cut -f 1 -d ".")

### for paragon, cadenza, renan, kariega
#panref_genomes='/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies'
#query=${panref_genomes}/${references[$reference]}.fasta
#out_file=$(echo ${query} | cut -f 8 -d "/" | cut -f 1 -d "_")


#===============================================================================
# Watkins genome inputs
#
# Define the directory containing the Watkins chromosome-scale assemblies.
# Construct the query FASTA path from the accession selected by the SLURM
# array task. Derive a concise accession identifier from the query filename
# for use in the BLAST output filename.
#===============================================================================

### For the watkins genomes
panref_genomes=/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies/2025-07-02_Watkins_assemblies
query=${panref_genomes}/IRN_${references[$reference]}_v1.fasta
out_file=$(echo ${query} | cut -f 9 -d "/" | cut -f 2 -d "_")


#===============================================================================
# Output directory
#
# Define and create the directory used to store tabular BLAST results. The
# `-p` option prevents an error when the directory already exists.
#===============================================================================

#===== main output directory for analysis data for the anchor sequencea analysis
out="${workDir}/anchorGeneSeqs/blastOut"
mkdir -p ${out}


# =============== check variables are okay ================================

#===============================================================================
# Logging function
#
# Write timestamped status messages to standard error. These messages are
# captured in the corresponding SLURM error log and document the inputs,
# outputs, and array-task identifier associated with each analysis.
#===============================================================================

#=== Initaite log line for debiging and tracking
function log_line() {
	echo $(date) "$1" >&2
}

# Record the execution date in YYYY-MM-DD format. This variable is available
# for output tracking, although it is not incorporated into the current output
# filename

# Report the principal variables used by the current array task.
log_line "Running from: ${workDir}"
log_line "blast database: ${database}"
log_line "current query: ${query}"
log_line "output dir: ${out}"
log_line "output file name: ${out_file}"
log_line "job array number: ${i}"
log_line "output file: ${out}/${out_file}_blast.out"


# =============== make the db from the Yr27 sequence  ===
#========================================================

# The following inactive command documents how the nucleotide BLAST database
# can be generated from the anchor-gene FASTA file. It is assumed that the
# required database files have already been created before this script is run.

#makeblastdb \
# -in ${database} -dbtype nucl


#===============================================================================
# BLASTN analysis
#
# Search each selected whole-genome assembly against the nucleotide database
# containing the chromosome 5B QTL anchor genes. Because the genome assembly
# is supplied as the query and the anchor genes form the database, the output
# reports genome-derived query regions that match the anchor-gene sequences.
#
# The analysis uses 32 threads, the standard `blastn` task, an E-value
# threshold of 0.001, and permits up to three high-scoring segment pairs per
# query–subject pair. Results are written in tabular BLAST format 6.
#===============================================================================

blastn \
 -num_threads 32 \
 -query ${query}\
 -task blastn \
 -db ${database}\
 -outfmt 6 \
 -max_hsps  3 \
 -evalue 0.001 \
 -out ${out}/${out_file}_hits.out
