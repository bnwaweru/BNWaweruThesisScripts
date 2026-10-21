#!/bin/bash

#===============================================================================
# SLURM job configuration
#
# This script is submitted as a SLURM array job. It requests one compute node,
# two CPU cores, 250 GB of memory, and a maximum runtime of 100 days and 12
# hours. Standard output and error messages are written to separate log files.
#===============================================================================

#SBATCH --partition=nbi-medium,jic-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 250G
#SBATCH -t 100-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/nucmerPG_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/nucmerPG_\%A_\%a.err # STDERR
#SBATCH --job-name=NUcmerPG
#SBATCH --array=0

# Assign the SLURM array-task identifier to a local variable. This value is
# subsequently used to select an entry from the reference-genome array.
i=$SLURM_ARRAY_TASK_ID


#===============================================================================
# Logging function
#
# Write timestamped status messages to standard error. These messages are
# captured in the SLURM error log and provide a record of job progress,
# selected inputs, and any errors encountered during execution.
#===============================================================================

# initiate log reporting
function log_line() {
	echo $(date) "$1" >&2
}


##==== Objective =====

# Extract the first 100 Mb of chromosome 5B from the selected pangenome
# assembly and compare it with the corresponding chromosome 5B interval from
# WATDE0042. Whole-region sequence alignment is performed using NUCmer to
# examine structural and sequence-level differences among the assemblies.

#-- Extract the first 100 Mb of the W-042 and W-100, align to chr5B of the pangenomes to see how they differ


## ==== load tools =============
#===========================================

# Load the software packages required for sequence extraction, indexing, and
# whole-genome alignment. Package identifiers correspond to modules available
# within the institutional computing environment.

source package eab121cb-2eb8-49c1-a9a5-a33754ea9fea  ## MUMmer - 4.0.0beta2
source package 222eac79-310f-4d4b-8e1c-0cece4150333  ## minimap2
source package c92263ec-95e5-43eb-a527-8f1496d56f1a   ## samtools 1.18


## =======  CONFIG ==========
#===========================================

# Define the principal working directory and change to this location before
# running subsequent commands.
workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis"
cd ${workDir}

# Define the directory in which NUCmer alignment outputs will be stored.
nucmerOutDir="${workDir}/nucmerResults"

# Create the NUCmer output directory if it does not already exist.
mkdir -p ${nucmerOutDir}

# Define directories containing the chromosome 5B query and reference
# sequences used for comparative analysis.

# query 100 Mb seq dir
queryDir=${workDir}/chr5BSeqs/queries
refDir=${workDir}/chr5BSeqs/refs

# Alternative directory containing the original pangenome reference
# assemblies. This path is retained for documentation but is not currently
# used by the script.
#pgRefDir="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation"

# Directory containing the downloaded and locally prepared genome assemblies.
pgRefDir="/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies"   # for the downloaded assemblies


#===============================================================================
# Reference-genome selection
#
# Define the set of pangenome assemblies available for comparison. Assemblies
# that are currently excluded from the analysis are retained as commented
# entries to document the complete set of potential references. In the present
# run, only the Renan v2.1 assembly is active.
#===============================================================================

# submit batch jobs based on the pangenome refs

declare -a references=(\
#"WATDE0100" \
#"arinaLrFor" \
#"chinesespring" \
#"julius" \
#"jagger" \
#"lancer" \
#"landmark" \
#"mace" \
#"norin61" \
#"stanley" \
#"symattis" \
#"paragon_v3" \
#"cadenza_v2" \
#"fielder_v1" \
#"kariega_v1" \
"renan_v2.1"
)

# Use the SLURM array-task identifier to calculate the index of the reference
# assembly to be processed by the current job.
reference=reference=$(($i%3))

# Retrieve the reference-assembly name from the references array.
ref=${references[$reference]}

# Set WATDE0042 as the accession against which the selected pangenome assembly
# will be aligned.
query="WATDE0042"


#===============================================================================
# Reference-genome input
#
# Define the FASTA file containing the selected pangenome assembly. The
# alternative naming convention for previously used assemblies is retained as
# a commented command for provenance.
#===============================================================================

#refGenome=${pgRefDir}/${ref}.fasta  # .fa for the 10+ genomes and .fasta for the downloaded assemblies
refGenome=${pgRefDir}/${ref}_renamed.fasta


#=== Extract Chr5B of pangenome =====

# Build a coordinate-sorted CSI index for the selected reference FASTA file.
# The index enables samtools faidx to retrieve a specified genomic interval.
samtools index --csi ${refGenome}

# Extract bases 1–100,000,000 from chromosome 5B of the selected pangenome
# assembly and write the interval to an individual FASTA file.
samtools faidx ${refGenome} chr5B:1-100000000 > ${refDir}/${ref}_chr5B_1-100Mb.fasta


#==== Nucmer inputs =====

# Define the chromosome 5B FASTA file from the selected pangenome assembly.
# The full-chromosome input used in an earlier version of the analysis is
# retained as a commented line for provenance.
#ref5BFasta=${refDir}/${ref}_chr5B.fasta
ref5BFasta=${refDir}/${ref}_chr5B_1-100Mb.fasta

# Define the WATDE0042 chromosome 5B query sequence corresponding to the first
# 100 Mb of the chromosome.
queryFasta=${queryDir}/${query}_chr5B_1-100Mb.fasta


##------- Monitor variables -------------
#------------------------------------------

# Record the principal job variables in the SLURM log. These messages document
# the working directory, array-task identifier, selected genome assembly,
# query accession, reference name, and output directory used for each run.
log_line "current working dir: $PWD"
log_line "slum job ID: ${i}"
log_line "current genome ref: ${refDir}/${ref}_chr5B.fasta"
log_line "current query: ${query}"
log_line "current reference: ${ref}"
log_line "output directory: ${nucmerOutDir}"


#===============================================================================
# Input-file validation
#
# Confirm that both FASTA inputs exist and contain data before initiating the
# NUCmer alignment. The job terminates with a non-zero exit status when either
# input is absent or empty.
#===============================================================================

# Check input files
if [[ ! -s "$ref5BFasta" ]]; then
    log_line "ERROR: Reference file $ref5BFasta missing or empty" >&2
    exit 1
fi

if [[ ! -s "$queryFasta" ]]; then
    log_line "ERROR: Query file $queryFasta missing or empty" >&2
    exit 1
fi

# Record the start of the NUCmer alignment in the job log.
log_line "NUcmer run starting ....."


#========== run NUcmer =======================
#==============================================

# Previous alignment orientation in which the selected pangenome assembly was
# supplied as the NUCmer reference and WATDE0042 was supplied as the query.
# These commands are retained as comments to document the earlier workflow.

#time nucmer --prefix=${nucmerOutDir}/${ref}_${query} ${ref5BFasta} ${queryFasta}

# Check nucmer delta output
#if [[ -s "${nucmerOutDir}/${ref}_${query}.delta" ]]; then
#    log_line "NUCmer alignments completed !!"
#else
#    log_line "ERROR: Nucmer failed, no delta output for ${ref} vs ${query}" >&2
#    exit 1
#fi


##====== switching it so that WATDE0042 is the ref
#==============================================

# Perform the chromosome 5B alignment using WATDE0042 as the first NUCmer
# input and the selected pangenome assembly as the second input. The `time`
# command records the elapsed time and computational resource use of the
# alignment. NUCmer writes the resulting alignment to a delta file with the
# prefix WATDE0042_<reference>.
time nucmer --prefix=${nucmerOutDir}/${query}_${ref} ${queryFasta} ${ref5BFasta}


#===============================================================================
# Output validation
#
# Confirm that NUCmer produced a non-empty delta file. Successful completion is
# reported in the log; otherwise, the job terminates with a non-zero status.
#===============================================================================

# Check nucmer delta output
if [[ -s "${nucmerOutDir}/${query}_${ref}.delta" ]]; then
    log_line "NUCmer alignments completed !!"
else
    log_line "ERROR: NUCmer failed, no delta output for ${query} vs ${ref} " >&2
    exit 1
fi
