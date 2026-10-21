#!/bin/bash

#===============================================================================
# SLURM job configuration
#
# Submit this script as a SLURM array job for NLR annotation of chromosome 5B
# sequences. Each array task processes one sequence listed in sampleList.txt.
#
# The job requests:
#   - one compute node;
#   - 20 CPU cores;
#   - 100 GB of memory;
#   - a maximum runtime of 30 days and 12 hours.
#
# Standard output and standard error are written to separate files containing
# the SLURM job and array-task identifiers.
#===============================================================================

#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem 100G
#SBATCH -t 30-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/NLRAnnotator_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/NLRAnnotator_\%A_\%a.err # STDERR
#SBATCH --job-name=NLRAnnot
#SBATCH --array=1-4


#===============================================================================
# Logging function
#
# Write timestamped status messages to standard error. These messages are
# captured in the SLURM error log and provide a record of the selected query,
# input sequence, output prefix, and progress through the annotation workflow.
#===============================================================================

#=== initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}


#===============================================================================
# Load software
#
# Load the NLR-Annotator software package and the Java Runtime Environment
# required to execute the Java archive files used by ChopSequence,
# NLR-Parser, and NLR-Annotator.
#===============================================================================

#===== load package =======
#NLR Annotator
source package 64ec168a-e76a-455f-92fa-8cb0e1c06644
source jre-7.21 #java


#===================================================================================
## ==========Defining the working directory  ============================
#===================================================================================
#!!!!!!!!!!!!!!!!!!!!!!! Change Here

# Define the principal working directory for the chromosome 5B comparative
# analysis and change to this directory before running the workflow.
work_dir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis"

cd ${work_dir}

# Define the main output directory for NLR-Annotator results.
annotator_dir=${work_dir}/NLRAnnotatorOut

# Alternative directory containing chromosome 5B query sequences. This path is
# retained as a comment to document an earlier or alternative input layout.
#queries_dir=${work_dir}/chr5BSeqs/queries

# Directory containing the chromosome 5B reference sequences analysed by the
# current version of the script.
queries_dir=${work_dir}/chr5BSeqs/refs


#===============================================================================
# Select the sequence processed by the current SLURM array task
#
# sampleList.txt contains one sample or assembly identifier per line. The
# SLURM array-task identifier is used as the line number, allowing each task
# to process a different chromosome 5B FASTA file.
#===============================================================================

#========== Define variables required ======
sampleList=${work_dir}/scripts/sampleList.txt
query=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${sampleList} )

# Create a sample-specific output directory.
# make query specific out dir
mkdir -p ${annotator_dir}/${query}

# Define the chromosome 5B input FASTA file and the common prefix used for all
# outputs produced for the selected sample.

# input and output
input=${queries_dir}/${query}_chr5B.fasta
outPrefix=${annotator_dir}/${query}/${query}_chr5B

# Recreate the sample-specific output directory if it is not already present.
# The -p option prevents an error when the directory already exists.

# make query specif out dir
mkdir -p ${annotator_dir}/${query}

# Historical variables retained as comments from an earlier workflow version.
#mot=${annotator_dir}/mot.txt
#store=${annotator_dir}/store.txt


#===============================================================================
# Report the variables used by the current array task
#===============================================================================

##------- Monitor variables -------------
#------------------------------------------

# Record the working directory, SLURM task identifier, selected query,
# chromosome FASTA input, and output prefix in the job log.
log_line "current working dir: $PWD"
log_line "slum job ID: ${SLURM_ARRAY_TASK_ID}"
log_line "query: ${query}"
log_line "input fasta: ${input}"
log_line "ouput file prefix: ${outPrefix}"


#===============================================================================
# NLR-Annotator software and resource paths
#
# Define the paths to:
#   - the MEME motif model used to identify NLR-associated motifs;
#   - the MAST executable used for motif searches;
#   - ChopSequence, which divides the input sequence into overlapping chunks;
#   - NLR-Parser, which identifies and organises NLR motif combinations;
#   - NLR-Annotator, which predicts and reports NLR loci.
#===============================================================================

## NLR annotator variables

meme="/jic/scratch/projects/training/steuernb/meme.xml"
mast="/nbi/software/testing/meme/4.9.1/x86_64/bin/mast"
chopper="/jic/scratch/projects/training/steuernb/ChopSequence.jar"
parser="/jic/scratch/projects/training/steuernb/NLR-Parser.jar"
annotator="/jic/scratch/projects/training/steuernb/NLR-Annotator.jar"

# Define the FASTA file that will contain the overlapping sequence fragments
# generated from the chromosome 5B input.
CHOPPED=${outPrefix}_chopped.fasta


#------------------------------------------------------
##------- Step 1: Chopping up the genome -------------
#------------------------------------------------------

# Divide the chromosome 5B FASTA sequence into 20-kb fragments with a 5-kb
# overlap between adjacent fragments. Overlapping fragments allow candidate
# NLR loci located near fragment boundaries to be represented in more than one
# input segment.

# chops the genome into 20kb chuncks with a 5kb overlap

log_line "chopping up input genome: ${input}"

java -XX:ParallelGCThreads=20 -jar ${chopper} -i ${input} -o ${CHOPPED} -l 20000 -p 5000


#---------------------------------------------------------------
##------- Step 2: Run NLR Parser to identify motifs -------------
#-----------------------------------------------------------------

# Run NLR-Parser on the chopped chromosome sequence. MAST searches each
# sequence fragment against the NLR-associated motif model in meme.xml.
# Candidate motif combinations and their positions are written to an XML
# intermediate file.

log_line "running NLR-Parser"

java -XX:ParallelGCThreads=2 -jar ${parser} -t 20 -y ${mast} -x ${meme} -i ${CHOPPED} -c ${outPrefix}_nlr.xml


#---------------------------------------------------------------
##------- Step 3: Run NLR Annotator  -------------
#-----------------------------------------------------------------

# Process the NLR-Parser XML output with NLR-Annotator to identify candidate
# NLR loci and generate tabular, BED, motif, alignment, and FASTA outputs.
#
# Output files include:
#   ${outPrefix}.txt
#       Tabular summary of predicted NLR loci.
#
#   ${outPrefix}.bed
#       Genomic coordinates of predicted NLR loci.
#
#   ${outPrefix}.motifs.bed
#       Coordinates of individual NLR-associated motifs.
#
#   ${outPrefix}.nbarcmotifalignment.fasta
#       FASTA alignment of NB-ARC-associated motifs.
#
#   ${outPrefix}.nlr.fasta
#       FASTA sequences of predicted NLR loci, including the specified
#       flanking sequence length.
#
# The final positional argument of 1000 specifies the number of flanking bases
# included when extracting predicted NLR sequences.

log_line "running NLR-Annotator"

java -XX:ParallelGCThreads=10 -jar $annotator \
	-i ${outPrefix}_nlr.xml \
	-o ${outPrefix}.txt \
	-b ${outPrefix}.bed \
	-m ${outPrefix}.motifs.bed \
	-a ${outPrefix}.nbarcmotifalignment.fasta \
	-f ${input} ${outPrefix}.nlr.fasta 1000


#===============================================================================
# Compress and preserve intermediate output
#
# Change to the main NLR-Annotator output directory, compress the XML file
# produced by NLR-Parser, and create a copy labelled as an intermediate file.
#===============================================================================

cd ${annotator_dir}
gzip ${outPrefix}_nlr.xml
cp ${outPrefix}_nlr.xml ${outPrefix}_intermediate.xml.gz


#===============================================================================
# Validate final output
#
# Confirm that the principal tabular NLR-Annotator output exists and is not
# empty. A successful run reports completion; otherwise, the script exits with
# a non-zero status so that SLURM records the task as failed.
#===============================================================================

# check that outfile is generated and is not empty
if [[ -s "${outPrefix}.txt" ]]; then
    echo "NLR Annotator for ${query} Complete!!"
else
    echo "ERROR: Output file ${outPrefix}.txt not found or empty!"
    exit 1
fi
