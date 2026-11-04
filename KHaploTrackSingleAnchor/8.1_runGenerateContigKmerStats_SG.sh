#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 10G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/GenerateContigStatsPy_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/GenerateContigStatsPy_\%A_\%a.err # STDERR
#SBATCH --job-name=ContigStatsPy
#SBATCH --array=1

##======= Objective ==============
# run python script on the cluster using python overlay image with required libraries loaded
# the functions are defined in the python script 8.0_generateContigKmerStats_SG.py
# the variables are defined in this script and passed dynamically to the script

x=${SLURM_ARRAY_TASK_ID}


# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}
#===========================================================================================================================
## =========== Define relevant directories and variables =================
#===========================================================================================================================

### Define the working directory

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"

cd ${workDir}

### Define variables

#===== genome being used used with the bait sequence
genome="chinesespring"

## === region with locus to extract
region="chr4B:30860000-30910000"


#======= spades directory
spadesOutDir="${workDir}/spades_contigs_assembly"

#====list of kmer filters used to process reads
kmer_filter_list="140_kmers,130_kmers,120_kmers,110_kmers,100_kmers,90_kmers,80kmers,70kmers,60kmers,50kmers,40kmers,30kmers,20kmers" 

## to run this as an array job, one job per sample/ref, dec;are the sample/refs as an array

declare -a references=(
    "jagger" \
    "stanley"
    )

reference=$(($x%2))
ref=${references[$reference]}

# in bash this will be the list of genotyes/accessions that will be declared in an array to be used
# dont add space between the names as it will be used as is in the script
cultivars="jagger,stanley" ##chinese,arina


##check the varaiables 

log_line "Running from: $PWD"
log_line "slurm job ID: $x"
log_line "spades directory: ${spadesOutDir}"
log_line "sample: ${cultivars}"

singularity exec --overlay \
~/2024-02_kmers_from_fasta_analysis/2024-03-19_kmer-work-waweru.img \
~/2024-02_kmers_from_fasta_analysis/d6d23ee4-df9a-46b5-9130-46b98205fdba.img \
python3 ${workDir}/scripts/8.0_generateContigKmerStats_SG.py \
--spadesOutDir ${spadesOutDir} \
--kmer_filter_list "${kmer_filter_list}" \
--cultivars "${cultivars}"