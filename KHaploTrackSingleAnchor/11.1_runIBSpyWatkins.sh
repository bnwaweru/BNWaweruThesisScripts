#!/bin/bash
#SBATCH --partition=jic-long,nbi-long,jic-medium,nbi-medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 350G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/IBSpyWatkins_\%A_\%a.out  # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/IBSpyWatkins_\%A_\%a.err # STDERR
#SBATCH --array=2-1052
#SBATCH -J IBSpyWatkins

# S B A T C H --localscratch=ssd:200

i=$SLURM_ARRAY_TASK_ID

source IBSpy-0.3.1

function log_line() {
	echo $(date) "$1" >&2
}


#===================================================================================
## ==========Defining the working directory  ============================
#===================================================================================
#!!!!!!!!!!!!!!!!!!!!!!! Change Here

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"

cd ${workDir}

#spadesOutDir="spades_contigs_assembly"
spadesOutDir="spades_contigs_assembly"

#===== genome being used used with the bait sequence
genome="chinesespring"

## === region with locus to extract, here the region is composed of all the three regions we were looking at
region="chr4B_30860000-30910000"


# Objective

#!!!!!!!!!!!!!!!!!!!!!!! Change Here =======================
#===================================================================================
#=======Define the reference to be used for IBSpy ==================
#===================================================================================

ibspy_reference=${workDir}/${spadesOutDir}/ibspy_out/IBSpy_kmersRefSeqs_combined_file_no_dups.fasta

#this file is quite big because it has kmers unique from every pangenome cultivar, before running IBSpy we removed the duplicates
#from 18,186 down to 8666


window_size=50000

#===================================================================================
## ==========Defining the samples to be used for IBSpy ============================
#===================================================================================


# the kmc kmer matrix for the samples, this is for all the watkins, both he moderns and the landraces
# kmer-matrices for the pangenomes
#/Volumes/Project-Scratch/watseq/kmer_agis/kmc_wheat/pangenome/kmers/
# they are 1046 watkins accessions in total, use a map file submit the job for each
db_path="/jic/scratch/projects/wildwheat/kmers_watkins"

#=== for the watkins
mapfile WatkinsList < /jic/scratch/groups/Cristobal-Uauy/waweru/2024-02_IBSpy_NA/kmer_analysis_project/scripts/WatkinsNamesList.txt

db=${db_path}/$(echo ${WatkinsList[$i]} | sed 's/ //g')/${WatkinsList[$i]}


# check variable assignment is as expected
echo $PWD
echo ${db_path}
echo ${db}
echo ${WatkinsList[$i]}
echo $(echo ${WatkinsList[$i]} | sed 's/ //g')
echo ${reference}

#===================================================================================
## ==========Defining the IBSpy output directory  ============================
#===================================================================================
# output directory for particular sample_reference combination
out_dir=${workDir}/${spadesOutDir}/ibspy_out
mkdir -p $out_dir


#============= check paths and variables are being expanded correctly ===================
#===================================================================================

log_line "Running from: $PWD"
log_line "slurm job id: $i"
log_line "Reference: $ibspy_reference"
log_line "name: ${WatkinsList[$i]}"
log_line "db: $db"
log_line "output director: $out_dir"
log_line "window size in kb: $window_size"
log_line "IBSpy run started"

#======= run IBSpy =============================
# date of run
run_date=$(date +'%d-%m-%Y')

# check if the output file exists and is not empty, if it does not then proceed with IBSpy runs
file_path=${out_dir}/$(echo ${WatkinsList[$i]} | sed 's/ //g')_${genome}${region}_kmers.tsv.gz

#======= check output file name is correct
log_line "IBSpy output file name: $(echo ${file_path} | sed 's/.gz//g') "


if [[ -e "$file_path" && -s "$file_path" ]]; then
    echo ${file_path} "File exists and is not empty" >> ${workDir}/${run_date}_IBSPy_runs_info.txt
else
		echo ${file_path} "File does not exist or is empty, proceeding with IBSpy run" >> ${workDir}/${run_date}_IBSPy_runs_info.txt
		time IBSpy \
		--kmer_size 31 \
		--window_size ${window_size} \
		--reference ${ibspy_reference} \
		--database $db \
		--database_format kmc3 \
		--output $(echo ${file_path} | sed 's/.gz//g') \
		--compress
fi

#--database_format jellyfish

now=$(date)
log_line "IBSpy completed"
