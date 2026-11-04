#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem 10G
#SBATCH -t 5-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmerFreqRawReads_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmerFreqRawReads_\%A_\%a.err # STDERR
#SBATCH --job-name=kmerFreqRawReads
#SBATCH --array=0-1

# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}

#=== declare the array iterator
i=$SLURM_ARRAY_TASK_ID

# ==== Objective

# count the frequency kmers of interest in the raw reads kmer dbs
# using kmc intersect, only kmers that actually exist in the raw data will be in the output folder

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"
cd ${workDir}

#====directory with kmc dbs generated

watkins_kmcdbs="/jic/scratch/projects/wildwheat/kmers_watkins"
pangenomes_kmcdbs="/jic/scratch/projects/watseq/kmer_agis/kmc_wheat/pangenome/kmers/"


#===== genome being used used with the bait sequence
genome="chinesespring"

## === region with locus to extract
region="chr4B:30860000-30910000"
region=$(echo $region | sed 's/:/_/g')


## =========== perform kmc intersect ========
#==============================================================================

#=== load required tools

source package 88c30fea-a263-44d2-a0b6-07681251a41d # kmc
source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6 # seqkit

#===== sample being processed, determines the path to be used for the kmcdbs

declare -a references=(
    "jagger/jagger2" \
    "stanley/stanley3"
    )

reference=$(($i%2))
ref=${references[$reference]}

# name for the directories and sample file
out_name_dir=$(echo ${ref} | cut -f 2 -d '/') ## this is only used of sample is in the pangenome array
out_name=$(echo ${ref} | cut -f 1 -d '/')


#the kmc kmer library for the watkins and the pangenomes are in differnt folders, check which it is to determine which path to use
# Define names of pangenomes array 
panref_names=("arinaLrFor" \
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

        # Check with declare -p
if [[ "$(declare -p panref_names)" =~ "${out_name}" ]]; then
    log_line "${out_name} is in the array."
	#====for pangenomes
	pangenomes_kmcdbs="/jic/scratch/projects/watseq/kmer_agis/kmc_wheat/pangenome/kmers/"
	#relevant kmc db
	kmcdb_reads="/jic/scratch/projects/watseq/kmer_agis/kmc_wheat/pangenome/kmers/${out_name_dir}/${out_name}_c1_k31_def"	

else
	log_line "${out_name} is not in the array."
	#====for watkins
	watkins_kmcdbs="/jic/scratch/projects/wildwheat/kmers_watkins"
	#relevant kmc db
	kmcdb_reads="${watkins_kmcdbs}/${out_name}/${out_name}"

fi



## folder with the kmer fasta files of non-WGS ref kmers 

spadesOutDir="spades_contigs_assembly"

kmer_seqs=${workDir}/${spadesOutDir}

# concatenate all the fasta files then generate kmer db from it to use for the intersect
#========= remove the concatenated file before re-doing the merging again

rm ${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_non_WGS_kmer_seqs.fasta

#============= merge the fasta files

cat ${kmer_seqs}/${out_name}/kmer_fasta_files/*.fasta > ${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_non_WGS_kmer_seqs.fasta

#=======remove duplicates from the file

seqkit rmdup -s ${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_non_WGS_kmer_seqs.fasta > ${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_non_WGS_kmer_seqs_no_dups.fasta

kmer_seq_file=${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_non_WGS_kmer_seqs_no_dups.fasta

#======== temp directory for kmc

temp_dir=${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_temp
mkdir -p ${temp_dir}

log_line "Running from: $PWD"
log_line "slurm job id: $i"
log_line "no dups reads: $kmer_seq_file"
log_line "sample: ${out_name}"
log_line "kmcdb of sample reads: ${kmcdb_reads}"
log_line "out_dir of kmc: ${kmer_seqs}/${out_name}/kmer_fasta_files/"
log_line "temp dir for kmc: $temp_dir"


#=================================================================================
#=========== generate the kmer db for the non duplicated kmer sequences using all reads
#=================================================================================

## running kmc 
time kmc \
-k31 \
-ci1 \
-t10 \
-m10 \
-fm ${kmer_seq_file} ${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_k31_ci1 ${temp_dir}

uniq_kmers_db=${kmer_seqs}/${out_name}/kmer_fasta_files/${out_name}_k31_ci1

## dump kmers to text file
kmc_tools transform \
    ${uniq_kmers_db} \
    dump -s ${uniq_kmers_db}.txt

#============ run kmc intersect with all raw reads =========================
#=================================================================================

# run kmc intersect 

kmc_tools -t10 simple \
	${uniq_kmers_db} -ci1 ${kmcdb_reads} -ci1 \
	intersect ${uniq_kmers_db}_intersect_ocmax -ocmax

#dump kmers to text

kmc_tools transform \
    ${uniq_kmers_db}_intersect_ocmax -ci1 \
    dump -s ${uniq_kmers_db}_intersect_ocmax.txt