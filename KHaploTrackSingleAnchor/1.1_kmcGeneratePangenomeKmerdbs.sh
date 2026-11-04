#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=19
#SBATCH --mem 50G
#SBATCH -t 5-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/generatePangenomeKmerDb_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/generatePangenomeKmerDb_\%A_\%a.err # STDERR
#SBATCH --job-name=kmcDbGen
#SBATCH --array=0-9

# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}



# ==== Objective

# generate kmerdbs of pangenomes from their WGS 

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"

cd ${workDir}


## =========== generate kmc dbs and dump db to text file ========
#==============================================================================

#=== load kmc 3.1.1

source package 88c30fea-a263-44d2-a0b6-07681251a41d

# based on the reference samples, run them in a batch array job

#=== declare the array iterator
i=$SLURM_ARRAY_TASK_ID

#location of pangenome assemblies raw reads

refs_WGS_seq="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation"


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
#"spelta.genome" \
)


reference=$(($i%10))
ref=${references[$reference]}

# name for the filtred reads file
out_name=$(echo ${ref}-pg)

#====directory to store kmc dbs generated

##instead of storing this in my personal scratch, store them with the other pangenome kmc files
out_dir="/jic/scratch/projects/watseq/kmer_agis/kmc_wheat/pangenome/kmers"

# create if it does not exist
out_dir_kmc_dbs=${out_dir}/${out_name}
mkdir -p ${out_dir_kmc_dbs}


# kmc temp dir
temp_dir=${out_dir_kmc_dbs}/${out_name}_temp
mkdir ${temp_dir}

#====for pangenomes
#create a file with the absolute path of the reads that will be used by kmc
ls -d ${refs_WGS_seq}/${ref}.fa > ${out_dir_kmc_dbs}/${out_name}_reads.txt

### monitor

log_line "kmc using reads from: ${refs_WGS_seq}"
log_line "sample: ${out_name} "
log_line "out_dir of kmc: ${out_dir_kmc_dbs}"
log_line "prefix for WGS kmcdb generated: ${out_dir_kmc_dbs}/${out_name}"
log_line "prefix for WGS MDR kmcdb generated: ${out_dir_kmc_dbs}/${out_name}-MDR"
log_line "slurm array job number: ${i}"
log_line "temp_dir kmc: $temp_dir"

## before running kmc, we check whether the kmer matrix is already available, i.e, for the WGS should already exist, but the MDR one will not exist yet
## if it exists and is not empty, then dont generate it

#=======================================================================================
## ======== Generating WGS kmer matrix for the entire the genome sequence  ===========
#=======================================================================================
WGS_kmcdb=${out_dir_kmc_dbs}/${out_name}.kmc_suf


if [[ -e "$WGS_kmcdb" && -s "$WGS_kmcdb" ]]; then
    echo ${WGS_kmcdb} "File exists and is not empty"
else
	echo ${WGS_kmcdb} "File does not exist or is empty, proceeding with kmer matrix generation"
	## ==== running kmc ========
	time kmc \
	-k31 \
	-ci1 \
	-t19 \
	-m50 \
	-fm @${out_dir_kmc_dbs}/${out_name}_reads.txt $(echo ${WGS_kmcdb} | sed 's/.kmc_suf//g') ${temp_dir} \

fi

log_line "WGS kmc db generation completed"


#=======================================================================================
## ======== Generating MDR from the genome sequence of IWGSC 2G genome of CS ===========
#=======================================================================================


WGS_MDR_kmcdb=${out_dir_kmc_dbs}/${out_name}-MDR_ci4.kmc_suf


if [[ -e "$WGS_MDR_kmcdb" && -s "$WGS_MDR_kmcdb" ]]; then
    echo ${WGS_MDR_kmcdb} "File exists and is not empty"
else
	echo ${WGS_MDR_kmcdb} "File does not exist or is empty, proceeding with kmer matrix generation"
	## ==== running kmc ========
	time kmc \
	-k31 \
	-ci2 \
	-t19 \
	-m50 \
	-fm @${out_dir_kmc_dbs}/${out_name}_reads.txt $(echo ${WGS_MDR_kmcdb} | sed 's/.kmc_suf//g') ${temp_dir} 

fi

log_line "WGS MDR with ci2 kmc db generation completed"








