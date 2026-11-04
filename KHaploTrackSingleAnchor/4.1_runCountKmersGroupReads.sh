#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 10G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/readFilterPy_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/readFilterPy_\%A_\%a.err # STDERR
#SBATCH --job-name=ReadFilterPy


##======= Objective ==============
# run python script on the cluster using python overlay image with required libraries loaded
# the functions are defined in the python script 4.0_CountKmersInReadsGroupReadsMPltfrm_v5.py
# the variables are defined in this script and passed dynamically to the script


#===========================================================================================================================
## =========== Define relevant directories and variables =================
#===========================================================================================================================

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"
cd ${workDir}


#===== genome being used used with the bait sequence
genome="chinesespring"

## === region with locus to extract (three)
#=================================================================================================

region="chr4B_30860000-30910000"


#anchor kmer databases,
#======= before filtering of repeats
kmer_db_file="${workDir}/kmcOut_${region}/${genome}_${region}.txt"

#======= after filtering out repeats
## bait kmc db filtered for repeats used is below;
kmer_db_file_no_repeats="${workDir}/kmcOut_${region}/${genome}_${region}-target-filtered-no-repeats.txt"

## ========= filtered reads directory
reads_dir=${workDir}/kmcOut_${region}/${genome}-filtered_reads # illumina reads

#file to save statistics of reads retrieved per ref/kmer_filter combination as they are processed
read_stats_file=${workDir}/kmcOut_${region}/${genome}_${region}_read_statsMPltfrm.csv


# in bash this will be the list of genotyes/accessions that will be declared in an array to be used
# dont add space between the names as it will be used as is in the script
 #,jagger,lancer,landmark,mace,stanley,mattis,norin61,julius,arina
cultivars="jagger,stanley"

#4.0_CountKmersInReadsGroupReadsMPltfrm_v5.py
#4.0_countKmersInReadsGroupReads_Oldv3.py

singularity exec /hpc-home/waweru/2025-01_pythonImg/2025Janpython3.img \
python3 ${workDir}/scripts/4.0_CountKmersInReadsGroupReadsMPltfrm_v5.py \
--reads_dir ${reads_dir} \
--kmc_kmer_filter 10 \
--kmer_db_file ${kmer_db_file} \
--kmer_db_file_no_repeats ${kmer_db_file_no_repeats}  \
--genome ${genome} \
--cultivars "${cultivars}" \
--read_stats_csv "${read_stats_file}"