#!/bin/bash
#SBATCH --partition=jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 100G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/genMetadataIBSpy_\%A_\%a.out  # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/genMetadataIBSpy_\%A_\%a.err # STDERR
#SBATCH --array=1
#SBATCH -J genMetadataIBSpy

#======== Objective
# generate a metadata file to use for combining IBSpy outputs
## this version of the script works with the files on the cluster rather than on the pc

#===================================================================================
# ========== Define variables and file paths 
#===================================================================================

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"

cd ${workDir}

#spadesOutDir="spades_contigs_assembly"
spadesOutDir="spades_contigs_assembly"

#===== genome being used used with the bait sequence
genome="chinesespring"
region="chr4B_30860000-30910000"

reference="${genome}${region}"
## === region with locus to extract
#region="chr4B:30860000-30910000"
#region=$(echo $region | sed 's/:/-/g')


spadesOutDir="${workDir}/${spadesOutDir}"


# folder to store script and files based on ibspy data generated
out_dir=${spadesOutDir}/ibspy_files_analysis
mkdir -p ${out_dir}


#====================== folder with IBSpy output files =============================
#===================================================================================

file_dir="${spadesOutDir}/ibspy_out"

#========= list of IBSpy output files

date=$(date +'%Y-%m-%d')

file_list="${out_dir}/${date}_${reference}_IBSpy_kmers_files_list.txt"  

#========= file with metadata infor on IBSpy files from above file
metadata_file="${out_dir}/${date}_${reference}_IBSpy_kmers_metadata_file.txt"

#=== variables for the python script to generate combined file

reference="${genome}${region}"
window_size=31
outfile="${out_dir}/${reference}_presence_absence.csv"
kmersums_file="${out_dir}/${reference}_presence_absence+kmersums.csv"


#====================== generate file with the list of IBSpy output files =============================
#======================================================================================================


ls -d ${file_dir}/*.tsv.gz > ${file_list}

echo "file list save to " ${file_list}

#====================== generate file with infor needed to combine the output files =============================
#======================================================================================================

# file path we get from ls -d and the rest from the file name
# for now we are only working files from each reference, hence we choose files from all the folders in the file directory

rm ${metadata_file}

## ===== add the header line to the metadata file before adding the infor in the loop ===
echo "reference query file path"  >> ${metadata_file}

for path in $( cat ${file_list});
do
# count the number of '/' characters in the file path to use with cut to get the correct field
#
count=$(echo ${path} | tr -cd '/' | wc -c | sed 's/ //g')
count1=$((count + 1))
echo ${count} ":"${count1}

ref=$(echo ${path} | cut -f ${count1} -d "/" | cut -f 2 -d "_" | cut -f 1 -d ".")
query=$(echo ${path} | cut -f ${count1} -d "/" | cut -f 1 -d "_")
file=$(echo ${path} | cut -f ${count1} -d "/")
echo ${count}
file_path=$(echo ${path} | cut -f 1-"${count}" -d "/")
#eval file_path=$(echo "${path}" | cut -f 1-"${count}" -d "/")
echo ${ref} ${query} ${file} ${file_path}"/" >> ${metadata_file}
done

echo "file list save to "${file_list}
echo "metadata information of files saved to "${metadata_file}

#================= generate comined file with above generated information ======

echo ${metadata_file}
echo ${reference}
echo ${window_size}
echo ${outfile}
echo ${kmersums_file}

singularity exec /hpc-home/waweru/2025-01_pythonImg/2025Janpython3.img \
python3 ${workDir}/scripts/12.2_combineIBSpyOutput_v3.py \
--metadata_file ${metadata_file} \
--reference ${reference} \
--window_size ${window_size} \
--outfile ${outfile} \
--kmersums_file ${kmersums_file}