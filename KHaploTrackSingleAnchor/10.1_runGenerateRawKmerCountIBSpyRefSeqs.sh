#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 10G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/RawKmerCountIBSpyRefGen_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/RawKmerCountIBSpyRefGen_\%A_\%a.err # STDERR
#SBATCH --job-name=RawKmerCount


##======= Objective ==============
# run python script on the cluster using python overlay image with required libraries loaded
# the functions are defined in the python script Step10_generateRawKmerCount_IBSpyRefSeqs_v2.py
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

#regions=(chr4B:30860000-30910000 chr4B:30860000-30910000 chr4B:30860000-30910000)


#for region in "${regions[@]}"; 
#do echo "Processing ${region}"
#region=$(echo $region | sed 's/:/_/g')
region="chr4B_30860000-30910000"

#===spades directory
spades_dir=${workDir}/spades_contigs_assembly

# in bash this will be the list of genotyes/accessions that will be declared in an array to be used
# dont add space between the names as it will be used as is in the script
cultivars="jagger,stanley" #chinese,jagger,lancer,landmark,mace,stanley,mattis,norin61,julius,arina

#====== for the region that will used to name files, we combine the genome and region into one

fileIdentifier=${genome}-${region}

echo ${fileIdentifier}

singularity exec --overlay \
~/2024-02_kmers_from_fasta_analysis/2024-03-19_kmer-work-waweru.img \
~/2024-02_kmers_from_fasta_analysis/d6d23ee4-df9a-46b5-9130-46b98205fdba.img \
python3 ${workDir}/scripts/10.0_generateRawKmerCountIBSpyRefSeqs.py \
--spades_dir ${spades_dir} \
--cultivars "${cultivars}" \
--region ${region}


#done

## find all the combined IBSpy RefSeqs fasta files and save into one file
## because we used all the kmers combined, we will have an ibspy_out folder outside of the spdes directories this time.

mkdir -p ${workDir}/spades_contigs_assembly/ibspy_out

cat $(find ${workDir}/spades_contigs_assembly/ -type f -name '*IBSpy_kmersRefSeqs.fasta') > ${workDir}/spades_contigs_assembly/ibspy_out/IBSpy_kmersRefSeqs_combined_file.fasta