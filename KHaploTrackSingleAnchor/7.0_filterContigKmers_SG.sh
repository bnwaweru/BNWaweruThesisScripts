#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem 10G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmers-filter_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/kmers-filter_\%A_\%a.err # STDERR
#SBATCH --job-name=kmc_kmer_filt
#SBATCH --array=0-1


##===== Objective ==========
# Filter for kmers generated from assembled contigs to check how many MDR kmers they have and remove any kmer that is present in CS
# filtering for kmers in CS in this case means that we will remove any kmers that represent the wild type CS/Stanley allele for the RHT-B1 allele,
# for future we think about this and to handle this


x=$SLURM_ARRAY_TASK_ID


# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}

now=$(date)

#!!!!!!!!!!!!!!! ======= change here ===========================================================================

#=================================================================================================
#=== Define variables and paths
#=================================================================================================

# working directory inside the folder for the current target interval spades directory
# we are working with the spades generated contigs hence that is the working dir here
workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"

cd ${workDir}

#!!!!!!!!!!!!!!! ======= change here ===========================================================================
#===============================================================================================================
# spades directory with the assembled contigs fasta

spadesOutDir="spades_contigs_assembly"

#===== genome being used used with the bait sequence
genome="chinesespring"

## === region with locus to extract
region="chr4B:30860000-30910000"
region=$(echo $region | sed 's/:/_/g')


#==========  details of the kmc db to be used for filtering of kmers ======================
# whole genome MDR kmers
WGS_MDR_kmers="/jic/scratch/projects/watseq/kmer_agis/00_kmers/pangenomes/chinesespring/chinesespring_MDR_k31_ci2"

# whole genome all kmers
WGS_k31_kmers="/jic/scratch/projects/watseq/kmer_agis/00_kmers/pangenomes/chinesespring/chinesespring_k31_ci1"
#===============================================================================================================


#++++++++++ load kmc 3.1.1  ++++++++++
#=======================================================
source package 88c30fea-a263-44d2-a0b6-07681251a41d


#===== sample being processed, the name will also be used to generate the folder for outputs


declare -a references=(
    "jagger" \
    "stanley"
    )


reference=$(($x%2))
ref=${references[$reference]}
#sample=$(echo ${ref} | cut -f 1 -d '/')  # this is for the pangenomes
sample=${ref}


 # monitoring variables
log_line "Running from: $PWD"
log_line "slurm job ID: $x"
log_line "assembling sample: ${sample}"
log_line "assembling contigs for region: ${region}"
log_line "WGS MDR kmer library: ${WGS_MDR_kmers}"
log_line "WGS all k31 kmer library: ${WGS_k31_kmers}"

# list of the kmer filters used 

kmer_filter_list=("150_kmers_plus" "150_kmers" "140_kmers" "130_kmers" "120_kmers"  "110_kmers" "100_kmers" "90_kmers" "80kmers" "70kmers" "60kmers" "50kmers" "40kmers" "30kmers" "20kmers")
#==== process each category of spades contigs
#===============================================================================================================
for kmerfilt in ${kmer_filter_list[@]};
  do contig_file="${workDir}/${spadesOutDir}/${sample}/${kmerfilt}/${sample}_${region}_${kmerfilt}_cum_blast_contigs.fasta"
  
  # folder to store the files once generated
  contig_kmers=${workDir}/${spadesOutDir}/${sample}/contig_kmers
  mkdir -p ${contig_kmers}
  cd ${contig_kmers}

  ## check the contig file exists and and isn not empty before proceeding
         ## ====== add a line not to process the files if it does not exist
        # Check if the retrieved_reads file exists and is not empty
  if [[ ! -s "$contig_file" ]]; then
      echo "File $contig_file is missing or empty. Skipping..."
      continue  # Skip to the next iteration if the file is missing or empty
  fi
   

  # using awk to print out every sequence in the contig file with its header

  awk '
      BEGIN {RS=">"; ORS=""}
      NR > 1 {
          split($0, lines, "\n");
          header = lines[1];
          sequence = "";
          for (i=2; i<=length(lines); i++) {
              sequence = sequence lines[i];
          }
          #print "Header: " header "\n";
          #print "Sequence: " sequence "\n";
          # Remove any characters from the header that are not safe for filenames
          clean_header = header;
          gsub(/[^a-zA-Z0-9_-]/, "_", clean_header);
          filename = clean_header ".fasta";
  
          # Write the header and sequence to a file
          print ">" header "\n" sequence > filename;
      }
    ' ${contig_file}


done

###============== next to generate kmc files for each of the contig_fasta files
#===============================================================================================================
contig_kmers=${workDir}/${spadesOutDir}/${sample}/contig_kmers

log_line "sample output directory: ${contig_kmers}"

for filename in ${contig_kmers}/*.fasta;
    do
    #fist check that file exsits and is not empty
    if [[ -e "$filename" && -s "$filename" ]]; then
        #echo "File '$filename' exists and is not empty."

        # generate the filename to be used to name kmer file generated
        base_name="${filename%.*}"
        echo ${base_name}
        #count the number of / characters to use to extract the correct name to name the output files
        count=$(echo ${base_name}| grep -o "/" | wc -l)
        count=$((count + 1)) # add one to get the correct field
        kmc_out_name=$(echo ${base_name} | cut -f ${count} -d "/")
        echo ${kmc_out_name}

        # folder for temp files being generated during kmc
        temp=${contig_kmers}/temp_${kmc_out_name}
        mkdir -p ${temp}

        ##============ run kmc to generate the kmer library and dump to text
        # ci1 to keep all the kmers
        ## ==== running kmc ========
        time kmc \
        -k31 \
        -ci1 \
        -t10 \
        -m100 \
        -fm ${filename} ${contig_kmers}/${kmc_out_name} ${temp}

        # dump kmers to text file
        kmc_tools transform \
	      ${contig_kmers}/${kmc_out_name} -ci1 \
	      dump -s  ${contig_kmers}/${kmc_out_name}.txt


        #============ chek how many of the kmers exist in the MDR repeat kmer lib ================
        #===============================================================================================================
        # using kmc_tools intersect which outputs kmers present in both datasets
        # when we count the kmers from this set we know which kmers are present in the MDR set of repeats

        #kmc_tools -t10 simple \
	    #	${contig_kmers}/${kmc_out_name} -ci1 ${WGS_MDR_kmers} -ci1 \
	    #	intersect ${contig_kmers}/${kmc_out_name}_MDR -ci1

        #dump kmers in MDR to text

        #kmc_tools transform \
        #${contig_kmers}/${kmc_out_name}_MDR -ci1 \
        #dump -s ${contig_kmers}/${kmc_out_name}_MDR.txt

        #============ remove any kmer that exsits in the kmer library from WGS
        #===============================================================================================================
        # using kmc_tools substract
        #only keep kmers from the first kmer matrix that are absent in the second kmer matrix
         # when we count the kmers from this set we know which kmers are unique to the accession and not part of the CS library
         # because we used use the canonical kmers, the expectation is that kmers from either the forward/reverse directions will be filtered out.

        kmc_tools -t10 simple \
	    	${contig_kmers}/${kmc_out_name} -ci1 ${WGS_k31_kmers} -ci1 \
	    	kmers_subtract ${contig_kmers}/${kmc_out_name}_MDR_no_WGS_def_kmers -ci1

         #dump kmers with CS_WGS kmers filered out to text

        kmc_tools transform \
        ${contig_kmers}/${kmc_out_name}_MDR_no_WGS_def_kmers -ci1 \
        dump -s ${contig_kmers}/${kmc_out_name}_MDR_no_WGS_def_kmers.txt


    else
        log_line "File '$filename' does not exist or is empty."
    fi
done
