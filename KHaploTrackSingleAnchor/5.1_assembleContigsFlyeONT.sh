#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batch_logs/minimap2_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batch_logs/minimap2_\%A_\%a.err # STDERR
#SBATCH --job-name=minimap2
#SBATCH --array=0-1

# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}

now=$(date)

i=$SLURM_ARRAY_TASK_ID

#=================================================================================================
#=== Load required packages ================================
#=================================================================================================

#=== load spades
source package 3a579940-1cba-4d60-86ef-43b8705935fb # SPAdes - 3.13.1

source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6  ## seqkit 0.12.0 

source package ed35be8a-e06d-4686-b3f8-db1ad9315c7a   ## "Flye - 2.7.1"

source package 4b4051e0-bd0f-4c18-bd5e-962700bf4b28   ## Samtools 1.23

source package 222eac79-310f-4d4b-8e1c-0cece4150333   ## Minimap2 2.24-41122

#=================================================================================================
#=== Define variables and paths
#=================================================================================================

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"
cd ${workDir}


#=================================================================================================
#=== Define variables and paths
#=================================================================================================
genome="chinesespring"
region="chr4B_30860000-30910000"
REF="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation/chinesespring.fa"

# ====== array of accession from which we extracted raw reads

declare -a samples=(
    "jagger" \
    "stanley"
    )


#=======!!!!!!!!!!!!! Change here !!!!!!!!!!!===============

sample=$(($i%2))
sample=${samples[$sample]}

#!!!!!!!!!!!!!!! ======= change here ===========================================================================
##++++++++++ DIRECTORY WITH INPUT READS GROUPED BY KMER FILTER +++++++++++++++++++++++++++++++++++++++++++++++++++++
#===============================================================================================================
## ========= filtered reads directory
#retrieved_reads_dir=${workDir}/kmcOut/extrdReads/kmer_filtered_fastqs # illumina reads
retrieved_reads_dir=${workDir}/kmcOut/extrdReads/ulRreads/kmer_filtered_fastqs ## UL reads


#!!!!!!!!!!!!!!! ======= change here ===========================================================================
##++++++++++ OUTPUT DIRECTORY FOR SPADES GENERATED CONTIGS +++++++++++++++++++++++++++++++++++++++++++++++++++++
#===============================================================================================================
outDir=${workDir}/spadesOut/${sample}/${genome}

mkdir -p ${outDir}

#===============================================================================================================
#++++++++++++++++++++ CHECK VARIABLES BEING USED FOR ANAYLYSIS ARE ALL CORRECT +++++++++++++++++++++++++++++++++
#===============================================================================================================
log_line "Running from: $PWD"
log_line "slurm job ID: ${i}"
log_line "processing sample: ${sample}"
log_line "assembling sample: ${sample} with reads from reference: ${genome}-${region}"
log_line "sample output directory: ${outDir}"



#===========================================================================================================================
#+++++++++++ Process reads at each kmer filter ++++++++++++++
#===========================================================================================================================

# list of the kmer filters used to categorize the _reads

kmer_filter_list=("150_kmers_plus" "150_kmers" "140_kmers" "130_kmers" "120_kmers"  "110_kmers" "100_kmers" "90_kmers" "80kmers" "70kmers" "60kmers" "50kmers")  ## "40kmers" "30kmers" "20kmers"

#kmer_filter_list=("10kmers")
#==== process each category of _reads
for kmerfilt in ${kmer_filter_list[@]};
    # the retrieved reads for each kmer filter will be a combination of all the reads at that kmer filter from all the genomes
    do 
    #cat ${retrieved_reads_dir}/*_${kmerfilt}_RF.fastq > ${retrieved_reads_dir}/${sample}_${kmerfilt}_RF.fastq
    # rather than combine all the reads from all the genomes, we first generate individual contigs from reads of each genome anchor
    retrieved_reads="${retrieved_reads_dir}/${sample}_${kmerfilt}.fastq"

         ## ====== add a line not to process the files if it does not exist
        # Check if the retrieved_reads file exists and is not empty
    if [[ ! -s "$retrieved_reads" ]]; then
        echo "File $retrieved_reads is missing or empty. Skipping..."
        continue  # Skip to the next iteration if the file is missing or empty
    fi
 
    #=========================================================================================================================================
    #++++++++++++++ ASSEMBLE READS INTO CONTIGS WITH FLYE FROM RELEVANT KMER FILTER ++++++++++++++
     #=========================================================================================================================================


    #===================================================================================================
    #make output directory for the reads of the particular kmer filter
    #===================================================================================================
    out_dir_filt=${outDir}/${kmerfilt}
    mkdir -p ${out_dir_filt}
    mkdir -p ${out_dir_filt}/minimapOut

    #=================================
    # temp folder for spades
    #=================================
    mkdir -p ${out_dir_filt}/temp

    # setting the kmer coverage cut of to 5, the dafault is ten but only works when we have many reads like from chinese chinesespring

     log_line "Flye starting for ${sample}_${kmerfilt}"


      flye \
        --nano-raw ${retrieved_reads} \
        --genome-size 0.1m \
        --threads 24 \
        --out-dir ${out_dir_filt}


    sed "s/contig/${sample}-${kmerfilt}/g" ${out_dir_filt}/assembly.fasta > ${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_Flye_contigs.fasta

  
    log_line "Flye completed "


    ## mapping the reads to the reference to inspect the accuracy of the contigs generated
    log_line "Mapping reads to REF with minimap ... "
    
    minimap2 -t 24 \
    -I 16G \
    -ax map-ont \
    --cs=long \
    ${REF} ${retrieved_reads} \
    | samtools sort -@ 24 -o ${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}.bam

    bamFile="${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}.bam"

    samtools index --csi ${bamFile}

    # check output files
    if [[ -s ${bamFile} ]]; then
        log_line "DONE: mapping reads to ref  completed !!"
    else
        log_line "FAILED: mapping reads to ref failed, output file is missing or empty !!"
    fi
    #===========================================================================================================================
    ### Extract top five/ten contigs for each contig file assembled to be submitted for blast ======================================
    #===========================================================================================================================

    #log_line "counting contig lengths and extracting top 10 for blast"
    #### get the length of assembled contigs and save the top 10 for blast

    # ================ Input FASTA file
    #input_fasta=${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_Flye_contigs.fasta

    # ================ Output FASTA file
    #output_fasta=${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_Flye_BlastContigs.fasta

    # ================ Temporary file to store contig lengths
    #lengths_file=${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_Flye_contigs_lengths.txt

    # ================ Parse the FASTA file and calculate lengths of each contig
    #awk '/^>/ { if (seqlen) { print seqlen, seqname }; seqname=$0; seqlen=0; next } { seqlen += length($0) } END { print seqlen, seqname }' "$input_fasta" > "$lengths_file"

    # ================ Sort contigs by length and get the top 10
    #top_contigs=$(sort -k1,1nr "$lengths_file" | head -n 10 | awk '{print $2}')

    # ================ Extract the top 10 contigs from the original FASTA file
    #awk -v top_contigs="$top_contigs" '
    #BEGIN {
    #    n = split(top_contigs, contig_array, "\n")
    #    for (i = 1; i <= n; i++) {
    #        contig_hash[contig_array[i]] = 1
    #    }
    #}
    #{
    #    if ($0 ~ /^>/) {
    #        print_this = ($0 in contig_hash)
    #    }
    #    if (print_this) {
    #        print
    #    }
    #}
    #' "$input_fasta" > "$output_fasta"


    #log_line "The 10 largest contigs have been saved to: $output_fasta"

done

