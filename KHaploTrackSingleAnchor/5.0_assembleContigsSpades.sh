#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 50G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/SPAdes_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/SPAdes_\%A_\%a.err # STDERR
#SBATCH --job-name=SPAdes
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

source package 46a62eca-4f8f-45aa-8cc2-d4efc99dd9c6   ## seqkit 0.12.0 

source package ed35be8a-e06d-4686-b3f8-db1ad9315c7a   ## load Flye

source package 222eac79-310f-4d4b-8e1c-0cece4150333   ## Minimap2 2.24-41122

source package 4b4051e0-bd0f-4c18-bd5e-962700bf4b28   ## Samtools 1.23

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
retrieved_reads_dir=${workDir}/kmcOut/extrdReads/kmer_filtered_fastqs # illumina reads
#retrieved_reads_dir=${workDir}/kmcOut/extrdReads/ulRreads/kmer_filtered_fastqs ## UL reads


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


#==============================================================================================================================================
#++++++++++++++++++++ CHECK IF SAMPLE IS PANGENOME OR WATKINS TO DETERMINE SOURCE OF READS AND HOW THEY ARE PROCESSED +++++++++++++++++++++++++
#==============================================================================================================================================

#the reads from pangenomes and those from watkins have different headers, check which it is to determine which commands to separate into F and R reads
# Define names of pangenomes array 
panref_names=("arinaLrFor" \
    "chinesespring" \
    "chinese" \
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

#===========================================================================================================================
#+++++++++++ Extract forward and reverse reads from each fastq file with of each kmer file for each reference ++++++++++++++
#===========================================================================================================================

# list of the kmer filters used to categorize the _reads

kmer_filter_list=("150_kmers_plus" "150_kmers" "140_kmers" "130_kmers" "120_kmers"  "110_kmers" "100_kmers" "90_kmers" "80kmers" "70kmers" "60kmers") ## "50kmers" "40kmers" "30kmers" "20kmers"

#kmer_filter_list=("10kmers")
#==== process each category of reads, from each kmer_filter
for kmerfilt in ${kmer_filter_list[@]};
   
    do 
    retrieved_reads="${retrieved_reads_dir}/${sample}_${kmerfilt}_RF.fastq"

         ## ====== add a line not to process the files if it does not exist
        # Check if the retrieved_reads file exists and is not empty
    if [[ ! -s "$retrieved_reads" ]]; then
        echo "File $retrieved_reads is missing or empty. Skipping..."
        continue  # Skip to the next iteration if the file is missing or empty
    fi

        # Check with declare -p
    if [[ "$(declare -p panref_names)" =~ "${sample}" ]]; then
        log_line "${sample} is in the array."

        #create a file with the absolute path of the reads that will be used by kmc

        #====for pangenomes
        log_line "extracting forward reads for ${sample}_${kmerfilt}"
        ## ==== extract forward reads
        awk 'NR%4 == 1 {if (/\/1$/) {f=1; r=0} else {f=0; r=1}} f' ${retrieved_reads} > ${outDir}/${sample}_${kmerfilt}_${genome}_R1.fastq

        log_line "extracting reverse reads for ${sample}_${kmerfilt}"
        ## ==== extract reverse reads
        awk 'NR%4 == 1 {if (/\/2$/) {f=0; r=1} else {f=1; r=0}} r' ${retrieved_reads} > ${outDir}/${sample}_${kmerfilt}_${genome}_R2.fastq

    else
        log_line "${sample} is not in the array."

            #for modern Watseqs/ or samples not pangenomes
        log_line "extracting forward and reverse reads"

        while IFS= read -r header && IFS= read -r sequence && IFS= read -r plus && IFS= read -r quality;
        do
            # Process header, sequence, plus, and quality lines here
            echo "Header: $header"
            echo "Sequence: $sequence"
            echo "Plus: $plus"
            echo "Quality: $quality"
            R1=$(echo $header | cut -d " " -f 2 | cut -d ":" -f 1)
                if [[ $R1 == 1 ]]; then
                        printf "$header\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R1.fastq
                        printf "$sequence\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R1.fastq
                        printf "$plus\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R1.fastq
                        printf "$quality\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R1.fastq
                elif [[ $R1 == 2 ]]; then
                        printf "$header\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R2.fastq
                        printf "$sequence\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R2.fastq
                        printf "$plus\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R2.fastq
                        printf "$quality\n" >> ${outDir}/${sample}_${kmerfilt}_${genome}_R2.fastq
                fi
        done < ${retrieved_reads}
        
    fi

    # forwad and reverse reads to be used in the assembly
    R1_reads="${outDir}/${sample}_${kmerfilt}_${genome}_R1.fastq"
    R2_reads="${outDir}/${sample}_${kmerfilt}_${genome}_R2.fastq"


   
    #=========================================================================================================================================
    #++++++++++++++ ASSEMBLE READS INTO CONTIGS WITH SPADes USING EXTRACTED FORWARD AND REVERSE READS FROM RELEVANT KMER FILTER ++++++++++++++
     #=========================================================================================================================================
    log_line "SPAdes started for ${sample}_${kmerfilt}"


    #===================================================================================================
    #make spades directory for the reads of the particular kmer filter
    #===================================================================================================
    out_dir_filt=${outDir}/${kmerfilt}
    mkdir -p ${out_dir_filt}
    mkdir -p ${out_dir_filt}/minimapOut

    #=================================
    # temp folder for spades
    #=================================
    mkdir -p ${out_dir_filt}/temp

    # setting the kmer coverage 
    #the dafault is ten but only works when we have many reads like from chinese chinesespring
    # cut of to 5 fi doing a re-teration or very few reads are expected
    # adding the sc option worked to assemble contigs from super few reads

    #--only-assembler
    #  --careful \

    spades.py --sc --only-assembler \
        --cov-cutoff 3 \
        -o ${out_dir_filt} \
        -1 ${R1_reads} \
       -2 ${R2_reads} \
        -t 24 \
        -m 50 \
        --tmp-dir ${out_dir_filt}/temp \
        -k 111,121,127


    # rename assembled contigs appropriately
    # the contig will be renamed based on the region, as given in the region variable at the beginning of the script
    # in the file while renaming also replace the NODE on the contig name with the sample name for reference while using the contigs

    sed "s/NODE/${sample}-${kmerfilt}-NODE/g" ${out_dir_filt}/contigs.fasta > ${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_SPAdes_contigs.fasta

    log_line "SPAdes completed "


    ## mapping the reads to the reference to inspect the accuracy of the contigs generated
    ## option is to also run blast, if so follow commented steps below
    log_line "Mapping reads to REF with minimap ... "

    REF="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation/chinesespring.fa"

    minimap2 -t 24 \
    -I 16G \
    -ax sr \
    ${REF} ${R1_reads} ${R2_reads} \
    | samtools sort -@ 24 -o ${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}_reads.bam

        ## keep only mapping with MAPQ > 30 and primary alighments

    #4       unmapped
    #256     secondary alignment
    #2048    supplementary alignment
    #4 + 256 + 2048 = 2308

    samtools view -b -q 30 -F 2308 ${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}_reads.bam \
    | samtools sort -o ${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}_reads.sortedMQ30Sort2.bam

    samtools index --csi ${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}_reads.sortedMQ30Sort2.bam

    bamFile="${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}_reads.sortedMQ30Sort2.bam"


    # check output files
    if [[ -s ${bamFile} ]]; then
        log_line "DONE: mapping reads to ref  completed !!"
    else
        log_line "FAILED: mapping reads to ref failed, output file is missing or empty !!"
    fi
   
   
   
   #===========================================================================================================================
    ### Extract top five/ten contigs for each contig file assembled to be submitted for blast ======================================
    #===========================================================================================================================

   # log_line "counting contig lengths and extracting top 10 for blast"
    #### get the length of assembled contigs and save the top 10 for blast

    # ================ Input FASTA file
    #input_fasta=${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_SPAdes_contigs.fasta

    # ================ Output FASTA file
    #output_fasta=${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_mafft_contigs.fasta

    # ================ Temporary file to store contig lengths

    #lengths_file=${out_dir_filt}/${sample}_${genome}_${kmerfilt}_cum_SPAdes_contigs_lengths.txt
    # ================ Parse the FASTA file and calculate lengths of each contig
    #awk '/^>/ { if (seqlen) { print seqlen, seqname }; seqname=$0; seqlen=0; next } { seqlen += length($0) } END { print seqlen, seqname }' "$input_fasta" > "$lengths_file"

    # ================ Sort contigs by length and get the top 10,2 
    #top_contigs=$(sort -k1,1nr "$lengths_file" | head -n 2 | awk '{print $2}')

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