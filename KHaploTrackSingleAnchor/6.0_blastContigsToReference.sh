#!/bin/bash
#SBATCH --partition=jic-short,nbi-short,jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem 100G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/blast_auto_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/blast_auto_\%A_\%a.err # STDERR
#SBATCH --job-name=blast_auto
#SBATCH --array=1



#=================================================================================================
#=== Define arrays, needed for proper naming and fetching of files
#=================================================================================================


workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"

cd ${workDir}

#===== genome being used used with the bait sequence
genome="chinesespring"

panref_genomes='/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation'

## ====== region with locus to extract (three)
#=================================================================================================

declare -a regions=("chr4B:30860000-30910000")

# ======== Define references, these are the pangenomes cultivars, used as the source of the anchor sequences
#================================================================================================================
ref_names=(
    #"arinaLrFor" \
    "chinesespring" \
    #"jagger" \
    #"julius" \
    #"lancer" \
    #"landmark" \
    #"mace" 
    #"stanley" \
    #"symattis" \
    #"norin61" 
    )

#======= Define the queries, these are the samples used to generate the conitgs
#================================================================================================= 
declare -a queries=(
    "jagger" \
    "stanley"
    )



#======== List of the kmer filters used to categorize the reads
#=================================================================================================
kmer_filter_list=("150_kmers_plus" "150_kmers" "140_kmers" "130_kmers" "120_kmers" "110_kmers" "100_kmers" "90_kmers" "80kmers" "70kmers" "60kmers" "50kmers" "40kmers" "30kmers" "20kmers")

# Loop through reference,queries and regions to generate the blast script
for ref in "${ref_names[@]}"; do
    #alt=${ref_to_alt_map["$ref"]}
    #echo "Reference: $ref, Alternative: $alt"

    for query in "${queries[@]}"; do

        sample=${query}

        for region in "${regions[@]}"; do

            region=$(echo $region | sed 's/:/_/g')
            
            out_dir=${workDir}/spades_contigs_assembly/${sample}

            blast_out=${workDir}/spades_contigs_assembly/${sample}/blast 
            mkdir -p ${blast_out}
            
            # for the sake of saving space and not doing  this again, copy the references to be used by blast into one central repo
            ref_dir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/2024-07_pangenome_assemblies"
            #mkdir -p ${ref_dir}

            #cp ${panref_genomes}/${alt}.fa ${ref_dir}

            # Process each category of reads
            for kmerfilt in ${kmer_filter_list[@]}; do
                # Create the script for the current reference/kmer filter combination
            
                script_path="${workDir}/scripts/blast_scripts/${sample}_${kmerfilt}_blast.sh"
                mkdir -p ${workDir}/scripts/blast_scripts
            
                {
                    echo "#!/bin/bash"
                    echo "#SBATCH --partition=nbi-medium,jic-medium,nbi-long,jic-long"
                    echo "#SBATCH --nodes=1"
                    echo "#SBATCH --cpus-per-task=10"
                    echo "#SBATCH --mem 100G"
                    echo "#SBATCH -o ${blast_out}/${sample}_${kmerfilt}_blast.%N.%j.out # STDOUT"
                    echo "#SBATCH -e ${blast_out}/${sample}_${kmerfilt}_blast.%N.%j.err # STDERR"
                    echo "#SBATCH --job-name=${sample}_${kmerfilt}"
                    echo "#SBATCH --array=1"
                    echo ""
                    echo "# spades directory for the reads of the particular kmer filter"
                    echo "out_dir_filt=${out_dir}/${kmerfilt}"
                    echo ""
                    echo "# the reference genome"
                    echo "reference=${ref_dir}/${ref}.fa"
                    echo "# FASTA file with contigs to use for blast"
                    echo "query=\${out_dir_filt}/${sample}_${region}_${kmerfilt}_cum_blast_contigs.fasta"
                    echo ""
                    echo "# output directory for blast alignments is the spades output directory for the reference/sample"
                    echo "# we keep all the results for the sample in one file for easier viewing on ensemble"
                    echo "blast_out=${blast_out}"
                    echo ""
                    echo "#import blast"
                    echo "source blast+-2.12.0"
                    echo ""
                    echo "#======= index reference for blast ========="
                    echo "#starting blast indexing of reference genome"
                    echo "#makeblastdb -in \${reference} -dbtype nucl"
                    echo ""
                    echo "#====== starting blast run ============"
                    echo "blastn -num_threads 10 -query \${query} -db \${reference} -task blastn-short -max_hsps 1 -max_target_seqs 1 -evalue 0.001 -outfmt 6 -out \${blast_out}/${sample}_${kmerfilt}_blast.out"
                } > $script_path

                echo "Created script: $script_path"
                sbatch $script_path
            done
        done
    done

done
