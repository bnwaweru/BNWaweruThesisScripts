#!/bin/bash
#SBATCH --partition=nbi-medium,jic-medium,nbi-long,jic-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem 10G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/refKmersBlast_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000/scripts/batchLogs/refKmersBlast_\%A_\%a.err # STDERR
#SBATCH --job-name=refKmersBlast
#SBATCH --array=1

#=== initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}
#===================================================================================
## ==========Defining the working directory  ============================
#===================================================================================
#!!!!!!!!!!!!!!!!!!!!!!! Change Here

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2023_11_chr2D_CS_635_651/kmers_from_fasta/chinesespring_chr4B_30860000-30910000"
cd ${workDir}

#===== genome being used used with the bait sequence
genome="chinesespring"

## === region with locus to extract (three)
#=================================================================================================

region="chr4B_30860000-30910000"

bait_seq_dir="baitSequences"


#===spades directory
spadesOutDir=${workDir}/spades_contigs_assembly

#============== the reference sequence
#this is the sequence that was used as bait/anchor
reference=${workDir}/${bait_seq_dir}/chinesespring_chr4B_30860000-30910000.fa

#====== the query sequence
#+++++ this the file with the unique kmers generated and used to check for presence/absence with IBSpy
query=${spadesOutDir}/ibspy_out/IBSpy_kmersRefSeqs_combined_file_no_dups.fasta


#============ output file for th alignmnets

blast_out=${spadesOutDir}/ibspy_files_analysis/${genome}${region}_ContigkmersSeqsToAnchorSeq_blast.out

#import blast
source blast+-2.12.0


#============= check paths and variables are being expanded correctly ===================
#===================================================================================

log_line "Running from: $PWD"
log_line "reference: ${reference}"
log_line "query: ${query}"
log_line "output file: ${blast_out}"
log_line "blast run started"

#======= index reference for blast =========
#starting blast indexing of reference genome
makeblastdb -in ${reference} -dbtype nucl

#====== starting blast run ============
blastn \
-num_threads 10 \
-query ${query} \
-db ${reference} \
-task blastn-short \
-max_hsps 1 \
-max_target_seqs 1 \
-evalue 0.001 \
-outfmt 6 \
-out ${blast_out}


log_line "blast run completed"