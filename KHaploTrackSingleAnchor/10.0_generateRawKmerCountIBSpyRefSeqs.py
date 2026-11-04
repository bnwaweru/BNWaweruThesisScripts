## ============ 12 December 2024 =============
# Script generates the raw kmer counts of kmers generated
# filters to keep kmers with an observed raw kmer count > 1
# kmers are saved as fasta sequences used reference sequence to run IBSpy for their presence absence
# variables required are passed dynamically to the functions in the bash script

#===========================
### import libraries
#============================

import sys
from Bio import SeqIO
import os
import gzip
import pandas as pd
import openpyxl
import argparse
import csv
from datetime import datetime


#===========================================================================================================================
## =========== Define parameters and functions to be used for analysis =================
#===========================================================================================================================
def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('--work_dir')
    parser.add_argument('--spades_dir', type=str, help="Path to directory with spades result files which includes kmer sequences and their files")
    parser.add_argument('--cultivars', type=str, help="Comma-separated list of cultivars")
    parser.add_argument('--region', type=str, help="target loci from which kmers were derived")
    parser.add_argument('--MDR', type=str, help="MDR filter used to remove repeats from anchor sequence")
    args = parser.parse_args()
    return args


#===========================================================================================================================
### Main function that generates results, no predifined functions here =====================================================================================
#===========================================================================================================================

def main():
	args = parse_arguments()

	print("Script is starting...")
	# Print to check arguments are expanded correctly
	print(f"Cultivars: {args.cultivars}, Region: {args.region}")

	spades_dir=args.spades_dir
	cultivars=args.cultivars.split(',')
	region=args.region
	MDRFilter=args.MDR
	


	for ref in cultivars:
		
		print(f"{ref}")
		
		## dicionary to keep the kmer name and seq
		uniq_kmers_no_dups_dic = {}
		kmers_raw_read_freq_dic = {}
		
		### file paths
		## once we have the non-WGS kmers to be used for the presence absence are names ina fasta file, they will be stored here
		#kmer_fasta_dir = f"{spades_dir}/{ref}/kmer_fasta_files"
		kmer_fasta_dir = f"{spades_dir}/kmer_fasta_files"
		sample_uniq_kmers_seqs=f"{kmer_fasta_dir}/{ref}_non_WGS_kmer_seqs_no_dups.fasta"
		kmer_raw_read_counts=f"{kmer_fasta_dir}/{ref}_k31_ci1_intersect_ocmax.txt"
		#kmer_raw_read_counts=f"{kmer_fasta_dir}/{ref}_k31_ci1_RF_intersect_ocmax.txt"

	    # generate the dictionary with the kmer seq and headers
		for record in SeqIO.parse(f"{sample_uniq_kmers_seqs}", 'fasta'):
			seqname=record.id
			kmer_seq=str(record.seq)
			uniq_kmers_no_dups_dic[seqname]=kmer_seq

	        # Read the kmer raw read counts into a dictionary
			with open(f"{kmer_raw_read_counts}", 'r') as counts_file:
				for line in counts_file:
					# Strip any leading/trailing whitespace or newline characters
					line = line.strip()
					
					# Split the line into key and value using '\t' (tab) as the delimiter
					if '\t' in line:  # Ensure there is a tab character in the line
						sequence, raw_count = line.split('\t', 1)  # Split at the first tab only
						# Add the key-value pair to the dictionary
						kmers_raw_read_freq_dic[sequence] = int(raw_count)
			
			## in the first dictionary the sequence is in value
	        ## in the second dictionary the sequence is the key

	        # how to combine the three into a dataframe

	        # first dataframe of kmer ID with sequence, sorted based on th kmer_seq

			kmer_seq_ID_df = pd.DataFrame(uniq_kmers_no_dups_dic.items(), columns= ['seqname', 'kmer_seq'])
			kmer_seq_ID_df_srtd = kmer_seq_ID_df.sort_values(by='kmer_seq')
			
			# second dataframe of kmer_seq with raw count, sorted based on the kmer_seq
			kmer_seq_raw_count_df = pd.DataFrame(kmers_raw_read_freq_dic.items(), columns= ['kmer_seq', 'raw_kmer_count'])
			kmer_seq_raw_count_df_srtd = kmer_seq_raw_count_df.sort_values(by='kmer_seq')
			
			# merge the two dfs based on the kmer_seq column
			kmer_ID_seq_raw_count_df = pd.merge(kmer_seq_ID_df_srtd,kmer_seq_raw_count_df_srtd, on='kmer_seq' )
			
			#save the merged dataframe to csv file
			kmer_ID_seq_raw_count_df.to_csv(f"{kmer_fasta_dir}/{ref}_{region}-{MDRFilter}_kmer_ID_seq_raw_count_df.csv")
			
			# Filter rows where raw_kmer_count > 1
			kmer_ID_seq_raw_count_df_fltrd =  kmer_ID_seq_raw_count_df[kmer_ID_seq_raw_count_df['raw_kmer_count'] > 1]
			
			# keep only sequences with raw kmer count < 255, anything with more will be repeat sequences
			kmer_ID_seq_raw_count_df_fltrdF =  kmer_ID_seq_raw_count_df_fltrd[kmer_ID_seq_raw_count_df_fltrd['raw_kmer_count'] < 255]
			
			
			#=== file where kmers for IBSpy will be saved
			IBSpy_fasta_ref=f"{kmer_fasta_dir}/{ref}_{region}-{MDRFilter}_IBSpy_kmersRefSeqs.fasta"
			
			with open(IBSpy_fasta_ref, "w") as ibspy_ref:
				
				# Write each sequence in FASTA format to the output file
				for index, row in kmer_ID_seq_raw_count_df_fltrdF.iterrows():
					#add the genome from which the kmer seq is from and the MDR filter used to generate it
					header = f">{row['seqname']}-{MDRFilter}"
					sequence = row['kmer_seq']
					ibspy_ref.write(f"{header}\n{sequence}\n")
					
			#print(kmer_seq_ID_df_srtd.head())
	    	#print(kmer_seq_raw_count_df_srtd.head())
	    	#print(kmer_ID_seq_raw_count_df.head)
			
			print(len(kmer_seq_ID_df_srtd))
			print(len(kmer_seq_raw_count_df_srtd))
			print(len(kmer_ID_seq_raw_count_df))
			print(len(kmer_ID_seq_raw_count_df_fltrd))


	print("Script complete, all raw kmer count files and IBSpy Ref file saved")


if __name__ == "__main__":
    main()

