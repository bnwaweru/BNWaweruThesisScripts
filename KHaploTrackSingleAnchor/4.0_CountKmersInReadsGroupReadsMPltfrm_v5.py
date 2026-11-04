## ============ 28 October 2024 =============
# Script to count the number of kmers in fastq files and to group reads based on kmers to use for assembly of reads
# this is version two of the script
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
    parser.add_argument('--reads_dir', type=str, help="Path to reads directory")
    parser.add_argument('--kmc_kmer_filter', type=int, help="minimum number of kmers used as a filter to extract reads with kmc")
    parser.add_argument('--kmer_db_file')
    parser.add_argument('--kmer_db_file_no_repeats', type=str, required=True, help="Path to kmer db txt file without repeats used to extract read using kmc")
    parser.add_argument('--cultivars', type=str, help="Comma-separated list of cultivars")
    parser.add_argument('--genome', type=str, required=True, help="name of genome used to extract the reads for the cultivar")
    parser.add_argument('--read_stats_csv', type=str, help="path to csv file to save statistics of reads retrieved per ref/kmer_filter combination as they are processed")
    args = parser.parse_args()
    return args

##========= Funtion 1====================
#===========================================================================================================================

def read_kmers_list(kmer_file):
    """read kmer text file line by line removing any characters, i.e white spaces and end of line characters"""
    with open(kmer_file, "r")as file:
        kmer_list=[line.split("\t")[0] for line in file.readlines()]
    return kmer_list

##========= Funtion 2 ====================
#===========================================================================================================================

def count_reads_fastq(file_path):
    """count the number of reads in a fastq file"""
    if file_path.endswith('.gz'):
        read_count = 0
        with gzip.open(file_path, 'rt') as file:
            for record in SeqIO.parse(file, "fastq"):
                read_count += 1
        return read_count
    else:
        read_count = 0
        with open(file_path, 'rt') as file:
            for record in SeqIO.parse(file, "fastq"):
                read_count += 1
        return read_count

##========= Funtion 3 ====================
#===========================================================================================================================

def count_kmers_in_read(read_sequence, kmer_database):
    """count the number of kmers in a read"""
    #args = parse_arguments()
    #kmer_database=set(read_kmers_list(args.kmer_db_file_no_repeats))
    count = 0
    k = len(next(iter(kmer_database)))  # Assumes all k-mers in the database have the same length

    for i in range(len(read_sequence) - k + 1):
        kmer = read_sequence[i:i+k]
        if kmer in kmer_database:
            count += 1

    return count

##========= Function 4 ====================
#===========================================================================================================================

def kmer_count_file(fastq_file, kmer_database):
    """ generate a lists that will store the readIDs and the kmer counts that can then be saved in a database """
    
    args = parse_arguments()
    #kmer_database=set(read_kmers_list(args.kmer_db_file_no_repeats))
    #print(kmer_database[:5])

    #== make directory to save the kmer count files within the reads directory
    kmer_counts_dir=f"{args.reads_dir}/kmer_count_files_MPltfrm"
    os.makedirs(kmer_counts_dir, exist_ok=True)
    
    #== place holder for the ref variable so that the scrip will run

    # lists to be used
    readIDs=[]
    kmerCounts=[]
    
    with open(fastq_file, "r") as fastq_file:
        for record in SeqIO.parse(fastq_file, "fastq"):
            read_sequence = str(record.seq)
            kmer_count = count_kmers_in_read(read_sequence, kmer_database)

            #print(f"Read {record.id}: Number of k-mers from database found = {kmer_count}")
            readIDs.append(record.id)
            kmerCounts.append(kmer_count)
            
    df=pd.DataFrame(columns=['readID', 'kmer_count'])
    for read_id, kmer_count in zip(readIDs, kmerCounts):
        #df=df.append({'readID': read_id, 'kmer_count': kmer_count}, ignore_index=True) ## changed this to updated funtion from .append to pd.concat
        df = pd.concat([df, pd.DataFrame([{'readID': read_id, 'kmer_count': kmer_count}])], ignore_index=True)

        
    return df

##========= Function 5 ====================
#===========================================================================================================================

from Bio import SeqIO

def kmer_count_and_extract_reads_cum(fastq_file, kmer_database):
    
    args = parse_arguments()

    #kmer_database=set(read_kmers_list(args.kmer_db_file_no_repeats))
    #print(list((kmer_database)[:5]))

    # Initialize lists to store read IDs based on k-mer count ranges
    readIds_20kmers = []
    readIds_30kmers = []
    readIds_40kmers = []
    readIds_50kmers = []
    readIds_60kmers = []
    readIds_70kmers = []
    readIds_80kmers = []
    readIds_90kmers = []
    readIds_100kmers = []
    readIds_110kmers = []
    readIds_120kmers = []
    readIds_130kmers = []
    readIds_140kmers = []
    readIds_150kmers_plus = []
    
    # Open the FASTQ file and parse each record
    with open(fastq_file, "r") as fastq:
        for record in SeqIO.parse(fastq, "fastq"):
            # Convert sequence to string
            read_sequence = str(record.seq)
            
            # Count k-mers in the read sequence
            kmer_count = count_kmers_in_read(read_sequence, kmer_database)
            #print(f"Read {record.id}: Number of k-mers from database found = {kmer_count}")
            # Determine the range of k-mer count and append read ID accordingly
            if kmer_count <= 20:
                readIds_20kmers.append(record.id)
                #print(f"Read {record.id}: Number of k-mers from database found = {kmer_count}")
            elif kmer_count >= 21 and  kmer_count <= 30:
                readIds_30kmers.append(record.id)
                #print(f"Read {record.id}: Number of k-mers from database found = {kmer_count}")
            elif kmer_count >= 31 and  kmer_count <= 40:
                readIds_40kmers.append(record.id)
                #print(f"Read {record.id}: Number of k-mers from database found = {kmer_count}")
            elif kmer_count >= 41 and  kmer_count <= 50:
                readIds_50kmers.append(record.id)
                #print(f"Read {record.id}: Number of k-mers from database found = {kmer_count}")
            elif kmer_count >= 51 and  kmer_count <= 60:
                readIds_60kmers.append(record.id)
            elif kmer_count >= 61 and  kmer_count <= 70:
                readIds_70kmers.append(record.id)
            elif kmer_count >= 71 and  kmer_count <= 80:
                readIds_80kmers.append(record.id)
            elif kmer_count >= 81 and  kmer_count <= 90:
                readIds_90kmers.append(record.id)
            elif kmer_count >= 91 and  kmer_count <= 100:
                readIds_100kmers.append(record.id)
            elif kmer_count >= 101 and  kmer_count <= 110:
                readIds_110kmers.append(record.id)
            elif kmer_count >= 111 and  kmer_count <= 120:
                readIds_120kmers.append(record.id)
            elif kmer_count >= 121 and  kmer_count <= 130:
                readIds_130kmers.append(record.id)
            elif kmer_count >= 131 and  kmer_count <= 140:
                readIds_140kmers.append(record.id)
            elif kmer_count >= 141 :
                readIds_150kmers_plus.append(record.id)
                #print(f"Read {record.id}: Number of k-mers from database found = {kmer_count}")
                # Add additional conditions for other k-mer count ranges if necessary

            ## generate list of cumulative reads starting with the hghest kmer filter.
            # Subsequent lower kmer filter set will have its reads plus the reads for the all kmer filters above it

                # Initialize lists to store read IDs based on k-mer count ranges
        readIds_20kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers+readIds_80kmers+readIds_70kmers+readIds_60kmers+readIds_50kmers+readIds_40kmers+readIds_30kmers+readIds_20kmers
        readIds_30kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers+readIds_80kmers+readIds_70kmers+readIds_60kmers+readIds_50kmers+readIds_40kmers+readIds_30kmers
        readIds_40kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers+readIds_80kmers+readIds_70kmers+readIds_60kmers+readIds_50kmers+readIds_40kmers
        readIds_50kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers+readIds_80kmers+readIds_70kmers+readIds_60kmers+readIds_50kmers
        readIds_60kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers+readIds_80kmers+readIds_70kmers+readIds_60kmers
        readIds_70kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers+readIds_80kmers+readIds_70kmers
        readIds_80kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers+readIds_80kmers
        readIds_90kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers+readIds_90kmers
        readIds_100kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers+readIds_100kmers
        readIds_110kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers+readIds_110kmers
        readIds_120kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers+readIds_120kmers
        readIds_130kmers_cum = readIds_150kmers_plus+readIds_140kmers+readIds_130kmers
        readIds_140kmers_cum = readIds_150kmers_plus+readIds_140kmers
        readIds_150kmers_plus = readIds_150kmers_plus
            
    # Return lists of read IDs
        return readIds_150kmers_plus,readIds_140kmers_cum,readIds_130kmers_cum,readIds_120kmers_cum,readIds_110kmers_cum,readIds_100kmers_cum,readIds_90kmers_cum,readIds_80kmers_cum,readIds_70kmers_cum,readIds_60kmers_cum,readIds_50kmers_cum,readIds_40kmers_cum,readIds_30kmers_cum,readIds_20kmers_cum

##========= Function 6 ====================
#===========================================================================================================================
## function to retrieve reads by reads IDs

def retrieve_records_by_id(file_path, read_ids):
    """ from a fastq file retrieve/extract reads based on a list or read_IDs"""
    records = []
    with open(file_path, 'rt')  as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in read_ids:
                records.append(record)
    return records


##========= Function 7 ====================
#===========================================================================================================================
## function to save stats of reads as they are generated into a csv file


# Function to write stats into the CSV file
def write_stats_to_csv(file_name, data):
    file_exists = False
    try:
        # Check if file exists to determine whether to write a header
        with open(file_name, 'r') as f:
            file_exists = True
    except FileNotFoundError:
        file_exists = False

    # Open the file in append mode
    with open(file_name, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Write header if file doesn't already exist
        if not file_exists:
            writer.writerow(["Reference","Total_Read_IDs", "Total_FR_reads","kmer_filt", "Reads_without_mates", "Reads_with_mates", "AverageReadlength"])

        # Write the data row
        writer.writerow(data)



def average_read_length(records) -> float:
    if not records:
        return 0.0
    return sum(len(rec.seq) for rec in records) / len(records)


#===========================================================================================================================
### Main function that generates results with functions =====================================================================================
#===========================================================================================================================



def parse_arguments():
    parser = argparse.ArgumentParser()
    #parser.add_argument('--work_dir')
    parser.add_argument('--reads_dir')
    parser.add_argument('--kmc_kmer_filter', type=int, help="minimum number of kmers used as a filter to extract reads with kmc")
    parser.add_argument('--kmer_db_file')
    parser.add_argument('--kmer_db_file_no_repeats')
    parser.add_argument('--cultivars', type=str, help="List of cultivars")
    parser.add_argument('--genome', type=str, required=True, help="name of genome used to extract the reads for the cultivar")
    parser.add_argument('--read_stats_csv', type=str, help="path to csv file to save statistics or reads retrieved per ref/kmer_filter combination as they are processed")
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    print("Script is starting...")
    print(f"Cultivars: {args.cultivars}")
    print(f"Genome: {args.genome}")

    # make output directories
    kmer_counts_dir = f"{args.reads_dir}/kmer_count_files_MPltfrm"
    os.makedirs(kmer_counts_dir, exist_ok=True)

    kmer_filtrd_fastqs = f"{args.reads_dir}/kmer_filtered_fastqs_MPltfrm"
    os.makedirs(kmer_filtrd_fastqs, exist_ok=True)

    # reset stats file
    file_path = args.read_stats_csv
    if os.path.exists(file_path):
        print("Stats File exists. Deleting it.")
        os.remove(file_path)
    else:
        print("Stats File does not exist.")

    # load kmers
    kmer_list = read_kmers_list(args.kmer_db_file)
    print(f"Number of kmers in bait sequence: {len(kmer_list)}")

    kmer_list_no_repeats = read_kmers_list(args.kmer_db_file_no_repeats)
    print(f"Number of kmers in bait sequence after filtering out repeat kmers: {len(kmer_list_no_repeats)}")

    # IMPORTANT: use a set/list of kmers, not the filename
    kmer_database = set(kmer_list_no_repeats)

    cultivars = args.cultivars.split(',')
    genome = args.genome

    # ---------------------------------------------------------
    # 1. Count kmers in each read and save count tables
    # ---------------------------------------------------------
    for ref in cultivars:
        #ref_fastq = f"{args.reads_dir}/{ref}-{genome}_{args.kmc_kmer_filter}kmers_non_repeat_target_filtered.fastq"
        ref_fastq=f"{args.reads_dir}/{ref}_{args.kmc_kmer_filter}kmers_non_repeat_target_filtered.fastq"

        if not os.path.exists(ref_fastq):
            print(f"WARNING: missing input FASTQ for {ref}: {ref_fastq}")
            continue

        num_reads = count_reads_fastq(ref_fastq)
        print(f"number of reads retrieved for {ref} is : {num_reads}")

        kmer_count_df = kmer_count_file(ref_fastq, kmer_database)
        kmer_count_df.to_csv(f"{kmer_counts_dir}/{ref}_reads_kmer_counts.csv", index=False)

    # ---------------------------------------------------------
    # 2. Build cumulative groups and extract reads
    #    - paired branch if RF file exists
    #    - long-read branch if RF file missing
    # ---------------------------------------------------------
    for ref in cultivars:
        #ref_fastq = f"{args.reads_dir}/{ref}-{genome}_{args.kmc_kmer_filter}kmers_non_repeat_target_filtered.fastq"
        ref_fastq=f"{args.reads_dir}/{ref}_{args.kmc_kmer_filter}kmers_non_repeat_target_filtered.fastq"
        #ref_fastq_RF = f"{args.reads_dir}/{ref}-{genome}_{args.kmc_kmer_filter}kmers_non_repeat_target_filtered_sequences_RF.fastq"
        ref_fastq_RF=f"{args.reads_dir}/{ref}_{args.kmc_kmer_filter}kmers_non_repeat_target_filtered_sequences_RF.fastq"

        if not os.path.exists(ref_fastq):
            print(f"WARNING: missing input FASTQ for {ref}: {ref_fastq}")
            continue

        num_reads = count_reads_fastq(ref_fastq)

        (
            readIds_150kmers_plus,
            readIds_140kmers_cum,
            readIds_130kmers_cum,
            readIds_120kmers_cum,
            readIds_110kmers_cum,
            readIds_100kmers_cum,
            readIds_90kmers_cum,
            readIds_80kmers_cum,
            readIds_70kmers_cum,
            readIds_60kmers_cum,
            readIds_50kmers_cum,
            readIds_40kmers_cum,
            readIds_30kmers_cum,
            readIds_20kmers_cum
        ) = kmer_count_and_extract_reads_cum(ref_fastq, kmer_database)

        readIDsList = [
            readIds_150kmers_plus,
            readIds_140kmers_cum,
            readIds_130kmers_cum,
            readIds_120kmers_cum,
            readIds_110kmers_cum,
            readIds_100kmers_cum,
            readIds_90kmers_cum,
            readIds_80kmers_cum,
            readIds_70kmers_cum,
            readIds_60kmers_cum,
            readIds_50kmers_cum,
            readIds_40kmers_cum,
            readIds_30kmers_cum,
            readIds_20kmers_cum
        ]

        kmer_filt = [
            "150_kmers_plus", "140_kmers", "130_kmers", "120_kmers",
            "110_kmers", "100_kmers", "90_kmers", "80kmers",
            "70kmers", "60kmers", "50kmers", "40kmers",
            "30kmers", "20kmers"
        ]

        previous_reads = 0

        # ---------------------------
        # PAIRED-END branch
        # ---------------------------
        if os.path.exists(ref_fastq_RF):
            print(f"{ref} has paired-end reads.")

            num_readsRF = count_reads_fastq(ref_fastq_RF)
            print(f"number of reads retrieved for {ref} is : {num_reads} and in the file with mates is : {num_readsRF}")

            for kmer, read_id_list in zip(kmer_filt, readIDsList):
                reads_without_mates = len(read_id_list)

                if reads_without_mates == previous_reads:
                    print(f"No new reads to retrieve with {kmer} ({reads_without_mates})")
                    continue

                print(f"retrieving reads with {kmer} ({reads_without_mates}) and their mates")

                records = retrieve_records_by_id(ref_fastq_RF, read_id_list)
                reads_with_mates = len(records)
                avg_len = round(average_read_length(records), 2)

                print(f"Number of reads without mates retrieved with {kmer} is {reads_without_mates} and with mates is {reads_with_mates}")

                output_file = f"{args.reads_dir}/kmer_filtered_fastqs/{ref}_{kmer}_RF.fastq"
                with open(output_file, "w") as handle:
                    SeqIO.write(records, handle, "fastq")

                row_data = [
                    ref,                 # Reference
                    num_reads,           # total selected read IDs
                    num_readsRF,         # total reads in RF file
                    kmer,                # kmer filter
                    reads_without_mates, # selected read IDs
                    reads_with_mates,    # records written
                    avg_len              # average read length
                ]
                write_stats_to_csv(args.read_stats_csv, row_data)

                previous_reads = reads_without_mates

        # ---------------------------
        # LONG-READ branch
        # ---------------------------
        else:
            print(f"{ref} has long reads.")

            print(f"number of reads retrieved for {ref} is : {num_reads}")

            for kmer, read_id_list in zip(kmer_filt, readIDsList):
                reads_selected = len(read_id_list)

                if reads_selected == previous_reads:
                    print(f"No new reads to retrieve with {kmer} ({reads_selected})")
                    continue

                print(f"retrieving reads with {kmer} ({reads_selected})")

                # IMPORTANT: pass the actual list, not f"{read_id_list}"
                records = retrieve_records_by_id(ref_fastq, read_id_list)
                reads_written = len(records)
                avg_len = round(average_read_length(records), 2)

                print(f"Number of reads retrieved with {kmer} is {reads_written}")

                output_file = f"{args.reads_dir}/kmer_filtered_fastqs/{ref}_{kmer}.fastq"
                with open(output_file, "w") as handle:
                    SeqIO.write(records, handle, "fastq")

                row_data = [
                    ref,            # Reference
                    num_reads,      # total selected reads
                    "None",         # no RF file
                    kmer,           # kmer filter
                    reads_selected, # selected read IDs
                    reads_written,  # records written
                    avg_len         # average read length
                ]
                write_stats_to_csv(args.read_stats_csv, row_data)

                previous_reads = reads_selected

    print("Script complete")

if __name__ == "__main__":
    main()