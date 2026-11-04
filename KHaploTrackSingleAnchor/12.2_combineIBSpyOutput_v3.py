## Combine output from IBSpy into one file 

#===========================
### import libraries
#============================

import os
import pandas as pd
import glob
import argparse
from datetime import datetime

pd.options.display.max_colwidth = 100



#===========================================================================================================================
## =========== Define functions to be used for analysis =================
#===========================================================================================================================

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--work_dir')
    parser.add_argument('--metadata_file')
    parser.add_argument('--reference')
    parser.add_argument('--window_size', type=int, help="window size used to run IBSpy")
    parser.add_argument('--outfile')
    parser.add_argument('--kmersums_file')
    args = parser.parse_args()
    return args



##========= Funtion 1====================
#===============================================
# define function to read the meatadata file 

pd.options.display.max_colwidth = 100

def read_lines_id(samples_id, reference):
    samples_id  = pd.read_csv(samples_id, delimiter=' ')
    samples_id = samples_id[samples_id['reference'] == reference]
    return samples_id


##========= Funtion 2 ==========================
#===============================================


def combine_db(samples_id, reference):
    import pandas as pd

    # Read the sample IDs
    samples_id = pd.read_csv(samples_id, delimiter=' ')
    samples_by_reference = samples_id[samples_id['reference'] == reference]

    # These will store all sample columns and the last row_names seen
    sample_columns = {}
    row_names = None

    # Loop through filtered samples
    for _, row in samples_by_reference.iterrows():
        variations_db = pd.read_csv(row['path'] + row['file'], delimiter='\t')
        sample_name = row['query']
        
        # Store the relevant column
        sample_columns[sample_name] = variations_db['variations']  # variations/observed_kmers Modify this line accordingly
        

        # Only need to store row names once
        if row_names is None:
            row_names = variations_db[['seqname', 'start', 'end']]

    # If no valid files were read
    if not sample_columns:
        raise ValueError("No valid samples were found for the given reference. Check if the reference matches file names.")

    # Combine everything efficiently
    samples_combined_db = pd.concat([row_names, pd.DataFrame(sample_columns)], axis=1)

    return samples_combined_db


def count_by_windows(combined_samples, window_size):
    in_db = combined_samples
    # in_db = pd.read_csv(combined_samples, delimiter='\t')
    window_size = window_size
    chrLen = in_db['end'].max() # longest chromosome
    
    w_pos = 0
    db_byChr = pd.DataFrame()
    
    while w_pos <= chrLen:
        # create a midpoint column that will be use for the aggregation instead of the end
        #especially import if using smaller window sizes
        in_db['midpoint'] = (in_db['start'] + in_db['end']) // 2
        by_windows_df = in_db[(in_db['midpoint'] > w_pos) & (in_db['midpoint'] <= w_pos + window_size)]
        
        #by_windows_df = in_db[(in_db['end'] > w_pos) & (in_db['end'] <= w_pos + window_size)]
        #by_windows_df = by_windows_df.drop(['start','end'], axis=1)
        by_windows_df = by_windows_df.drop(['start','end'], axis=1).copy()
        by_windows_df['start'] = w_pos + 1
        by_windows_df['end'] = w_pos + window_size
        w_pos += window_size
        db_byChr = pd.concat([db_byChr, by_windows_df], ignore_index=True) # replaced deprecated append function with pd.concat 
        
    by_windows_db = db_byChr.groupby(['seqname','start','end']).sum().reset_index()
    return by_windows_db



def main():
    args = parse_arguments()

    print("Combining files.....")

    metadata_file = args.metadata_file
    reference = args.reference
    window_size = args.window_size
    outfile = args.outfile
    kmersums_file = args.kmersums_file

    variations_combined = combine_db(metadata_file, reference)
    by_window_db = count_by_windows(variations_combined, window_size)
    
    # Save output to CSV file
    print(outfile)
    by_window_db.to_csv(outfile)

    #======== Generate kmer-sums file
    print("Combining files complete")


    #kmer_pre_abs = pd.read_csv(outfile,  usecols=lambda column: column != 'Unnamed: 0')
    #kmer_pre_abs['kmer_sum'] = kmer_pre_abs.iloc[:, 3:].sum(axis=1)
    #kmer_pre_abs = kmer_pre_abs.sort_values(by='kmer_sum', ascending=True)
    #kmer_pre_abs.to_csv(kmersums_file)

    #print("Generating file with the kmersums done")



if __name__ == "__main__":
    main()





