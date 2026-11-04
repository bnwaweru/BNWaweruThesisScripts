# this script converts the the output of blast (output format 6) to bed file for visualization of the alignments 
# it finds all the blast outputs for each reference (per kmer filter) and saves them into one file before converting to bed file

## ======================= 
### import libraries ======
## ======================= 
import os
import argparse

#===========================================================================================================================
## ======================= Function to convert blast output to bed file and its arguments  ======================
#===========================================================================================================================

def parse_arguments():
    parser =argparse.ArgumentParser()
    parser.add_argument('--spades_dir', type=str, help="")
    parser.add_argument('--cultivars', type=str, help="Comma-separated list of pangenome cultivars")
    parser.add_argument('--kmerfilterlist', type=str, help="comma separated list of kmer filters used")
    args = parser.parse_args()
    return args


##=============================================================
##========= Funtion 1====================
#=========================================

def getBedLines(fh):
    for line in fh:
        splitline = line.split("\t")
        query_id = splitline[0]
        hit_id = splitline[1]
        perc_ID = splitline[2]
        qstart = splitline[6]
        qend = splitline[7]
        hit_start = splitline[8]
        hit_end = splitline[9]
        strand = "+"
        if int(splitline[8]) > int(splitline[9]):
            strand = "-"
        rgb = "0,255,0"

        if float(perc_ID)<99:
            rgb = "255,125,0"
       
        yield(hit_id,hit_start,hit_end,query_id,hit_start,hit_end,perc_ID,strand,rgb)



#===========================================================================================================================
### Main function that generates results with functions =====================================================================================
#===========================================================================================================================

def parse_arguments():
    parser =argparse.ArgumentParser()
    parser.add_argument('--spades_dir', type=str, help="")
    parser.add_argument('--cultivars', type=str, help="Comma-separated list of pangenome cultivars")
    parser.add_argument('--kmer_filter_list',type=str, help="comma separated list if kmer filters used")
    args = parser.parse_args()
    return args



def main():
    args = parse_arguments()

    print("Script is starting...")

    cultivars=args.cultivars.split(',')
    kmer_filter_list=args.kmer_filter_list.split(',')
    spades_dir=args.spades_dir

    for ref in cultivars:

        print(f"finding blast output files for {ref}")
        blast_out_all=[]
        blast_out_all_file=f"{spades_dir}/{ref}/{ref}_all_blast.txt"

        for kmer in kmer_filter_list:

            blast_file=f"{spades_dir}/{ref}/blast/{ref}_{kmer}_blast.out"

             ## check that the file exists and is not empty before processing
            if not os.path.isfile(blast_file) or os.path.getsize(blast_file) == 0:
                print(f"{blast_file} is missing or empty. Skipping...")
                continue  # Skip to the next iteration if the file is missing or empty

            if os.path.getsize(blast_file) > 0:
                print(f"{blast_file}")
                with open(blast_file, 'r', encoding='utf-8') as file:
                # Read the content of the file and add it to the list
                    blast_out_all.append(file.read())

            # Write the combined content to the output file
        with open(blast_out_all_file, 'w', encoding='utf-8') as output:
            for content in blast_out_all:
                output.write(content)
                #output.write("\n")  # Optional: Add a newline between file contents

        print(f"combining blast outputs for {ref} complete")

    for ref in cultivars:

        print(f"converting blast outputs for {ref} to bed file")

        blast_out_all_file=f"{spades_dir}/{ref}/{ref}_all_blast.txt"

        bed_file=f"{spades_dir}/{ref}/{ref}_all_blast.bed"

        with open(bed_file,"w") as out:
            with open(blast_out_all_file, "r") as fh:
                for hit_id,hit_start,hit_end,query_id,hit_start,hit_end,perc_ID,strand,rgb in getBedLines(fh):
                    out.write(f"{hit_id}\t{hit_start}\t{hit_end}\t{query_id}\t{perc_ID}\t{strand}\t{hit_start}\t{hit_end}\t{rgb}\n")

        print(f"converting blast outputs for {ref} to bed complete")

    print("Script complete!!")


if __name__ == "__main__":
    main()