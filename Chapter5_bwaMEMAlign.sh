#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem 150G
#SBATCH -t 10-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batchLogs/bowtieAlign_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batchLogs/bowtieAlign_\%A_\%a.err # STDERR
#SBATCH --job-name=bowtie
#SBATCH --array=0


taskID=${SLURM_ARRAY_TASK_ID:-0}
SLURM_CPUS_PER_TASK=32



# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}

#===== Objective =======

# Align reads to WATDE0042 whole genome then call SNPs
# use SNPs to design primers for genotyping the population


## ==== load tools =============
#===========================================

source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  ## bedtools 2.31.0
source package 9117dc0d-5bb0-48ce-83c3-bc0db08e0cd3  # gffread 0.12.8
source package 1041444f-cd25-4107-a5c7-5e86cb1728fe  #  seqkit 2.8.2
source package c92263ec-95e5-43eb-a527-8f1496d56f1a   ## samtools 1.18
source package fa33234e-dceb-4a58-9a78-7bcf9809edd7 ## bwa 0.7.17-r1188
source package 222eac79-310f-4d4b-8e1c-0cece4150333   ## Minimap2 2.24-41122

## =======  CONFIG ==========

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis"
cd ${workDir}

#============= location of pangenome assemblies raw reads
refs_raw_reads="/jic/scratch/projects/watseq/kmer_agis/raw_reads"

#============= location of raw reads of Watkins, both moderns and landraces
moderns_raw_reads="/jic/archive/projects/watseq"

# reference genome

REF="/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies/2025-07-02_Watkins_assemblies/IRN_WATDE0042_v1.fasta"

outDir=${workDir}/SNPCalling
bamDir=${outDir}/outBam
samtoolsOut=${outDir}/outSamtools

mkdir -p ${outDir} ${bamDir} ${samtoolsOut}


declare -a samples=(\
#"WATDE0042"
#"paragon" \
#"cadenza" \
#"norin61" \
#"julius" \
#"arinaLrFor" \
#"chinesespring" \
#"jagger" \
#"lancer" \
#"landmark" \
#"mace" \
#"stanley" \
#"symattis" \
#"WATDE0187" \
#"WATDE0076" \
#"WATDE0075" \
#"WATDE0776" \
#"WATDE0898" \
#"WATDE0552" \
#"WATDE0549" \
#"WATDE0095" \
#"WATDE0057" \
#"WATDE0031" \
#"WATDE0695" \
#"WATDE0100" \
"WATDE0228" \
#"WATDE0639" \
#"WATDE0077" \
#"WATDE0175" \
#"WATDE0062" \
#"WATDE0059"
)


sample=${samples[$taskID]}
#sample=$(($i%6))


# Define names of pangenomes array 
panref_names=("arina" \
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
    "mattis" \
	"cadenza" \
	"claire" \
	"paragon" \
	"weebil" \
    "kariega" \
    "spelt"
	)



# Check with declare -p
if [[ "$(declare -p panref_names)" =~ "${sample}" ]]; then
    log_line "${sample} is in the array."

	#create a file with the absolute path of the reads that will be used by kmc

	#====for pangenomes
	log_line "using reads from: ${refs_raw_reads}"
	ls -d ${refs_raw_reads}/${sample}/${sample}2/*.fastq.gz > "${outDir}/${sample}_reads.txt"	
else
    log_line "${sample} is not in the array."

	#for modern Watseqs/ or samples not pangenomes
	log_line "using reads from: ${moderns_raw_reads}"
	find "${moderns_raw_reads}" -type f -name "*${sample}*.gz" -print > "${outDir}/${sample}_reads.txt"
	
fi


# find the raw read files the files to save in a list to pass to bwa

#R1_ARGS=$(sort ${outDir}/${sample}_reads.txt | xargs)
#R2_ARGS=$(sort ${outDir}/${sample}_reads.txt | xargs)

# re-doing how to read the files into bwa

# Read + split
readarray -t R1_FILES < <(grep -E '.R?1([^0-9]|$)|_1\.' ${outDir}/${sample}_reads.txt | sort)
readarray -t R2_FILES < <(grep -E '.R?2([^0-9]|$)|_2\.' ${outDir}/${sample}_reads.txt | sort)

echo "R1 count: ${#R1_FILES[@]}"
echo "R2 count: ${#R2_FILES[@]}"


# Optional: print first few to verify
echo "First R1/R2 pairs:"
for i in $(seq 0 $(( ${#R1_FILES[@]} - 1 ))); do
  echo "  ${R1_FILES[$i]}  <->  ${R2_FILES[$i]}"
done

## final output bam and the temp dir for intermediate sam amd bam files

OUT="${bamDir}/${sample}.markdup.bam"
TMPDIR="${outDir}/tmp/${sample}"  # all intermidiate files with be out here, onlt the last file with duplicated marked will be placed in the bam folder
mkdir -p "$TMPDIR"


## conctenate the reads into a tmep folder to use for mapping then dlete after 

# Confirm all input files exist
for f in "${R1_FILES[@]}" "${R2_FILES[@]}"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing file: $f"; exit 1; }
done

# ---- concatenate into temp gz fastqs ----
R1_CAT="${TMPDIR}/${sample}.R1.concat.fastq.gz"
R2_CAT="${TMPDIR}/${sample}.R2.concat.fastq.gz"

log_line "Concatenating R1 -> $R1_CAT"
#cat "${R1_FILES[@]}" > "$R1_CAT"

log_line "Concatenating R2 -> $R2_CAT"
#cat "${R2_FILES[@]}" > "$R2_CAT"

# Optional quick sanity: show sizes
ls -lh "$R1_CAT" "$R2_CAT"


## monitor variables

log_line "Running from: $PWD"
log_line "job ID : $taskID"
log_line "sample : ${sample} "
log_line "sample temp folder : ${TMPDIR}"
log_line "final out bam : ${OUT}"
log_line "bwa-mem run starting "


## index reference for bwa and samtools

# BWA-MEM2 index
#bwa index ${REF}

# FASTA index for samtools/bcftools
#samtools faidx  ${REF}

log_line " Indexing complete, alignment starting"

# Read group: important if you later do multi-lane merges or GATK
RG="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA"

# Align -> sort (coordinate) -> markdup -> index
# -M marks shorter split hits as secondary (common best practice)
# Use -t for threads; samtools sort -@ for threads

#bwa mem -t ${SLURM_CPUS_PER_TASK} -M -R ${RG} ${REF} $R1_CAT $R2_CAT \
#  | samtools sort -n -@ ${SLURM_CPUS_PER_TASK} -T "${TMPDIR}/qsort" -o "${TMPDIR}/${sample}.qname.bam"


bamFile="${out_dir_filt}/minimapOut/${sample}_${genome}_${kmerfilt}.bam"

samtools index --csi ${bamFile}

# Check bwa  output
if [[ -s "${TMPDIR}/${sample}.qname.bam" ]]; then
    log_line "bwa alignment of reads for ${sample} completed !!"
else
    log_line "ERROR: bwa alignment for ${sample} failed, no bam output " >&2
    exit 1
fi


## samtools fixmate required the query to be sorted by name not coordinate

samtools sort -n -@ "${SLURM_CPUS_PER_TASK}" -T "${TMPDIR}/qsort" -o "${TMPDIR}/${sample}.qname.bam" ${TMPDIR}/${sample}.sorted.bam

# Mark duplicates (samtools markdup needs name-collated intermediate)
samtools fixmate -@ "${SLURM_CPUS_PER_TASK}" -m "${TMPDIR}/${sample}.qname.bam" "${TMPDIR}/${sample}.fixmate.bam"

# now coordinate sort
samtools sort -@ "${SLURM_CPUS_PER_TASK}" -T "${TMPDIR}/csort" -o "${TMPDIR}/${sample}.coord.bam" "${TMPDIR}/${sample}.fixmate.bam"

samtools markdup -@ "${SLURM_CPUS_PER_TASK}" -r "${TMPDIR}/${sample}.coord.bam" "$OUT"


samtools index --csi -@ "${SLURM_CPUS_PER_TASK}"  "$OUT"

# Basic stats
samtools flagstat -@ "${SLURM_CPUS_PER_TASK}" "$OUT" > "${bamDir}/${sample}.flagstat.txt"
samtools idxstats "$OUT" > "${bamDir}/${sample}.idxstats.txt"


# Check samtools
if [[ -s "$OUT" ]]; then
    log_line "samtools fixmate, sort and mark dups for ${sample} completed !!"
#    # ---- cleanup ONLY after success ----
    log_line  "Cleaning up temp files..."
    rm -f "${TMPDIR}/${sample}.fixmate.bam" "${TMPDIR}/${sample}.coord.bam" "${TMPDIR}/${sample}.qname.bam"
    rm -f "$R1_CAT" "$R2_CAT"
else
    log_line "ERROR: samtools fixmate, sort and mark dups  for ${sample} failed, final bam output missing" >&2
    exit 1
fi


## align long reads for sample WATDE0228
log_line "Aligning ${sample} long read with minimap2" 


minimap2 -t 32 \
-I 16G \
-ax map-ont \
--cs=long \
${REF} /ei/projects/a/adbcc1a6-6f74-4827-8d11-010eeea9210f/data/raw/watkins/watkins_fastq/C2_WATDE0228.fastq.gz \
| samtools sort -@ 32 -o "${bamDir}/${sample}_UL.coord.bam"


# filter and subset final bam file to read aligning to chrom5B:1-100Mb

  #4       unmapped
    #256     secondary alignment
    #2048    supplementary alignment
    #4 + 256 + 2048 = 2308

log_line "Minimap2 complete, filtering bam file."

#OUT="${bamDir}/${sample}.markdup.bam"
OUT="${bamDir}/${sample}_UL.coord.bam"

#filter and sort 
samtools view -b -q 30 -F 2308 ${OUT}  | samtools sort -o ${bamDir}/${sample}_UL.sortedMQ30.bam

#index

samtools index --csi -@ "${SLURM_CPUS_PER_TASK}" ${bamDir}/${sample}_UL.sortedMQ30.bam

#subset
samtools view -b ${bamDir}/${sample}_UL.sortedMQ30.bam "chr5B:1-100000000" > ${bamDir}/${sample}_UL.chr5B_1-100Mb.markdup.sortedMQ30.bam

# index
samtools index --csi -@ "${SLURM_CPUS_PER_TASK}" ${bamDir}/${sample}_UL.chr5B_1-100Mb.markdup.sortedMQ30.bam

fltrdFile=${bamDir}/${sample}_UL.chr5B_1-100Mb.markdup.sortedMQ30.bam.csi


# Check samtools filtering output
if [[ -s "$fltrdFile" ]]; then
    log_line "samtools filtering for ${sample} completed !!"
else
    log_line "ERROR: samtools filtering  for ${sample} failed, filtred bam output missing" >&2
    exit 1
fi