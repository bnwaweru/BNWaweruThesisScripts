#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 150G
#SBATCH -t 10-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/freebayes_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/freebayes_\%A_\%a.err # STDERR
#SBATCH --job-name=callVariants
#SBATCH --array=0

i=$SLURM_ARRAY_TASK_ID
slurmCPUs=24



# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}

#===== Objective =======

# Call Variants from reads aligned to Resistant parent genome



## ==== load tools =============
#===========================================

source package 4b4051e0-bd0f-4c18-bd5e-962700bf4b28   ## Samtools 1.23
source package 43f9ae9c-a231-41d7-9d07-142ea9637a4d   ## freebayes 1.2.0
# ======== CONFIG ============

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis"
cd ${workDir}

## variables and paths

REF="/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies/2025-07-02_Watkins_assemblies/IRN_WATDE0042_v1.fasta"

bamDir=${workDir}/SNPCalling/outBam
bcfOut=${bamDir}/freebayesOut
OUTDIR=${bcfOut}/5B_1-100Mb

mkdir -p ${bcfOut}  ${OUTDIR}


# make a list of the samples
#for file in $(ls ${bamDir}/*.markdup.bam); 
#do
#    sample=$(basename "$file")
#    sample=${sample%.markdup.bam}

    #delete the sample and bam list file its not just adding lines to an laready existing file
    #rm ${bcfOut}/samples.txt. ${bcfOut}/bams.txt

#    echo "${sample}" >> ${bcfOut}/samples.txt 
#    echo "${bamDir}/${sample}.markdup.bam" >> ${bcfOut}/bams.txt

    # check if index file exists else create it
#    idxFile=${bamDir}/${sample}.markdup.bam.csi 
#    if [[ -s ${idxFile} ]]; then
#	    log_line "index file for ${sample} exists!!"
#    else
#        log_line "index file for ${sample} does not exist, proceeding with indexing ....!!"
#        samtools index --csi -@ "${SLURM_CPUS_PER_TASK}"  ${file}
#    fi
#done 

samples=${bcfOut}/samples.txt 
bamFilesList=${bcfOut}/bams.txt
REGION="chr5B:1-100000000"

## Output file variables

RAW_VCF="${OUTDIR}/5B_1_100Mb.raw.vcf"
NORM_VCF="${OUTDIR}/5B_1_100Mb.norm.vcf.gz"
BSA_VCF="${OUTDIR}/5B_1_100Mb.bsa.vcf.gz"
KASP1_VCF="${OUTDIR}/5B_1_100Mb.kasp.stage1.vcf.gz"
KASP1_VCFTags="${OUTDIR}/5B_1_100Mb.kasp.tags.vcf.gz"
KASP_VCF="${OUTDIR}/5B_1_100Mb.kasp.final.vcf.gz"



#============= check paths and variables are being expanded correctly ===================
#===================================================================================

log_line "Running from: $PWD"
log_line "slurm job id: $i"
log_line "Reference: ${REF}"
log_line "calling variants for: ${REGION}"
log_line "output directory: ${OUTDIR}"
log_line "Freebayes run started....."

# set variable for running freebayes

### freebyaes allows for settng of ploidy 6 for wheat 
PLOIDY=6   

# Tune these to your data:
# -q 30 helps avoid homeolog/paralog multi-mapping in wheat subgenomes
# -Q 20 is a good baseline for Illumina
MAPQ=30
BASEQ=20

## running freebayes

# FreeBayes polyploid call
#freebayes \
#  -f ${REF} \
#  -L ${bamFilesList} \
#  -r ${REGION} \
#  -p ${PLOIDY} \
#  --genotype-qualities \
#  --strict-vcf \
#  --min-base-quality ${BASEQ} \
#  --min-mapping-quality ${MAPQ} \
#  --min-alternate-count 3 \
#  --min-alternate-fraction 0.05 \
#  --use-best-n-alleles 4 \
#  -g 500
#  > ${RAW_VCF}


log_line "compressing ${RAW_VCF}" 

bgzip -f "$RAW_VCF"
tabix -f -p vcf "${RAW_VCF}.gz"


# check output files
if [[ -s "${RAW_VCF}.gz" ]]; then
	log_line "DONE: Freebayes call variants successfully completed !!"
else
	log_line "FAILED: Freebayes call variants failed, output file is missing or empty !!"
    exit 1
fi

# next is to decompose the SNPs into individual ones where there are complex SNPs and normalise split multi-alleles


log_line " Normalising and splitting multi-alleles ..."

# normailze variants with bcftools, these will ensure consistent variant corrdinated wand also split mutliallelic alleles

# -Oz \                                           # write compressed VCF
bcftools norm --threads ${slurmCPUs} \
  -f ${REF} \
  -m -any \
  ${RAW_VCF}.gz \
  -Oz -o ${NORM_VCF}

# index the output file for easy referencing 

tabix -f -p vcf "${NORM_VCF}"

# check output files
if [[ -s "${NORM_VCF}" ]]; then
	log_line "DONE: Normalizing variants successfully completed !!"
else
	log_line "FAILED: Normalizing variants failed, output file is missing or empty !!"
    exit 1
fi


# -------------------------------
# 1) Stricter SNP set for KASP
# -------------------------------
  # -v snps \                 # keep only SNPs
  # -m2 \                     # minimum of 2 alleles total
  # -M2 \                     # maximum of 2 alleles total = biallelic only
  # --SnpGap 20 \             # filter SNPs if they are within 20 bp of an indel
 # -e 'QUAL<30' \            # exclude sites with QUAL < 30 (stricter than BSA set)

bcftools view --threads ${slurmCPUs} \
  -v snps \
  -m2 \
  -M2 \
  ${NORM_VCF} \
| bcftools filter  --threads ${slurmCPUs}\
  --SnpGap 20 \
| bcftools filter  --threads ${slurmCPUs}\
  -e 'QUAL<30' \
  -Oz \
  -o ${KASP1_VCF}


## index the vcf file

bcftools index --threads ${slurmCPUs}\
  -f \
  -t \
  ${KASP1_VCF}

# check output files
if [[ -s "${KASP1_VCF}" ]]; then
	log_line "DONE: Filtering for KASP genotyping SNPs completed !!"
else
	log_line "FAILED: Filtering for KASP genotyping SNPs failed, output file is missing or empty !!"
    exit 1
fi

# -------------------------------
# 2) Add cohort-level summary tags
# ------------------------------

  bcftools +fill-tags --threads ${slurmCPUs} "${KASP1_VCF}" \
  -Oz \
  -o "${KASP1_VCFTags}" \
  -- -t AN,AC,AF,MAF,NS,F_MISSING                           # tags to calculate and add

bcftools index --threads ${slurmCPUs} \
  -f \
  -t \
  "${KASP1_VCFTags}"


# check output files
if [[ -s "${KASP1_VCFTags}" ]]; then
	log_line "DONE: Adding cohort level summaries to KASP VCF completed !!"
else
	log_line "FAILED: Adding cohort level summaries to KASP VCF failed, output file is missing or empty !!"
    exit 1
fi

# -------------------------------
# 3) Final KASP-friendly filter
# ------------------------------

# -i 'MAF>=0.05 && F_MISSING<=0.20' \   # keep only sites meeting both conditions

bcftools view --threads ${slurmCPUs} \
  -i 'MAF>=0.05 && F_MISSING<=0.20' \
  -Oz \
  -o "$KASP_VCF" \
  "${KASP1_VCFTags}"


bcftools index --threads ${slurmCPUs} \
  -f \
  -t \
  "$KASP_VCF"

bcftools stats "$KASP_VCF" > "${OUTDIR}/5B_1_100Mb.kasp.stats.txt"


# check output files
if [[ -s ${KASP_VCF} ]]; then
	log_line "DONE: Final KASP VCF generation completed !!"
else
	log_line "FAILED: Final KASP VCF generation failed, output file is missing or empty !!"
fi



### remov any SNPs that have a nother SNP 50bp up or dowstream of it


awk 'BEGIN{OFS="\t"}
{
    chrom[NR]=$1
    pos[NR]=$2
    line[NR]=$0
}
END{
    for(i=1;i<=NR;i++){
        keep=1
        if(i>1 && chrom[i]==chrom[i-1] && (pos[i]-pos[i-1])<=50) keep=0
        if(i<NR && chrom[i]==chrom[i+1] && (pos[i+1]-pos[i])<=50) keep=0
        if(keep) print line[i]
    }
}' ${KASP_VCF} > snps.no_nearby_50bp.tsv


awk '{print $1"\t"$2}' snps.no_nearby_50bp.tsv > keep_sites.tsv

# then using bcf tool to get the final set of SNPs

bcftools view -T keep_sites.tsv snps.only.vcf.gz -Oz -o snps.kasp_50bp_filtered.vcf.gz
bcftools index snps.kasp_50bp_filtered.vcf.gz


## convert SNPs to bed intervals so we cn extract flamking sequences

bcftools query -f '%CHROM\t%POS\n' snps.kasp_50bp_filtered.vcf.gz > snps.pos
