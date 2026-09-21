#!/bin/bash
#SBATCH --partition=nbi-medium,jic-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 250G
#SBATCH -t 100-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/SNPsFilterBlast_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batch_logs/SNPsFilterBlast_\%A_\%a.err # STDERR
#SBATCH --job-name=SNPBlast
#SBATCH --array=1

i=$SLURM_ARRAY_TASK_ID
slurmCPUs=24


# initiate log reporting

function log_line() {
	echo $(date) "$1" >&2
}

#===== Objective =======

# filter SNPs for KASPs, removing any occuring within 50 bp of another
# extract flanks then blast to remove any off target SNPs


## ==== load tools =============
#===========================================

source package 4b4051e0-bd0f-4c18-bd5e-962700bf4b28   ## Samtools 1.23
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478   ## bedtools 2.31.0
source package 59f39368-c0c2-4511-94d6-85fdcd062755   ## blast+ 2.17.0 

# ======== CONFIG ============

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis"
cd ${workDir}

## variables and paths

REF="/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies/2025-07-02_Watkins_assemblies/IRN_WATDE0042_v1.fasta"

outDir=${workDir}/SNPCalling/outBam/freebayesOut/5B_1-100Mb
KASP_VCF="${outDir}/5B_1_100Mb.kasp.final.vcf.gz"

cd ${outDir} 

#============= check paths and variables are being expanded correctly ===================
log_line "Running from: $PWD"
log_line "slurm job id: $i"
log_line "Reference: ${REF}"
log_line "output directory: ${OUTDIR}"
log_line "input VCF file: ${KASP_VCF}"
log_line "awk started....."

#===================================================================================
#+++++++++ remove any SNPs that have a nother SNP 50bp up or dowstream of it
#===================================================================================

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${KASP_VCF} > 5B_1_100Mb.kasp.snps.final.tsv


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
}' 5B_1_100Mb.kasp.snps.final.tsv > 5B_1_100Mb.kaspSNPs.no_nearby_50bp.tsv


awk '{print $1"\t"$2}' 5B_1_100Mb.kaspSNPs.no_nearby_50bp.tsv > 5B_1_100Mb.kaspSNPsKeep_sites.tsv

# then using bcf tool to get the final set of SNPs

bcftools view -T 5B_1_100Mb.kaspSNPsKeep_sites.tsv ${KASP_VCF} -Oz -o 5B_1_100Mb.kaspSNPsKeep_sites.vcf.gz
bcftools index 5B_1_100Mb.kaspSNPsKeep_sites.vcf.gz


## convert SNPs to bed intervals so we cn extract flanking sequences

bcftools query -f '%CHROM\t%POS\n' 5B_1_100Mb.kaspSNPsKeep_sites.vcf.gz > 5B_1_100Mb.kaspSNPsKeep_sites.pos


# extract flanks, 50 bp for primer design, 100 bp for blast

awk 'BEGIN{OFS="\t"} {start=$2-51; end=$2+50; if(start<0) start=0; print $1, start, end, $1"_"$2}' 5B_1_100Mb.kaspSNPsKeep_sites.pos > 5B_1_100Mb.kaspSNPsKeep_sites50bpFlanks.bed

awk 'BEGIN{OFS="\t"} {start=$2-101; end=$2+100; if(start<0) start=0; print $1, start, end, $1"_"$2}' 5B_1_100Mb.kaspSNPsKeep_sites.pos > 5B_1_100Mb.kaspSNPsKeep_sites100bpFlanks.bed

# extract with bedtools 

bedtools getfasta -fi ${REF} -bed 5B_1_100Mb.kaspSNPsKeep_sites100bpFlanks.bed -fo 5B_1_100Mb.kaspSNPsKeep_sites100bpFlanks.fasta -name

bedtools getfasta -fi ${REF} -bed 5B_1_100Mb.kaspSNPsKeep_sites50bpFlanks.bed -fo 5B_1_100Mb.kaspSNPsKeep_sites50bpFlanks.fasta -name

# run blast 

blastn -num_threads ${slurmCPUs}  \
 -query 5B_1_100Mb.kaspSNPsKeep_sites100bpFlanks.fasta \
 -db ${REF} \
 -task blastn-short \
 -max_hsps 10 \
 -max_target_seqs 10 \
 -outfmt 6 \
 -out 5B_1_100Mb.kaspSNPsKeep_sites100bpFlanks_v2.txt



# to filter, we SNP queries with non-slf-hits were discarded
## after filtering the SNPs we want to keep after the blast alignment
## generate a file for feeding into polymarker with a 50 bp on each side of the snp


## first generate a file with the REF/SNP values

awk 'NR==FNR {keep[$1 FS $2]; next} ($1 FS $2) in keep' 5B_8-15Mb_blast_hitsKeepSitesFiltered.txt 5B_1_100Mb.kasp.snps.final.tsv >> 5B_8-15Mb_KeepSitesBlastFiltSNPs.tsv


# make flank beds
awk 'BEGIN{OFS="\t"}{
    id=$1"_"$2;
    s=$2-51; if(s<0) s=0;
    e=$2-1;
    print $1,s,e,id,$3,$4
}' 5B_8-15Mb_KeepSitesBlastFiltSNPs.tsv > kasp.left.bed

awk 'BEGIN{OFS="\t"}{
    id=$1"_"$2;
    s=$2;
    e=$2+50;
    print $1,s,e,id,$3,$4
}' 5B_8-15Mb_KeepSitesBlastFiltSNPs.tsv > kasp.right.bed

# extract flanks
bedtools getfasta -fi ${REF} -bed kasp.left.bed -fo kasp.left.fa -nameOnly
bedtools getfasta -fi ${REF} -bed kasp.right.bed -fo kasp.right.fa -nameOnly

# combine into bracketed fasta
paste \
  <(awk '/^>/{if(seq){print name"\t"seq}; name=substr($0,2); seq=""; next} {seq=seq$0} END{print name"\t"seq}' kasp.left.fa) \
  <(awk '/^>/{if(seq){print name"\t"seq}; name=substr($0,2); seq=""; next} {seq=seq$0} END{print name"\t"seq}' kasp.right.fa) \
  <(awk 'BEGIN{OFS="\t"}{print $1"_"$2,$1,$2,$3,$4}' 5B_8-15Mb_KeepSitesBlastFiltSNPs.tsv) \
| awk 'BEGIN{OFS="\n"}
{
    left_id=$1; left_seq=$2;
    right_id=$3; right_seq=$4;
    snp_id=$5; chrom=$6; pos=$7; ref=$8; alt=$9;

    if(left_id != snp_id || right_id != snp_id){
        print "ID mismatch: " left_id " " right_id " " snp_id > "/dev/stderr";
        next
    }

    print ">" snp_id;
    print left_seq "[" ref "/" alt "]" right_seq;
}' > 5B_8-15Mb_KeepSitesBlastFiltKaspSNPs_bracketed.fasta

# fasta -> csv
awk '
BEGIN{
    OFS=","
    print "marker,chr,sequence"
}
/^>/{
    marker=substr($0,2)
    split(marker,a,"_")
    chr=a[1]
    next
}
{
    print marker, chr, $0
}' 5B_8-15Mb_KeepSitesBlastFiltKaspSNPs_bracketed.fasta > 5B_8-15Mb_KeepSitesBlastFiltKaspSNPs_bracketed.csv


## now to subset the final set of SNPs get the samples and allele calls


bcftools view -R 5B_8-15Mb_PolymarkerSelectedSNPpos.txt 5B_1_100Mb.kasp.final.vcf.gz -Oz -o 5B_8-15Mb_PolymarkerSelectedSNPs.vcf.gz

# index
bcftools index 5B_8-15Mb_PolymarkerSelectedSNPs.vcf.gz


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%QUAL[\t%SAMPLE=%GT]\n' 5B_8-15Mb_PolymarkerSelectedSNPs.vcf.gz | \
awk 'BEGIN{
    OFS="\t";
    print "chrom","snp_pos","ref","alt","vcf_id","qual","n_called","n_alt_carriers","n_hom_alt","n_het","alt_carriers"
}
{
    chrom=$1;
    pos=$2;
    ref=$3;
    alt=$4;
    vid=$5;
    qual=$6;

    n_called=0;
    n_alt_carriers=0;
    n_hom_alt=0;
    n_het=0;
    carriers="";

    for(i=7; i<=NF; i++){
        split($i, a, "=");
        sample=a[1];
        gt=a[2];

        if(gt != "./." && gt != ".|."){
            n_called++;
        }

        if(gt == "1/1" || gt == "1|1"){
            n_alt_carriers++;
            n_hom_alt++;
            if(carriers=="") carriers=sample;
            else carriers=carriers ";" sample;
        }
        else if(gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0"){
            n_alt_carriers++;
            n_het++;
            if(carriers=="") carriers=sample;
            else carriers=carriers ";" sample;
        }
    }

    print chrom,pos,ref,alt,vid,qual,n_called,n_alt_carriers,n_hom_alt,n_het,carriers
}' > 5B_8-15Mb_PolymarkerSelectedSNPs.genotypes.summary.tsv