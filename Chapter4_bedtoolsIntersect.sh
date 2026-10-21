#!/bin/bash
#SBATCH --partition=jic-medium,nbi-medium,jic-long,nbi-long
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 10G
#SBATCH -t 10-12:00 # time (D-HH:MM)
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batchLogs/bedtools_\%A_\%a.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis/scripts/batchLogs/bedtools_\%A_\%a.err # STDERR
#SBATCH --job-name=bedtoolsInt
#SBATCH --array=0-1

i=${SLURM_ARRAY_TASK_ID:-0}

# initiate log reporting
function log_line() {
	echo $(date) "$1" >&2
}

##==== Objective =====

#-- Extract genes of WATDE0042/WATDE0100 predicted by helixer that ovelap with NLRs predicted by NLR annotator

## ==== load tools =============
#===========================================

source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  ## bedtools 2.31.0
source package 9117dc0d-5bb0-48ce-83c3-bc0db08e0cd3  # gffread 0.12.8
source package 1041444f-cd25-4107-a5c7-5e86cb1728fe  #  seqkit 2.8.2
source package 4b4051e0-bd0f-4c18-bd5e-962700bf4b28  ## Samtools 1.23
source package 4c883633-af2d-4fac-ab67-a1574f7fe079  ## Agat 0.9.2

## =======  CONFIG ==========

workDir="/jic/scratch/groups/Cristobal-Uauy/waweru/2025-09_5B_comparativeAnalysis"
cd ${workDir}

#watGenomeDir="/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies/2025-07-02_Watkins_assemblies"
pg2GenomeDir="/jic/scratch/groups/Cristobal-Uauy/waweru/assemblies"
helixergffs=${workDir}/helixerOut
NLRAnnotgffs=${workDir}/NLRAnnotatorOut

pgGff3s="/jic/research-groups/Cristobal-Uauy/assemblies/10wheat_annotation"
pgNLRs="/jic/scratch/groups/Cristobal-Uauy/waweru/10wheatgenomes_annotations_bnw/pangenome_references_NLR_annotations_Jan_2023"

bedtoolsOutDir=${workDir}/bedtoolsOut
nlrGenesFastas=${workDir}/bedtoolsFastaOut

mkdir -p ${bedtoolsOutDir} ${nlrGenesFastas}

#declare -a queries=("WATDE0042" "WATDE0100")
declare -a queries=(\
#"arinaLrFor" \
#"chinesespring" \
#"julius" \
#"jagger" \
#"lancer" \
#"landmark" \
#"norin61" \
#"mace" \
#"stanley" \
#"symattis" \
#"cadenza_v2" \
#"fielder_v1" \
#"paragon_v3" \
"renan_v2.1" \
#"kariega_v1" 
)


#query=$(($i%2))
query=${queries[$i]}


#----------------------------------------------------------------------------
## ------ Define input and output files based on query name variable --------
#----------------------------------------------------------------------------

### for the watkins assemblies
#++++++++++++++++++++++++++++++++++++

#genesGFF=${helixergffs}/${query}_chr5B_helixer_v2.gff3
#genesGFF=${helixergffs}/${query}_WG_helixer.gff3   ## also cadenza, fielder and paragon
genesGFF=${helixergffs}/${query}_chr.gff3 # for renan, and without the chr for kariega
#genesBED=${helixergffs}/${query}_chr5B_helixer.bed
genesBED=${helixergffs}/${query}_chr5B.bed
nlrsBED=${NLRAnnotgffs}/${query}/${query}_chr5B.bed
#pgFasta=${watGenomeDir}/IRN_${query}_v1.fasta
pgFasta=${pg2GenomeDir}/${query}_renamed.fasta


### for the pangenome assemblies
#++++++++++++++++++++++++++++++++++++
#genesGFF=${pgGff3s}/${query}.consolidated.Oct2020.v12.gff3
#genesBED=${helixergffs}/${query}_chr5B_genes.bed
#nlrsTxt=${pgNLRs}/${query}.nlr.txt
#nlrsBED=${NLRAnnotgffs}/${query}_chr5B_NLRs.bed
#pgFasta=${pgGff3s}/${query}.fa

genesNLRsBed=${bedtoolsOutDir}/${query}_chr5B_genes-NLRs_overlap.bed
getFastaBed=${bedtoolsOutDir}/${query}_chr5B_0-100Mb_genesNLRs_getFasta.bed
nlrGenesFasta=${nlrGenesFastas}/${query}_chr5B_0-100Mb_NLRsFullGene.fasta
nlrTargetGenesTxt=${bedtoolsOutDir}/${query}_chr5B_0-100Mb_NLRsgenes.txt
nlrTargetGenesGff3=${bedtoolsOutDir}/${query}_chr5B_0-100Mb_NLRsgenes.gff3
nlrGenesCDSFasta=${nlrGenesFastas}/${query}_chr5B_0-100Mb_NLRsCDS.fasta
nlrGenesCDSPep=${nlrGenesFastas}/${query}_chr5B_0-100Mb_NLRsCDSpep.fasta
nlrGenesCDSTrasc=${nlrGenesFastas}/${query}_chr5B_0-100Mb_NLRsCDStranscripts.fasta

mkdir -p ${bedtoolsOutDir}/geneSpaceFiles
geneSpaceBedFile1=${bedtoolsOutDir}/geneSpaceFiles/${query}_temp.bed
geneSpaceBedFile2=${bedtoolsOutDir}/geneSpaceFiles/${query}.bed
geneSpacePepFile1=${bedtoolsOutDir}/geneSpaceFiles/${query}_temp.fa
geneSpacePepFile2=${bedtoolsOutDir}/geneSpaceFiles/${query}_withDups.fa
geneSpacePepFile3=${bedtoolsOutDir}/geneSpaceFiles/${query}.fa 



## grep for genes on chr5B from the genes gff3 (chr,start,end, geneID, strand)
#======================================================================================

#grep -w "gene" ${genesGFF} | \
#awk 'BEGIN{OFS="\t"} { 
#    match($9, /Name=([^;]+)/, a); 
#    if (a[1] == "") match($9, /ID=([^;]+)/, a); 
#    print $1, $4-1, $5, a[1], $7
#}' > ${genesBED}


## using this instead so that we are ectracting chr5B genes from the whole genome gff3
awk 'BEGIN{OFS="\t"}
    $1 == "chr5B" && $3 == "gene" {
        match($9, /Name=([^;]+)/, a)

        if (a[1] == "") {
            match($9, /ID=([^;]+)/, a)
        }

        print $1, $4-1, $5, a[1], $7
    }
' "${genesGFF}" > "${genesBED}"


## grep for NLRs on chr5B and write out a bedfile with required fields (chr,start,end, name, completeness)
# this is only for the PG queries
#======================================================================================

## for pangenome NLR files 
#sed 's/partial (pseudogene)/partial(pseudogene)/g' ${nlrsTxt} | \
#    sed 's/complete (pseudogene)/complete(pseudogene)/g' | \
#    awk '$1 == "chr5B" {print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3}' > ${nlrsBED}


##------- Monitor variables -------------
#------------------------------------------
log_line "current working dir: $PWD"
log_line "slum job ID: ${i}"
log_line "current genes bed file: ${genesBED}"
log_line "current NLR Annot bed file: ${nlrsBED}"
log_line "output directory: ${bedtoolsOutDir}"
log_line "bedtools run starting ....."


#===== Check main input files ==================
#======================================================================================
if [[ ! -s "$genesBED" ]]; then
    log_line "ERROR: helixer genes bed file $genesBED missing or empty" >&2
    exit 1
fi
if [[ ! -s "$nlrsBED" ]]; then
    log_line "ERROR: NLR Annot bed file $nlrsBED missing or empty" >&2
    exit 1
fi


#===== run bedtools intersect ====================
#====================================================

## -f 1.0 , reports only 100% alignemnts with the gene, to report any overlap remove it, or add lower threshold (0.1 for 10%)
## -u reports each overlap only once
## wa -report entry in a
## -wb - report entry in B

bedtools intersect -wa -wb -a ${genesBED} -b ${nlrsBED} > ${genesNLRsBed}

# extract bed file that will be used to get the fasta, that will be strand aware and with gene being a combination of gene ID and NLR id
# (chr,start,end, name, strand)
# genes within between 0-100 Mb
#======================================================================================================================

#awk '{OFS="\t"; print $1, $2, $3, $4"|"$(NF-1), $5}' ${genesNLRsBed} > ${getFastaBed}
# only extract for the first 100 Mb, this we to extract the whole gene sequence
# $2 < 100000000 && $3 < 100000000 → filters for entries with both start and end < 100 Mb
# OFS="\t" → keeps output tab-separated

# for pangenomes
#awk '$2 < 100000000 && $3 < 100000000 {OFS="\t"; print $1, $2, $3, $4"|"$(NF-1), $5}' ${genesNLRsBed} > ${getFastaBed}

# for watkins, cadenza, paragon and fielder
awk '$2 < 100000000 && $3 < 100000000 {OFS="\t"; print $1, $2, $3, $4"|"$9, $5}' ${genesNLRsBed} > ${getFastaBed}

#============= Check bedtools intersect worked ====================

if [[ -s "${getFastaBed}" ]]; then
    log_line "bedtools overlap successfully completed"
else
    log_line "ERROR: bedtools failed, no output" >&2
    exit 1
fi


## Extract the full NLR gene fasta sequence, this will be used for gene to gene alignemnts
#==========================================================================================
# -name -add gene name
# -s - strand aware

bedtools getfasta -fi ${pgFasta} -bed ${getFastaBed} -fo ${nlrGenesFasta} -name -s

#============= Check bedtools getfasta output ====================

if [[ -s "${nlrGenesFasta}" ]]; then
    log_line "DONE: Fasta sequences successfully extracted for ${query} !!"
else
    log_line "ERROR: Fasta sequences extraction for ${query} failed, no output" >&2
    exit 1
fi


## ---------- Extract the CDS for the NLRs in the first 100 Mb of the chromosome ------
#=======================================================================================

## get the gene names of NLR genes (between 0-100 Mb)
#------------------------------------------------------ 

awk '$2 < 100000000 && $3 < 100000000 {OFS="\t"; print $4}' ${genesNLRsBed} > ${nlrTargetGenesTxt}

## get gff3 infor (all features) for the NLR genes only , listed in the target txt file
#---------------------------------------------------------------------------------------
# -f Obtain patterns from FILE, one per line.
# -F  Interpret PATTERNS as fixed strings, not regular expressions.

grep -F -f ${nlrTargetGenesTxt} ${genesGFF} > ${nlrTargetGenesGff3}

## get a bed file to use for GeneSpace from the final gff3 file

awk '$3=="gene" {
    split($9,a,"ID=");
    split(a[2],b,";");
    print $1, $4, $5, b[1];

}' OFS="\t" ${nlrTargetGenesGff3} > ${geneSpaceBedFile1}

## -------- run gffread ---------------
#---------------------------------------
## use gffread to extract the target gene CDS sequence from the genome file translate them into peptide then also get a tracripts file which is optional
#-g genome.fa → genome reference

## gffread options
# -x cds.fa → spliced CDS sequences (nucleotide)
# -y peptides.fa → translated peptides
# -w transcripts.fa → full transcript sequences (optional)
# -F → keep full model (includes subfeatures)
# -K → skip genes not found in GFF3
# -O   process other non-transcript GFF records (by default non-transcript records are ignored)
# --add-hasCDS : add a "hasCDS" attribute with value "true" for transcripts that have CDS features
# --cluster-only: same as -M/--merge but without discarding any of the "duplicate" transcripts, only create "locus" features

gffread -g ${pgFasta} \
    -x ${nlrGenesCDSFasta} \
    -y ${geneSpacePepFile1} \
    -w ${nlrGenesCDSTrasc} \
    -F --cluster-only \
    ${nlrTargetGenesGff3} 

   


## clean up the gene names for the genespace files


awk '{
  if ($0 ~ /^>/) {
    split($0,a," ");   # split by space
    id = a[1];         # take only the first field
    print id
  } else {
    print
  }
}' ${geneSpacePepFile1} > ${geneSpacePepFile2}

# Step 1: Select longest transcript per gene
awk -v sp="$query" '
  /^>/ {
    if(seq!=""){
      if(length(seq)>max_len[gene]){
        max_seq[gene]=seq
        max_len[gene]=length(seq)
        max_id[gene]=gene
      }
    }
    header=$0
    sub(/^>/,"",header)
    split(header,a,"[ .]")   # gene = first field before .1/.2
    gene=a[1]
    seq=""
    next
  }
  /^[^>]/ { seq=seq $0 }
  END {
    if(seq!=""){
      if(length(seq)>max_len[gene]){
        max_seq[gene]=seq
        max_len[gene]=length(seq)
        max_id[gene]=gene
      }
    }
    for(g in max_seq){
      print ">" sp "_" g
      print max_seq[g]
    }
  }
' "${geneSpacePepFile2}" > "${geneSpacePepFile3}"

log_line "[INFO] Longest isoform FASTA created: ${geneSpacePepFile3}"

# Step 2: Clean sequences (replace invalid chars)
sed -i '' 's/\./X/g; s/\*/X/g; s/[a-z]/\U&/g' "${geneSpacePepFile3}"

log_line "[INFO] Sequences cleaned for Diamond/OrthoFinder."

# Step 3: Update BED file with species prefix
awk -v sp="$query" 'BEGIN{FS=OFS="\t"} {$4=sp"_"$4}1' "${geneSpaceBedFile1}" > "${geneSpaceBedFile2}"

log_line "[INFO] BED file updated with species prefix: ${geneSpaceBedFile2}"

log_line "[DONE] Files ready for GeneSpace."


## then use seqkit to keep the longest gene transcript, that is remove the duplicate and just keep the longest one, the above one removes the .1 or .2 in the names hence the names will be duplicated if a gen has more than one transcript

seqkit rmdup -s ${geneSpacePepFile2} -o ${geneSpacePepFile3}


# Check gff3 output
if [[ -s "${geneSpacePepFile3}" ]]; then
    log_line "✅ DONE: CDS Fasta sequences successfully extracted and translated for ${query} !!"
else
    log_line "❌ ERROR: CDS Fasta sequences extraction and translation for ${query} failed, no output" >&2
    exit 1
fi


### for renan, we first converted the gft to gff3 uwing the tool agat

#agat_convert_sp_gxf2gxf.pl --gff renan_v2.1.gtf --gff_version_input 2.5 --gff_version_output 3 --output renan_v2.1.gff3 