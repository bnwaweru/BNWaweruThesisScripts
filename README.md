# BNWaweru, PhD thesis scripts

Analysis code supporting my PhD thesis on pangenome-informed haplotype analysis for disease
resistance in wheat.

## Overview

Breeding for disease resistance depends on knowing which haplotype an accession carries at a
resistance locus. The difficulty is that most cultivars and landraces of interest have no genome
assembly, so the sequence underlying their haplotype is not directly available, and without that
sequence, diagnostic markers cannot be designed. This thesis combines the wheat pangenome with
*k*-mer-based methods to characterise haplotypes at resistance loci without needing an assembly for
every accession.

Chapter 2 surveys haplotype structure at known resistance loci using identity-by-state variation
counts. Chapter 3 develops a *k*-mer framework that reconstructs local haplotype sequence from reads
retrieved against an assembled anchor genome and converts it into candidate markers. Chapter 4
characterises structural variation at the *Yr47* locus across pangenome assemblies and Watkins
landrace genomes. Chapter 5 tests the genetics of *Yr47* resistance in a biparental mapping
population.

## Using this repository

Scripts are grouped by thesis chapter. Each is the code as it was run rather than a packaged tool:
file paths, module names and job-scheduler directives are specific to the John Innes Centre HPC,
and the input data which include genome assemblies, raw reads, precomputed *k*-mer databases, are held in
institutional storage rather than here.

## Pipelines

| | |
|---|---|
| [`KHaploTrackSingleAnchor/`](KHaploTrackSingleAnchor/) | The *k*-mer haplotype pipeline, single-anchor module 

## Scripts by chapter

### Chapter 2: Haplotype structure of known resistance loci using identity by state in Python (IBSpy)

- [`Chapter2_dataExtractionSyntenicWindowsYr7LinePlot.Rmd`](Chapter2_dataExtractionSyntenicWindowsYr7LinePlot.Rmd)
- [`Chapter2_syntenicWindowsNLRsandnonNLRBlocls.qmd`](Chapter2_syntenicWindowsNLRsandnonNLRBlocls.qmd)

### Chapter 3: A *k*-mer-based framework for haplotype reconstruction and marker discovery

- [`Chapter3_IBSpyHelperPipeline.ipynb`](Chapter3_IBSpyHelperPipeline.ipynb)
- [`Chapter3_dotPlotsNucmerAlign_v2.rmd`](Chapter3_dotPlotsNucmerAlign_v2.rmd)
- [`Chapter3_extractSyntenicWindowsPangenomes_Aug172026.qmd`](Chapter3_extractSyntenicWindowsPangenomes_Aug172026.qmd)
- [`Chapter3_mummerAlignChr2DPangenomes.sh`](Chapter3_mummerAlignChr2DPangenomes.sh)
- [`Chapter3_mummerDotplotFunctions.R`](Chapter3_mummerDotplotFunctions.R)

### Chapter 4: Structural variation at the *Yr47* locus in wheat pangenomes and representative resistant and susceptible Watkins landrace genomes

- [`Chapter4_NLRAnnotatorThesis.sh`](Chapter4_NLRAnnotatorThesis.sh)
- [`Chapter4_addGeneCoordsToNLRTrackerExtractNBARCSeqs.ipynb`](Chapter4_addGeneCoordsToNLRTrackerExtractNBARCSeqs.ipynb)
- [`Chapter4_bedtoolsIntersect.sh`](Chapter4_bedtoolsIntersect.sh)
- [`Chapter4_blastPangenomesToAnchorGeneSeqsThesis.sh`](Chapter4_blastPangenomesToAnchorGeneSeqsThesis.sh)
- [`Chapter4_customRiparianPlotThesis.Rmd`](Chapter4_customRiparianPlotThesis.Rmd)
- [`Chapter4_dotPlotFiguresWatkinsPGsChr5B.rmd`](Chapter4_dotPlotFiguresWatkinsPGsChr5B.rmd)
- [`Chapter4_nucmer5BPangenomesThesis.sh`](Chapter4_nucmer5BPangenomesThesis.sh)
- [`Chapter4_pafGenSimilarityMatAllvsAllGenes.ipynb`](Chapter4_pafGenSimilarityMatAllvsAllGenes.ipynb)
- [`Chapter4_plot_yr47_orthogroup_tree.py`](Chapter4_plot_yr47_orthogroup_tree.py)

### Chapter 5: Genetics of *Yr47* resistance using an F5:6 derived Paragon x WATDE0042 biparental mapping population

- [`Chapter5_F5_PxW308_RIL_DataAnalysis.Rmd`](Chapter5_F5_PxW308_RIL_DataAnalysis.Rmd)
- [`Chapter5_FilterSNPsBlast.ipynb`](Chapter5_FilterSNPsBlast.ipynb)
- [`Chapter5_bwaMEMAlign.sh`](Chapter5_bwaMEMAlign.sh)
- [`Chapter5_fastqcTrim.sh`](Chapter5_fastqcTrim.sh)
- [`Chapter5_filterBlastKASP-SNPs.py`](Chapter5_filterBlastKASP-SNPs.py)
- [`Chapter5_filterKaspSNPsRunBlast.sh`](Chapter5_filterKaspSNPsRunBlast.sh)
- [`Chapter5_freebayesCallVariants.sh`](Chapter5_freebayesCallVariants.sh)
- [`Chapter5_starMapping.sh`](Chapter5_starMapping.sh)
- [`Chapter5_stringtie.sh`](Chapter5_stringtie.sh)

## Software

Scripts call KMC, samtools, SeqKit, SPAdes, Flye, BLAST+, minimap2, MUMmer, STAR, StringTie,
freebayes, bedtools, NLR-Annotator and [IBSpy](https://github.com/Uauy-Lab/IBSpy), with analysis
in Python 3 and R. Versions are recorded in the individual scripts and in each chapter's methods
section in the thesis.

## Citing

<!-- add the thesis citation, and a Zenodo DOI if you mint one from the v1.0-thesis tag -->
