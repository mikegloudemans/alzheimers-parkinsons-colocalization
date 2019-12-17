# Alzheimers / Parkinsons Colocalization Analysis

Analysis performed by Mike Gloudemans

Updated 12/16/2019

## Summary

These project contains the scripts required to perform the colocalization analyses described in the paper.

The top-level project directory should contain the folders `bin`, `data`, `output`, `tmp`, and `scripts`. All scripts should be
run from this top-level directory, or they'll be unable to locate the required files.

### Components within `scripts` folder

#### `scripts/pipe.sh`

The full workflow for all colocalization-related analyses.

#### `scripts/auxiliary`

A few other miscellaneous scripts for analyses that may be of interest but were not explicitly described in the paper.

#### `scripts/config`

Config files for selection of SNP lists and for the colocalization tests themselves.

#### `scripts/enumerate_candidate_snps`

Scripts for assembling a broad list of candidate SNPs and loci for further follow-up using the other
data sources described in the paper.

Also, scripts for annotating these candidates with GWAS, colocalization, and ATAC-seq peak overlaps.

And, importantly, QC checks for the above steps!

##### `scripts/enumerate_candidate_snps/ld_annotations`

Selection of "LD buddies" near GWAS SNPs.

#### `scripts/munge_gwas`

Scripts for formatting GWAS files consistently.

#### `scripts/test_mapt_haplotypes`

Scripts for differential testing of expression between the two possible MAPT haplotypes.


## Required Tools

The following analyses were performed using other publicly available tools. To fully complete the
analyses, you will have to install or link these tools within the `bin` subfolder, or modify the
`pipe.sh` script to include the paths to the directories where these tools are installed.

* The tools for downloading and munging the publicly available GWAS summary statistics
are available at https://github.com/mikegloudemans/gwas-download/.
* The colocalization pipeline is available in a basic form at https://bitbucket.org/mgloud/production_coloc_pipeline/src. 
An extended and hopefully easier-to-use pipeline with a greater variety of options and analyses, closer to what was used for this
paper, will soon be available at https://github.com/mikegloudemans/ensemble_coloc_pipeline.
* Graphical visualization of colocalizations was performed using [LocusCompare](https://locuscompare.com) (Liu et al. 2019)
* The analysis performed in this paper uses an integration of the publicly available tools [FINEMAP](http://www.christianbenner.com/) (Benner et al. 2016)
and [eCAVIAR](http://zarlab.cs.ucla.edu/tag/ecaviar/) (Hormozdiari et al. 2016)

## Required Data

This project makes uses of a variety of data tables, some already publicly available
and some internally generated. I'm currently exploring ways to just link all or most of the required
data files here as a single download. Until then, please contact me directly (see _Contact_ section below) and I'll share 
the relevant files directly, ASAP.

90% of the analysis, including the core colocalization analysis can be completed using just the following
data files:

### Getting started

* An hg38-formatted version of the [1000 Genomes VCF](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/)
is required for computing allele frequencies in a reference population.
* GWAS summary statistics are publicly available; consistently-formatted versions of these and other GWAS can be [downloaded directly](https://github.com/mikegloudemans/gwas-download).
* GTEx v8 brain eQTL association statistics can be downloaded from the [GTEx Portal](https://gtexportal.org/home/datasets). Some minor pre-processing will be required to run
these scripts; this process is described [here](https://bitbucket.org/mgloud/production_coloc_pipeline/src).

### All required data

To run all the scripts listed here, the `data` folder technically needs to contain all of the following files:

* `data/1KG`: 1KG VCF for hg38, publicly available for download as described above.
* `data/atacseq`: Single-cell and bulk ATAC-seq peaks. Availability TBA.
* `data/dbsnp`: Downloadable from dbSNP; versions for hg19 and hg38 will be needed. We recommend version 150 because this is what we used.
* `data/eqtls`: GTEx eQTLs for v8, downloadable from GTEx Portal as described above.
* `data/gencode`: Gencode annotations 
* `data/gtex`: GTEx folder containing raw RNA-seq counts, for the MAPT analysis
* `data/gwas`: All GWAS, downloadable as described above. If already re-formatted for colocalization, place them in a `formatted` subdirectory; if not, place them in a
   `raw` subdirectory.
* `data/gwas_top_hits`:	
* `data/hgnc_conversion`: Mapping of Ensembl gene IDs to HGNC names, obtained through Biomart.
* `data/indexed-dbsnp`: Sorted and tabixed version of dbsnp for fast access
* `data/ld`: Pairwise LD scores from 1K Genomes Phase 1, used for selecting candidate SNP lists.
* `data/mapt`: List of SNPs that almost always differ between the two possible MAPT haplotypes, used for labeling
GTEx samples by haplotype status

## Contact

For any questions about this pipeline, obtaining specific data files, or any colocalization-related analyses for this project, 
please contact Mike Gloudemans (mgloud@stanford.edu). I'll be glad to help you get these analyses up and running!
