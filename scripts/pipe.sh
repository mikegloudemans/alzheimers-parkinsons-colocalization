# Author: Mike Gloudemans
# Full pipeline for AD/PD GWAS analysis

##########################################################
# Colocalization testing and functional annotation
##########################################################

mkdir data
mkdir output
mkdir tmp
mkdir bin

mkdir -p output/ld_buddies
mkdir -p output/mapt/plots
mkdir -p output/auxiliary
mkdir -p output/test-snps
mkdir -p output/ld_buddies
mkdir -p output/candidate_snp_lists
mkdir -p output/colocalization
mkdir -p output/coloc_summary

mkdir -p output/mapt/plots/full-locus
mkdir -p output/mapt/plots/boxplots
mkdir -p output/mapt/plots/linearplots

#
# Preprocessing
#

# Munge all GWAS files
bash scripts/munge_gwas/munge_pipe.sh

#
# Colocalization
#

# Get list of all SNPs to test for colocalization
python bin/overlap/list_snps_to_test.py scripts/config/overlap/ad-pd.overlap.config

snp_list=output/test-snps/ad-pd-v8_all-gwas_gtex-brain_gwas-pval1e-05_eqtl-pval1e-05_gwas-window500000_eqtl-window0_coloc-tests.txt

cut -f1,2,8 $snp_list | tail -n +2 | sort | uniq > output/test-snps/final-test-list.txt

# Run colocalization on all of these SNPs
python bin/ensemble_coloc_pipeline/dispatch.py scripts/config/colocalization/ad-pd.coloc.config 20

#
# Get LD buddies per GWAS (for significant hits and for colocalized hits)
#


# First generate the list of ALL SNPs passing the (loose) threshold, even if not lead SNP
# Will filter these later for the final list
python bin/gwas-download/overlap/list_snps_to_test.py scripts/config/overlap/ad-pd-all-passing.overlap.config

# Split it up by trait
for g in `tail -n +2 output/test-snps/ad-pd-v8-all-passing_all-gwas_gwas-pval1e-05_gwas-window0_snps-considered.txt | cut -f4 | sort | uniq`; do
	gwas=`basename $g`
	head -n 1 output/test-snps/ad-pd-v8-all-passing_all-gwas_gwas-pval1e-05_gwas-window0_snps-considered.txt > output/test-snps/$gwas.tests.txt
	grep $gwas output/test-snps/ad-pd-v8-all-passing_all-gwas_gwas-pval1e-05_gwas-window0_snps-considered.txt > output/test-snps/$gwas.tests.txt
done

# Now need to get LD buddies for all of them
# Convert hg19 LD matrix to hg38 first
python scripts/enumerate_candidate_snps/ld_annotations/lift_ld_maps.py

# Then get LD buddies for these SNPs (not) all of these will be included for 
# additional analysis, but we'll filter the list down afterwards
python scripts/enumerate_candidate_snps/ld_annotations/get_ld_buddies.py

# Annotate SNPs on how they were selected and other metadata 
# Perform QC checks to make sure everything looks good
Rscript scripts/enumerate_candidate_snps/step1_assemble_ld_annotation_matrix.R
Rscript scripts/enumerate_candidate_snps/step1_qc_check_render.R
Rscript scripts/enumerate_candidate_snps/step2_gwas_coloc_addition.R
Rscript scripts/enumerate_candidate_snps/step2_qc.R
Rscript scripts/enumerate_candidate_snps/step3_add_atac_peaks.R

# Aggregate colocalization information by locus
Rscript scripts/post_coloc/aggregate_coloc_table.R

#####################################################
# MAPT analysis (independent of above analyses)
#####################################################


# Sort individuals by which haplotype they have
python scripts/test_mapt_haplotypes/test_mapt_haplotypes.py

# Test differences in expression between haplotypes
Rscript scripts/test_mapt_haplotypes/plot_mapt_comparisons.R

#####################################################
# Additional postprocessing scripts
#####################################################

# Count numbers of colocalizations and extent of overlap between
# different traits
bash scripts/auxiliary/analyze_results.sh > output/auxiliary/colocalization_quick_metrics.txt

