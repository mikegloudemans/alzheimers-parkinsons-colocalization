#!/bin/bash

# DANGER DANGER DANGER DANGER DANGER DANGER DANGER DANGER DANGER DANGER
# WARNING: Ref and alt alleles are not correct here!
# 23andme specifies alleles but these alleles are just in
# alphabetical order. I assume that the ref / alt for these
# SNPs are actually just the ones designated in dbSNP (not
# sure which version)
# DANGER DANGER DANGER DANGER DANGER DANGER DANGER DANGER DANGER DANGER

echo -e "rsid\tchr\tsnp_pos\tref\talt\tpvalue\tbeta\tse\tpass" > data/gwas/prepared/23andme_PD_hg19.txt
join -t $'\t' <(tail -n +2 data/gwas/raw/23andme/snp_anno/all_snp_info-5.0.txt | cut -f1,4,5,6,7 | sed s/\\//\\t/g) <(tail -n +2 data/gwas/raw/23andme/sum_stats/PD_all_xmeta.dat | cut -f1,3,4,5,6) | cut -f1 --complement | awk '{if ($9 == "Y") print $0}' | sort -k2,2 -k3,3n >> data/gwas/prepared/23andme_PD_hg19.txt

bgzip -f data/gwas/prepared/23andme_PD_hg19.txt
tabix -f -s 2 -b 3 -e 3 -S 1 data/gwas/prepared/23andme_PD_hg19.txt.gz

# Do it for hg38
echo -e "rsid\tchr\tsnp_pos\tchr_hg19\tsnp_pos_hg38\tref\talt\tpvalue\tbeta\tse\tpass" > data/gwas/prepared/23andme_PD_hg38.txt
join -1 3 -2 1 <(zcat data/dbsnp/sorted_1kg_matched_hg38_snp150.txt.gz | sort -k3,3) <(zcat data/gwas/prepared/23andme_PD_hg19.txt.gz | tail -n +2 | sort -k1,1) | sort -k2,2 -k3,3n | sed s/\ /\\t/g >> data/gwas/prepared/23andme_PD_hg38.txt

bgzip -f data/gwas/prepared/23andme_PD_hg38.txt
tabix -f -s 2 -b 3 -e 3 -S 1 data/gwas/prepared/23andme_PD_hg38.txt.gz
