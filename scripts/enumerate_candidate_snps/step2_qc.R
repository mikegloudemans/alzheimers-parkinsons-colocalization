# QC for GWAS and effect direction annotation step
#
# Load step 2 data
require(tidyverse)
data = read.table("output/candidate_snp_lists/ld_buddies_table_stage2.tsv", header=TRUE, stringsAsFactors = FALSE)

# All SNPs selected for “GWAS” should have GWAS p < 5e-8.
range(data[grepl("GWAS", data$source),]$pvalue)

# All SNPs not tagged as “GWAS” should either belong to Nalls-Chang group, be in LD or should have p > 5e-8.
bool = !grepl("GWAS", data$source) & data$pvalue < 5e-8
bool[is.na(bool)] = FALSE
table(data[bool,]$source_gwas)

# All SNPs that do not appear in the GWAS at all should be in the LD category or Nalls-Chang
table(data[is.na(data$pvalue),]$source)
# Well, it turns out there are a few coloc SNPs that are not in the GWAS at all; this is because
# they were selected on the basis of being significant in a different trait. Weirder still, for
# some of them this means the "seed" coloc SNP wasn't even tested in that particular GWAS!
#
# We'll test this instead during the final step below, the manual visualization of colocalization.

# All SNPs tagged “LD” should be tagging a SNP from the same GWAS that was selected for one of the above reasons.
all_ld_tags = unique(paste(data$ld_tag_locus, data$source_gwas, sep="_"))
table(data[paste(data$snp_id, data$source_gwas, sep="_") %in% all_ld_tags,]$source)

# Direction should be roughly 50/50 between positive and negative for each GWAS
table(data[c("direction", "source_gwas")])

# And directions of SNPs in multiple GWAS should be consistent
# Just for QC, swap the ones with negative direction for
# comparability
data$tmp = ""
data$effect_allele = toupper(data$effect_allele)
data$noneffect_allele = toupper(data$noneffect_allele)
neg_snps = (data$direction == "-")
neg_snps[is.na(neg_snps)] = FALSE
data$tmp[neg_snps] = data$effect_allele[neg_snps]
data$effect_allele[neg_snps] = data$noneffect_allele[neg_snps]
data$noneffect_allele[neg_snps] = data$tmp[neg_snps]
mismatch_detection = data %>% group_by(snp_id) %>% summarize(types = length(unique(paste(effect_allele, noneffect_allele))))
table(mismatch_detection$types) # Should be 1 for all SNPs


# TODO:
# All coloc SNPs should verifiably be actually colocalized, when looking at the associated GWAS and eQTL plots
colocs = data[grepl("coloc", data$source),]
colocs = colocs[order(colocs$source_gwas),]
print(paste(colocs$snp_id, colocs$source_gwas))


