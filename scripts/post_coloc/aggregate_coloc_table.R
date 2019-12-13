require(tidyverse)

# Make a table summarizing all critical info about colocalizations

ld_buddies_table = read.table("output/candidate_snp_lists/ld_buddies_table_stage3.tsv", header=TRUE)
ld_buddies_table$snp_id = as.character(ld_buddies_table$snp_id)

header_stream = gsub("-", "", system("cat output/colocalization/2019-08-24_19-43-32.358274_ad-pd/*clpp* | head -n 1", intern=TRUE))
coloc_stream = system("cat output/colocalization/2019-08-24_19-43-32.358274_ad-pd/*clpp* | sort -k6,6gr | uniq", intern=TRUE)
coloc_table = read_delim(file=paste(coloc_stream, collapse="\n"), delim="\t", col_names = strsplit(header_stream, "\t")[[1]])
coloc_table$short_ensembl = substring(coloc_table$feature, 1, 15)

genes = read.table("data/hgnc_conversion/mart_export.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
genes = genes %>% filter(HGNC.symbol != "") %>% select(-Transcript.stable.ID)
genes = genes[!duplicated(genes$Gene.stable.ID),]
coloc_table = left_join(coloc_table, genes, by=c("short_ensembl" = "Gene.stable.ID"))
coloc_table[is.na(coloc_table$HGNC.symbol),]$HGNC.symbol = coloc_table[is.na(coloc_table$HGNC.symbol),]$feature


coloc_table$base_gwas_file = gsub("_hg38_txt_gz", "", coloc_table$base_gwas_file)
coloc_table$base_gwas_file = gsub("_2019_txt_gz", "", coloc_table$base_gwas_file)
coloc_table$base_gwas_file = gsub("_2013_txt_gz", "", coloc_table$base_gwas_file)
coloc_table$base_gwas_file = gsub("_2012_txt_gz", "", coloc_table$base_gwas_file)
coloc_table$base_gwas_file = gsub("_2018_txt_gz", "", coloc_table$base_gwas_file)

ld_buddies_table$source_gwas = gsub("_2018", "", ld_buddies_table$source_gwas)
ld_buddies_table$source_gwas = gsub("_2012", "", ld_buddies_table$source_gwas)
ld_buddies_table$source_gwas = gsub("_2019", "", ld_buddies_table$source_gwas)
ld_buddies_table$source_gwas = gsub("_2013", "", ld_buddies_table$source_gwas)
ld_buddies_table$source_gwas = gsub("_hg38", "", ld_buddies_table$source_gwas)

coloc_table = coloc_table %>% filter(clpp > 0.01, n_snps >= 20, log_gwas_pval > 4, log_eqtl_pval > 4)
coloc_table = coloc_table[!grepl("Phosph", coloc_table$base_gwas_file),]
coloc_table = coloc_table[!grepl("Amyloid", coloc_table$base_gwas_file),]
coloc_table = coloc_table %>% arrange(-clpp)

# Map each ref_snp to its corresponding locus
locus_codes = ld_buddies_table[c("snp_id", "locus_num")]
locus_codes = locus_codes[!duplicated(locus_codes),]
coloc_table = left_join(coloc_table, locus_codes, by=c("ref_snp" = "snp_id"))

sqtl_traits = function(x) {temp = paste(x[grepl("clu", x)], collapse=";")}
eqtl_traits = function(x) {temp = paste(x[grepl("ENS", x)], collapse=";")}
num_eqtl_traits = function(x) {temp = length(unique(x[grepl("ENS", x)]))}
sqtl_tissues = function(x) {temp = paste(gsub("Brain_", "", gsub("_sQTLs_txt_gz", "", x[grepl("sQTL", x)])), collapse=";")}
eqtl_tissues = function(x) {temp = paste(gsub("Brain_", "",gsub("_allpairs_txt_gz_eQTLs_txt_gz", "", x[grepl("eQTL", x)])), collapse=";")}
coloc_stats = coloc_table %>% group_by(base_gwas_file, locus_num) %>% summarize(has_eqtl_coloc = sum(grepl("eQTL", eqtl_file)) > 0, has_sqtl_coloc = sum(grepl("sQTL", eqtl_file)) > 0, sqtl_traits = sqtl_traits(feature), eqtl_traits=eqtl_traits(feature), eqtl_tissues = eqtl_tissues(eqtl_file), sqtl_tissues = sqtl_tissues(eqtl_file), num_eqtl_traits=num_eqtl_traits(feature))

top_colocs = coloc_table %>% arrange(-clpp) %>% filter(grepl("eQTL", eqtl_file))
top_colocs = top_colocs[!duplicated(top_colocs[c("ref_snp", "base_gwas_file")]),]

coloc_stats = left_join(coloc_stats, top_colocs)

coloc_stats = coloc_stats %>% rename(top_coloc_tissue = eqtl_file) %>% select(-n_snps, -gwas_trait)
coloc_stats$trait = ""
coloc_stats$trait[grepl("Alz", coloc_stats$base_gwas_file)] = "Alzheimers"
coloc_stats$trait[grepl("23andme", coloc_stats$base_gwas_file)] = "Parkinsons"
coloc_stats$trait[grepl("Pank", coloc_stats$base_gwas_file)] = "Parkinsons"

coloc_stats$top_coloc_tissue = gsub("_allpairs_txt_gz_eQTLs_txt_gz", "", coloc_stats$top_coloc_tissue)
coloc_stats = coloc_stats %>% arrange(locus_num, clpp)

atac_overlaps_per_locus = ld_buddies_table %>% 
	group_by(source_gwas, locus_num) %>% 
	summarize(number_ld_buddies = max(number_ld_buddies),
		  num_ATAC_overlaps = sum(direct_atac_overlap_single_cell | direct_atac_overlap_broad_tissue_regions | direct_atac_overlap_narrow_tissue_regions), 
		  num_single_cell_overlaps = sum(direct_atac_overlap_single_cell), 
		  num_broad_tissue_overlaps = sum(direct_atac_overlap_broad_tissue_regions), 
		  num_narrow_tissue_overlaps = sum(direct_atac_overlap_narrow_tissue_regions))

coloc_stats = left_join(coloc_stats, atac_overlaps_per_locus, by = c("base_gwas_file" = "source_gwas", "locus_num"))

# Right now, the column order isn't the best AND
# some of the loci appear doubly because they have more than one ref_snp
# from different GWAS lead hits...
#
# Could fix that if it's an obstruction
#
# All in all though, it looks not too bad

write_delim(coloc_stats, path="output/coloc_summary/coloc_details_table1.tsv", quote_escape=FALSE, delim="\t")
