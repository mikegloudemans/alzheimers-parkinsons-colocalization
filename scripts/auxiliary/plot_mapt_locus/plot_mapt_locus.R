### Auxiliary script to plot PD GWAS and Brain Cortex eQTLs together
### for visualization.
###
### At least one large region is missing; follow-up checks determined
### that this region contained no SNPs mapped in either the eQTL or
### GWAS files; this is likely an unmappable repeat region.
###

require(dplyr)

header = strsplit(system("zcat data/eqtls/gtex_v7/Brain_Cortex.allpairs.txt.gz | head -n 1", intern=TRUE), "\\t")[[1]]

tabix_stream = system("tabix data/eqtls/gtex_v7/Brain_Cortex.allpairs.txt.gz 17:42601944-45642200", intern=TRUE)
eqtls = read.table(text=paste(tabix_stream, collapse="\n"), col.names=header, stringsAsFactors = FALSE)
eqtls$chr = paste0("chr", eqtls$chr)

header = strsplit(system("zcat data/gwas/prepared/23andme_PD.txt.gz | head -n 1", intern=TRUE), "\\t")[[1]]
tabix_stream = system("tabix data/gwas/prepared/23andme_PD.txt.gz chr17:42601944-45642200", intern=TRUE)
gwas = read.table(text=paste(tabix_stream, collapse="\n"), col.names=header, stringsAsFactors = FALSE)


plot(gwas$snp_pos, -log10(gwas$pvalue), main="GWAS", pch=16)


combo = inner_join(gwas, eqtls, by=c("chr", "snp_pos"), suffix=c("_gwas", "_eqtl"))

# Loop through plotting one gene at a time
for (gene in unique(combo$gene))
{
	pdf(paste0("output/mapt/plots/full-locus/", gene, ".pdf"), height=7, width=8)
		subcombo = combo[combo$gene == gene,]
		#X11()
		layout(matrix(c(1,2), nrow=2))
		plot(subcombo$snp_pos, -log10(subcombo$pvalue_gwas), main="GWAS", pch=16)
		plot(subcombo$snp_pos, -log10(subcombo$pvalue_eqtl), main=paste(gene, "eQTL"), pch=16)
	dev.off()
}
