require(dplyr)
require(readr)

########################################################
# Constants
########################################################

# Boundaries of MAPT region
mapt_start = 45510000
mapt_end = 46580000

# If log2 expression is below this value, we consider the gene unexpressed
expression_cutoff = 1

# Plotting options

plot_boxplot_pdfs = FALSE

########################################################
# Loading data
########################################################

# Load data matrix
data = read_delim("output/mapt/mapt_expression_matrix.txt", delim="\t")

data$log2_expression = log(data$expression + 1, 2)
data$gene = substring(data$gene, 1, 15)

# Get HGNC names and add them
genes = read.table("data/hgnc_conversion/mart_export.txt", header=TRUE, sep="\t")
genes = genes[,-2]
names(genes) = c("ensembl", "hgnc")
genes = genes[!duplicated(genes$ensembl),]
data = merge(genes, data, all.y=TRUE, by.x="ensembl", by.y="gene")
colnames(data)[1] = "gene"

data = data[data$haplotype %in% c("H1/H1", "H2/H1", "H2/H2"),]

########################################################
# Gene-level tests and box plots
########################################################

# For starters, just loop through each gene / tissue combo.
# Do a basic box plot and t-test

for (brain_region in unique(data$brain_region))
{
	print(brain_region)
	for (gene in unique(data$gene))
	{
		sub = data[(data$gene == gene) & (data$brain_region == brain_region),]
		het = sub[sub$haplotype=="H2/H1",]
		homo_h1 = sub[sub$haplotype=="H1/H1",]
		w = wilcox.test(het$log2_expression, homo_h1$log2_expression)
		if (is.na(w$p.value))
		{
			w$p.value = 1
		}

		# Plot boxplots only if requested in "constants" section
		if (plot_boxplot_pdfs)
		{
			pdf(paste0("output/mapt/plots/boxplots/", gsub(" ", "", brain_region), gene, ".pdf"), width=5, height=8)
				boxplot(log2_expression ~ haplotype, data = sub, main = paste(brain_region, "\n", gene, "\n", paste0("Wilcoxon p-value (H1/H1 vs. H2/H1): ", formatC(w$p.value, format="e", digits=2))))
			dev.off()
		}
	}
}

########################################################
# Plot differential expression across entire MAPT locus
########################################################
for (brain_region in unique(data$brain_region))
{
	bars = data.frame(list(diff=rep(0,length(unique(data$gene))),
			       pos=rep(0,length(unique(data$gene))),
			       start=rep(0,length(unique(data$gene))),
			       end=rep(0,length(unique(data$gene))),
			       homoh1=rep(0,length(unique(data$gene))),
			       het=rep(0,length(unique(data$gene))),
			       homoh2=rep(0,length(unique(data$gene))),
			       median=rep(0,length(unique(data$gene))),
			       signif=rep("",length(unique(data$gene))),
			       gene=rep("",length(unique(data$gene))),
			       gene_hgnc=rep("",length(unique(data$gene))),
			       var=rep(0,length(unique(data$gene)))), stringsAsFactors=FALSE)
	for (i in 1:length(unique(data$gene)))
	{
		gene = unique(data$gene)[i]
		sub = data[(data$gene == gene) & (data$brain_region == brain_region),]

		het = sub[sub$haplotype=="H2/H1",]
		homo_h1 = sub[sub$haplotype=="H1/H1",]
		homo_h2 = sub[sub$haplotype=="H2/H2",]

		w = wilcox.test(homo_h1$log2_expression, het$log2_expression)
		if (is.na(w$p.value))
		{
			w$p.value = 1
		}
		if ( w$p.value < (0.05 / length(unique(data$gene)) / length(unique(data$brain_region))) )
		{
			bars$signif[i] = "*"
		}

		bars$gene[i] = as.character(gene)
		bars$gene_hgnc[i] = as.character(sub$hgnc[1]) 
		bars$diff[i] = median(het$log2_expression) - median(homo_h1$log2_expression)
		bars$start[i] = het$start[1]
		bars$end[i] = het$end[1]
		bars$homoh1[i] = median(homo_h1$log2_expression)
		bars$het[i] = median(het$log2_expression)
		bars$homoh2[i] = median(homo_h2$log2_expression)
		bars$median[i] = median(c(homo_h1$log2_expression, het$log2_expression, homo_h2$log2_expression))
		bars$var[i] = var(c(homo_h1$log2_expression, het$log2_expression, homo_h2$log2_expression))
	}
	bars$gene_hgnc[bars$gene_hgnc == ""] = bars$gene[bars$gene_hgnc == ""]			# If no HGNC name, just use Ensembl ID
	bars$gene_hgnc[(abs(bars$diff) < 0.3) & (!(grepl("MAPT", bars$gene_hgnc)))] = ""	# Only show genes if large difference or if MAPT
	
	# Filter out if hardly any expression?
	#bars = bars[abs(bars$median) > 1,]
	color <- rgb(190, 190, 190, alpha=80, maxColorValue=255)

	pdf(paste0("output/mapt/plots/linearplots/", gsub(" ", "", brain_region), ".pdf"), width=15, height=8)
		plot(bars$start, bars$diff, type="h", ylab="Log-fold change, H1/H2 minus H1/H1", xlab="Position of TSS on chromosome 17", main = paste("Expression effects in", brain_region), lwd=5, col="gray70", ylim=c(-1.5, 2.2), xlim=c(mapt_start, mapt_end))
		abline(a=0, b=0, col="red")
		text(bars$start, bars$diff, bars$gene_hgnc, cex=0.8, pos=2, col="black")
		text(bars$start[bars$diff < 0], bars$diff[bars$diff < 0], bars$signif[bars$diff < 0], cex=2, pos=1, col="black")
		text(bars$start[bars$diff > 0], bars$diff[bars$diff > 0], bars$signif[bars$diff > 0], cex=2, pos=3, col="black")
	dev.off()
}


########################################################
# Additional calculations for the paper
########################################################
cortex_only = data[data$brain_region == "Brain - Cortex",]
mapt = cortex_only[cortex_only$hgnc == "MAPT",]
kansl1_as = cortex_only[cortex_only$hgnc == "KANSL1-AS1",]
mapk8 = cortex_only[cortex_only$hgnc == "MAPK8IP1P2",]

# Narrow it down to the actual region we're most interested in
cortex_only = cortex_only[pmax(cortex_only$start, cortex_only$end) > mapt_start,]
cortex_only = cortex_only[pmin(cortex_only$start, cortex_only$end) < mapt_end,]
mean_expression = cortex_only %>% group_by(gene) %>% summarize(mean_expression = mean(log2_expression), hgnc=hgnc[1])
print(mean_expression[mean_expression$mean_expression > expression_cutoff,])
expressed_genes = mean_expression$gene[mean_expression$mean_expression > expression_cutoff]
cortex_only = cortex_only[cortex_only$gene %in% expressed_genes,]
total_num_tests = length(unique(cortex_only$gene))

het_mapt = mapt[mapt$haplotype=="H2/H1",]
homo1_mapt = mapt[mapt$haplotype=="H1/H1",]
homo2_mapt = mapt[mapt$haplotype=="H2/H2",]
w_mapt_half = wilcox.test(het_mapt$log2_expression, homo1_mapt$log2_expression)$p.value
w_mapt_half_corrected = wilcox.test(het_mapt$log2_expression, homo1_mapt$log2_expression)$p.value * total_num_tests
w_mapt_full = wilcox.test(homo2_mapt$log2_expression, homo1_mapt$log2_expression)$p.value
w_mapt_full_corrected = wilcox.test(homo2_mapt$log2_expression, homo1_mapt$log2_expression)$p.value * total_num_tests

het_kansl1_as = kansl1_as[kansl1_as$haplotype=="H2/H1",]
homo1_kansl1_as = kansl1_as[kansl1_as$haplotype=="H1/H1",]
homo2_kansl1_as = kansl1_as[kansl1_as$haplotype=="H2/H2",]
w_kansl1_as_half = wilcox.test(het_kansl1_as$log2_expression, homo1_kansl1_as$log2_expression)$p.value
w_kansl1_as_half_corrected = wilcox.test(het_kansl1_as$log2_expression, homo1_kansl1_as$log2_expression)$p.value * total_num_tests
w_kansl1_as_full = wilcox.test(homo2_kansl1_as$log2_expression, homo1_kansl1_as$log2_expression)$p.value
w_kansl1_as_full_corrected = wilcox.test(homo2_kansl1_as$log2_expression, homo1_kansl1_as$log2_expression)$p.value * total_num_tests

het_mapk8 = mapk8[mapk8$haplotype=="H2/H1",]
homo1_mapk8 = mapk8[mapk8$haplotype=="H1/H1",]
homo2_mapk8 = mapk8[mapk8$haplotype=="H2/H2",]
w_mapk8_half = wilcox.test(het_mapk8$log2_expression, homo1_mapk8$log2_expression)$p.value
w_mapk8_half_corrected = wilcox.test(het_mapk8$log2_expression, homo1_mapk8$log2_expression)$p.value * total_num_tests
w_mapk8_full = wilcox.test(homo2_mapk8$log2_expression, homo1_mapk8$log2_expression)$p.value
w_mapk8_full_corrected = wilcox.test(homo2_mapk8$log2_expression, homo1_mapk8$log2_expression)$p.value * total_num_tests


