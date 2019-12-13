require(tidyverse)

data = read_delim("output/candidate_snp_lists/ld_buddies_table_stage1.tsv", delim="\t")
gwas_keys = c("Alzheimers_Kunkle_2019", "Nalls_23andMe", "23andme_PD_hg38", "Chang_23andMe_Parkinsons", "Alzheimers_Jansen_2018", "Alzheimers_Lambert_2013", "Parkinsons_Pankratz_2012")
gwas_values = c("data/gwas/prepared/Alzheimers_Kunkle_2019/Alzheimers_Kunkle_2019.txt.gz",
		"data/gwas/prepared/23andme_PD_hg38.txt.gz",
		"data/gwas/prepared/23andme_PD_hg38.txt.gz",
		"data/gwas/prepared/23andme_PD_hg38.txt.gz",
		"data/gwas/prepared/Alzheimers_Jansen_2018/Alzheimers_Jansen_2018.txt.gz",
		"data/gwas/prepared/Alzheimers_Lambert_2013/Alzheimers_Lambert_2013.txt.gz",
		"data/gwas/prepared/Parkinsons_Pankratz_2012/Parkinsons_Pankratz_2012.txt.gz")

data$file = gwas_values[match(data$source_gwas, gwas_keys)]

## For each SNP, based on the GWAS field, get the associated GWAS pvalue (if it exists)
effect = rep("NA", dim(data)[1])
noneffect = rep("NA", dim(data)[1])
direction = rep("NA", dim(data)[1])
pvalue = rep(-1, dim(data)[1])
for (i in 1:dim(data)[1])
{
	print(i)
	header = strsplit(system(paste0("zcat ", data$file[i], " | head -n 1"), intern=TRUE), "\t")[[1]]
	effect_index = match("effect_allele", header)
	noneffect_index = match("non_effect_allele", header)
	if (is.na(effect_index))
	{
		effect_index = match("alt", header)
		noneffect_index = match("ref", header)
	}
	beta_index = match("beta", header)
	pval_index = match("pvalue", header)
	# Just get the first one, hopefully that'll be good enough, if not we can investigate later
	out = system(paste0("tabix ", data$file[i], " ", as.character(data$chr[i]), ":", as.character(data$pos[i]), "-", as.character(data$pos[i])), intern=TRUE)
	if (length(out) == 0)
	{
		out = system(paste0("tabix ", data$file[i], " chr", as.character(data$chr[i]), ":", as.character(data$pos[i]), "-", as.character(data$pos[i])), intern=TRUE)
		if (length(out) == 0)
		{
			print(paste0("tabix ", data$file[i], " ", as.character(data$chr[i]), ":", as.character(data$pos[i]), "-", as.character(data$pos[i])))
			print(data$source_gwas[i])
			next
		}
	}
	out = strsplit(strsplit(out, "\n")[[1]], "\t")[[1]]
	noneffect[i] = out[noneffect_index]
	effect[i] = out[effect_index]
	if (as.numeric(out[beta_index]) > 0)
	{
		direction[i] = "+"
	}
	else if (as.numeric(out[beta_index]) < 0)
	{
		direction[i] = "-"
	}
	pvalue[i] = out[pval_index]
}

data$effect_allele = effect
data$noneffect_allele = noneffect
data$direction = direction
data$pvalue = pvalue

data[data$pvalue == -1,]$pvalue = NA
coloc_mat = data %>% group_by(locus_num, source_gwas) %>% summarize(has_coloc = sum(grepl("coloc", source)) > 0)
data = left_join(data, coloc_mat, by=c("locus_num", "source_gwas"))

write.table(data, "output/candidate_snp_lists/ld_buddies_table_stage2.tsv", sep="\t", quote=FALSE, row.names=FALSE)


## For each lead colocalized SNP, get the list of all tissue/gene/GWAS combinations in which it colocalized

