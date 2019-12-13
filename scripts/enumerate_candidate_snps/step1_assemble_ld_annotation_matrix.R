require(tidyverse)

#
# Constants
#

coloc_dir = "output/colocalization/2019-08-24_19-43-32.358274_ad-pd/"
coloc_file_pattern = "_finemap_clpp_status.txt"

gwas_dir = "output/test-snps"
gwas_file_pattern = "gz.tests.txt"

get_ld_buddies = function()
{
	#
	# Load LD buddies for all GWAS of interest
	#
	dir_files = dir("output/ld_buddies")
	out_files = dir_files[str_detect(dir_files, "with.ld.buddies.txt")]
	ld_buddy_files = paste0("output/ld_buddies/", out_files)
	info = lapply(ld_buddy_files, read_delim, delim="\t")

	#
	# Concatenate all of them into a single table
	#

	# For each SNP, indicate which GWAS it was from, and LD with the lead SNP
	for (i in 1:length(ld_buddy_files))
	{
		info[[i]]$source_gwas = out_files[i]
	}

	ld_buddies = do.call(rbind, info)

	# Get rid of unconventional chromosomes
	ld_buddies = ld_buddies[!grepl("_", ld_buddies$chr),]

	ld_buddies$source_gwas = gsub("_txt_gz_finemap_clpp_status", "", ld_buddies$source_gwas)
	ld_buddies$source_gwas = gsub(".with.ld.buddies.txt", "", ld_buddies$source_gwas)

	ld_buddies = ld_buddies[!duplicated(ld_buddies[c("chr", "pos", "source_gwas", "source")]),]

	# If they appear as both GWAS hits and coloc hits, then combine these into a single entry
	ld_buddies = ld_buddies %>% group_by(chr, pos, r2_with_ld_tag, ld_tag_chr, ld_tag_pos, source_gwas) %>% summarize(source = paste(sort(source), collapse=";"))

	# SNPs shouldn't be listed as LD buddies for another SNP, if they're already significant
	# SNPs of their own accord
	# This never seems to happen though, which is good I think
	lead_snps = ld_buddies %>% filter(source != "LD")
	lead_snp_ids = paste(lead_snps$chr, lead_snps$pos, lead_snps$source_gwas, sep="_")
	ld_buddies = ld_buddies[!((paste(ld_buddies$chr, ld_buddies$pos, ld_buddies$source_gwas, sep="_") %in% lead_snp_ids) & (ld_buddies$source == "LD")),]

	return(ld_buddies)
}

ld_buddies = get_ld_buddies()

#
# Add other metadata
#

# If we really want to add GWAS p-values, then just pull them with tabix...

# And get tissues in which it colocalized, I guess...

# Count the number of unique loci that we're testing
# Function inputs a vector of SNPs and clusters them into loci
# by distance
group_to_loci = function(x)
{
	ids = unique(as.character(x))
	ids = ids[order(ids)]
	chr = sapply(as.character(ids), function(x) {strsplit(x, "_")[[1]][1]})
	pos = as.numeric(sapply(as.character(ids), function(x) {strsplit(x, "_")[[1]][2]}))
	loc_nums = rep(0, length(ids))
	loc_nums[1] = 1
	for (i in 2:length(ids))
	{
		# Check if there's a SNP above in the list within 1 MB of this SNP]
		same = (chr[1:(i-1)] == chr[i]) & (abs(pos[1:(i-1)] - pos[i]) < 500000)
		if (length(which(same)) != 0)
		{
			loc_nums[i] = loc_nums[which(same)[1]]
		}
		else
		{
			loc_nums[i] = max(loc_nums) + 1
		}
	}

	mapped_loci = sapply(x, function(j)
	{
	     loc_nums[which(ids == j)]
	})

	return(mapped_loci)
}
ld_buddies$snp_id = paste(ld_buddies$chr, ld_buddies$pos, sep="_")
ld_buddies$locus_num = group_to_loci(ld_buddies$snp_id)

ld_buddy_counts = ld_buddies %>% group_by(locus_num, source_gwas) %>% summarize(number_ld_buddies = length(source))
ld_buddies = left_join(ld_buddies, ld_buddy_counts, by=c("locus_num", "source_gwas"))

ld_buddies$ld_tag_locus = paste(ld_buddies$ld_tag_chr, ld_buddies$ld_tag_pos, sep="_")

# Get rsid by pulling from 1000 Genomes with tabix...
get_rsid_by_hg19 = function (x)
{
       chrom = strsplit(x, "_")[[1]][1]
       pos = strsplit(x, "_")[[1]][2]

       if (chrom == "X")
       {
	       rsid = system(sprintf("tabix data/indexed-dbsnp/hg19/common_all_20170710.vcf.gz X:%s-%s | cut -f3", pos, pos), intern=TRUE)
       }
       else if (chrom == "Y")
       {
	       rsid = system(sprintf("tabix data/indexed-dbsnp/hg19/common_all_20170710.vcf.gz Y:%s-%s | cut -f3", pos, pos), intern=TRUE)
       }
       else
       {
	       rsid = system(sprintf("tabix data/indexed-dbsnp/hg19/common_all_20170710.vcf.gz %s:%s-%s | cut -f3", chrom, pos, pos), intern=TRUE)
       }

       rsid = rsid[grepl("rs", rsid)]
       if (length(rsid) == 0)
       {
	       return(x)
       }
       else
       {
	       return(rsid[1])
       }

}

# Get rsid by pulling from 1000 Genomes with tabix...
get_rsid_by_hg38 = function (x)
{
       chrom = strsplit(x, "_")[[1]][1]
       pos = strsplit(x, "_")[[1]][2]

       if (chrom == "X")
       {
	       rsid = system(sprintf("tabix data/indexed-dbsnp/hg38/common_all_20170710.vcf.gz X:%s-%s | cut -f3", pos, pos), intern=TRUE)
       }
       else if (chrom == "Y")
       {
	       rsid = system(sprintf("tabix data/indexed-dbsnp/hg38/common_all_20170710.vcf.gz Y:%s-%s | cut -f3", pos, pos), intern=TRUE)
       }
       else
       {
	       rsid = system(sprintf("tabix data/indexed-dbsnp/hg38/common_all_20170710.vcf.gz %s:%s-%s | cut -f3", chrom, pos, pos), intern=TRUE)
       }

       rsid = rsid[grepl("rs", rsid)]
       if (length(rsid) == 0)
       {
	       return(x)
       }
       else
       {
	       return(rsid[1])
       }

}

ld_buddies$snp_id = paste(unlist(ld_buddies$chr), unlist(ld_buddies$pos), sep="_")
ld_buddies$rsid = sapply(ld_buddies$snp_id, get_rsid_by_hg38)

# See if R is smart enough to do this without totally breaking...

dbsnp = read_delim("data/dbsnp/sorted_hg19_snp150.txt.gz", col_names=FALSE, delim="\t")
colnames(dbsnp) = c("chrom_hg19", "snp_pos_hg19", "rsid")

# Throw away entries that aren't on a valid chromosome
dbsnp = dbsnp[!grepl("_", dbsnp$chrom_hg19),]

ld_buddies = left_join(ld_buddies, dbsnp, by="rsid")
ld_buddies$chrom_hg19 = as.numeric(sub("chr", "", ld_buddies$chrom_hg19))

#
# Save matrix stage 1
# Send to Ryan and other collaborators
#
write_delim(ld_buddies, "output/candidate_snp_lists/ld_buddies_table_stage1.tsv", delim="\t")

