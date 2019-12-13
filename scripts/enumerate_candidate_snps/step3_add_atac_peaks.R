#
# Stage 3: Add in annotations of which ones overlap peaks
#

require(tidyverse)
filtered_ld_buddies = read_delim("output/candidate_snp_lists/ld_buddies_table_stage2.tsv", delim="\t")

atac_groups = list()

atac_dir = "data/atacseq/kundaje_version/support30percent/"
atac_files = dir(atac_dir)
atac_groups[[atac_dir]] = atac_files[!grepl(".gz", atac_files) & !(grepl("allregions", atac_files))]

atac_dir = "data/atacseq/kundaje_version/support30percent_BroadRegions/"
atac_files = dir(atac_dir)
atac_groups[[atac_dir]] = atac_files[grepl(".merged.bed.gz", atac_files) & !(grepl("tbi", atac_files)) & !(grepl("allregions", atac_files))]

atac_dir = "data/atacseq/kundaje_version/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge"
atac_files = dir(atac_dir)
atac_groups[[atac_dir]] = atac_files[grepl(".narrowPeak.gz", atac_files) & !(grepl("tbi", atac_files)) & !(grepl("allregions", atac_files))]


atac_data = list()
for (group in names(atac_groups))
{
	atac_data[[group]] = list()
	for (tissue in atac_groups[[group]])
	{
		atac_data[[group]][[tissue]] = read.table(paste0(group, "/", tissue), col.names=c("chr", "start", "end"), stringsAsFactors=FALSE)
	}
}

names(atac_data) = c("narrow_tissue_regions", "broad_tissue_regions", "single_cell")

size = dim(filtered_ld_buddies)[1]
nearest_peaks_to_snps = list()
for (group in names(atac_data))
{
	nearest_peaks_to_snps[[group]] = data.frame(list(direct_atac_overlap=rep(FALSE, size), containing_atac_tissues=rep("", size), nearest_atac_tissue=rep("", size), start=rep(0, size), end=rep(0, size), dist=rep(Inf, size)),stringsAsFactors=FALSE)


	# For each SNP we want to track:
	# - Does it overlap an ATAC-seq peak
	# - How close is the closest peak
	# - How many cell types have overlapping peaks
	# - Which cell type(s) has overlapping peaks
	# - Which cell type peak is closest
	# 
	# Then repeat for all peaks in single cell clusters

	for (i in 1:size)
	{
		if (i %% 100 == 0)
		{
			print(i)
		}
		if (is.na(filtered_ld_buddies$chr[i]) | is.na(filtered_ld_buddies$pos[i]))
		{
			nearest_peaks_to_snps[[group]]$containing_atac_tissues[i] = NA
			nearest_peaks_to_snps[[group]]$nearest_atac_tissue[i] = NA
			nearest_peaks_to_snps[[group]]$start[i] = NA
			nearest_peaks_to_snps[[group]]$end[i] = NA
			nearest_peaks_to_snps[[group]]$dist[i] = NA
			nearest_peaks_to_snps[[group]]$direct_atac_overlap[i] = NA
			next
		}
		chrom = paste0("chr", filtered_ld_buddies$chr[i])
		pos = as.numeric(filtered_ld_buddies$pos[i])


		for (ad in names(atac_data[[group]]))
		{
			ad_mat = atac_data[[group]][[ad]]
			ad_short = strsplit(ad, "\\.")[[1]][1]
			sub = ad_mat[ad_mat$chr == chrom,]
			# We'll measure distance to the midpoint of the peak interval
			dist_to_peak = min(abs(((sub$start + sub$end) / 2) - pos))
			peak_loc = sub[which(abs(((sub$start + sub$end) / 2) - pos) == dist_to_peak),]
			
			# Test whether it's actually IN the peak
			containing_peak = sub[sub$start <= pos & sub$end >= pos,]
			if (dim(containing_peak)[1] != 0)
			{
				peak_loc = containing_peak
				nearest_peaks_to_snps[[group]]$containing_atac_tissues[i] = paste(nearest_peaks_to_snps[[group]]$containing_atac_tissues[i], ad_short, sep=",") 
				nearest_peaks_to_snps[[group]]$nearest_atac_tissue[i] = nearest_peaks_to_snps[[group]]$containing_atac_tissues[i]
				nearest_peaks_to_snps[[group]]$start[i] = peak_loc$start[1]
				nearest_peaks_to_snps[[group]]$end[i] = peak_loc$end[1]
				nearest_peaks_to_snps[[group]]$dist[i] = dist_to_peak
				nearest_peaks_to_snps[[group]]$direct_atac_overlap[i] = TRUE
			}
			
			else if (dist_to_peak < abs(nearest_peaks_to_snps[[group]]$dist[i]))
			{
				nearest_peaks_to_snps[[group]]$nearest_atac_tissue[i] = ad_short
				nearest_peaks_to_snps[[group]]$start[i] = peak_loc$start[1]
				nearest_peaks_to_snps[[group]]$end[i] = peak_loc$end[1]
				nearest_peaks_to_snps[[group]]$dist[i] = min(c(abs(sub$start - pos), abs(sub$end - pos)))
			}
		}
	}
	nearest_peaks_to_snps[[group]]$containing_atac_tissues = sapply(nearest_peaks_to_snps[[group]]$containing_atac_tissues, function(x)
									{
										if (is.na(x) | (nchar(x) == 0))
										{
											return(x)
										}
										if (substring(x,1,1) == ",")
										{
											return(substring(x, 2))
										}
										return(x)
									})
	nearest_peaks_to_snps[[group]]$nearest_atac_tissue = sapply(nearest_peaks_to_snps[[group]]$nearest_atac_tissue, function(x)
									{
										if (is.na(x) | (nchar(x) == 0))
										{
											return(x)
										}
										if (substring(x,1,1) == ",")
										{
											return(substring(x, 2))
										}
										return(x)
									})
	nearest_peaks_to_snps[[group]]$nearest_atac_tissue[nearest_peaks_to_snps[[group]]$nearest_atac_tissue == ""] = "none"
	nearest_peaks_to_snps[[group]]$containing_atac_tissues[nearest_peaks_to_snps[[group]]$containing_atac_tissues == ""] = "none"
}


# Add results to the matrix and then output it (part 2)
for (group in names(nearest_peaks_to_snps))
{
	group_prefix = strsplit(group, "/")[[1]]
	group_prefix = group_prefix[length(group_prefix)]
        names(nearest_peaks_to_snps[[group]]) = paste(names(nearest_peaks_to_snps[[group]]), group_prefix, sep="_")
	filtered_ld_buddies = cbind(filtered_ld_buddies, nearest_peaks_to_snps[[group]])
}

write_delim(filtered_ld_buddies, "output/candidate_snp_lists/ld_buddies_table_stage3.tsv", delim="\t")



# Then do any other meta-analysis

# For example...
# Do most loci have at least one SNP overlapping an LD SNP? (try plotting some of them)
# How close are the SNPs on average?
# Any specifically interesting loci where a tie might be "broken" by this LD pattern?

