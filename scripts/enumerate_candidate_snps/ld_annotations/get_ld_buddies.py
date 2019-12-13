import glob
import subprocess
import os
import sys

gwas_threshold = 5e-8

def main():
    # Get SNPs to test for each trait
    # For some of them we have full GWAS sumstats
    snps_to_test = get_snps_to_test()

    # Now, get LD for all variants of interest
    for trait in snps_to_test:
        print trait
        # Remove the previous LD collection file for this trait, if one exists
        with open("output/ld_buddies/{0}.with.ld.buddies.txt".format(trait), "w") as w:
            w.write("chr\tpos\tr2_with_ld_tag\tld_tag_chr\tld_tag_pos\ttrait\tsource\n")
            # We write a p-value only for those SNPs that we know it, which is OK because
            # all other SNPs at this point should have less significant p-values.
            for gs in snps_to_test[trait]:
                w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(gs[0], gs[1], 1.0, gs[0], gs[1], trait, gs[2]))
                snps_to_test[trait] = [(s[0], s[1]) for s in snps_to_test[trait]]

            # Store a list of all SNPs that have LD R^2 > 0.8 with at least one of the GWAS SNPs
            closest_ld_snps = {} 

            # Loop through GWAS significant SNPs
            for gs in snps_to_test[trait]:
                # Tabix to get all SNPs in LD with this SNP...although...is it bi-directional?
                # Answer: yes (thank you Peyton Greenside for assembling this file and saving me
                #   hours, days, possibly years of computational time with VCFtools)
                # If not, will need to get a slightly wider range
                output = subprocess.check_output("tabix data/ld/EUR_geno_hg38.txt.gz chr{0}:{1}-{1}".format(gs[0],gs[1],gs[1]), shell=True)
                if output.strip() == "":
                    # No LD pairs
                    continue
                for line in output.strip().split("\n"):
                    data = line.strip().split("\t")
                    
                    # Keep track of the closest LD SNP of the reference set (and which locus number that falls within)
                    # Will fix locus numbers later by merging nearby ones, I guess
                   
                    data[9] = data[9].replace("chr", "")
                    data[11] = data[11].replace("chr", "")

                    # Make sure the SNP's not in our list already due to passing GWAS cutoff
                    if (data[11], data[12]) in snps_to_test[trait]:
                        continue

                    # Add it to the list, if LD high enough
                    assert(float(data[8])) > 0.8
                    if (data[11], data[12]) in closest_ld_snps:
                        if float(data[8]) > closest_ld_snps[(data[11], data[12])][0]:
                            closest_ld_snps[(data[11], data[12])] = (float(data[8]), data[9], data[10]) # (Best LD, best LD snp)
                    else:
                        closest_ld_snps[(data[11], data[12])] = (float(data[8]), data[9], data[10]) # (Best LD, best LD snp)

            # Go back through the list of LD SNPs and output them
            for cls in closest_ld_snps:
                main_snp = closest_ld_snps[cls]
                w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(cls[0].replace("chr", ""), cls[1], main_snp[0], main_snp[1].replace("chr", ""), main_snp[2], trait, "LD"))
                


def get_snps_to_test():
    traits = glob.glob("output/test-snps/*.gz.tests.txt")
    snps_to_test = {}

    for trait in traits:
        short_trait = trait.replace("output/test-snps/", "").replace(".txt.gz.tests.txt", "")
        snps_to_test[short_trait] = set([])
        # Compute LD buddies for each of these SNPs using VCF tools in 1K genomes, up to R^2 > 0.8
        with open(trait) as f:
            for line in f:
                chrom = line.strip().split()[0]
                pos = line.strip().split()[1]
                pvalue = line.strip().split()[2]
                if float(pvalue) > 5e-8:
                    continue
                snps_to_test[short_trait].add((chrom, pos, "GWAS"))

    # For some traits, all we have is a list from Ryan
    other_traits = glob.glob("data/gwas_top_hits/*")
    for trait in other_traits:
        short_trait = trait.replace("_LeadSNPs.txt", "").replace("data/gwas_top_hits/", "")
        snps_to_test[short_trait] = set([])
        with open(trait) as f:
            f.readline()
            for line in f:
                chrom = line.strip().split()[1]
                pos = line.strip().split()[3]
                snps_to_test[short_trait].add((chrom, pos, "Nalls-Chang"))

    # And then there are the coloc lists...
    colocs = glob.glob("output/colocalization/2019-08-24_19-43-32.358274_ad-pd/*_finemap_clpp_status.txt")
    for coloc in colocs:
        short_coloc = coloc.replace("_txt_gz_finemap_clpp_status.txt", "").replace("output/colocalization/2019-08-24_19-43-32.358274_ad-pd/", "")
        with open(coloc) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                chrom = data[0].split("_")[0]
                pos = data[0].split("_")[1]
                clpp = float(data[5])
                if clpp < 0.01:
                    continue
                log_gwas_p = float(data[6])
                if log_gwas_p < 4:
                    continue
                n_snps = int(data[4])
                if n_snps < 20:
                    continue

                snps_to_test[short_coloc].add((chrom, pos, "coloc"))

    for trait in snps_to_test:
        snps_to_test[trait] = list(snps_to_test[trait])
    return snps_to_test

if __name__ == "__main__":
    main()
