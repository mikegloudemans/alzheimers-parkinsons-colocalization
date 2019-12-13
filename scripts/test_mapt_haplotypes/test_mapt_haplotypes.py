import pandas as pd
import subprocess
from StringIO import StringIO
import gzip
import sys

# Load GTEx metadata (probably with pandas)
metadata = pd.read_csv("data/gtex/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Get IDs of the individuals with any brain tissue
#brain = metadata[metadata['SMTSD'].apply(lambda x: "Brain" in str(x))]
brain = metadata
brain_samps = {}
for i in range(len(list(brain["SAMPID"]))):
    brain_samps[list(brain["SAMPID"])[i]] = list(brain["SMTSD"])[i]
brain_ids = list(set(["-".join(bi) for bi in [m.split("-")[:2] for m in list(brain["SAMPID"])]]))

# Subset VCF to region of interest, load it, subset to proper individuals
with gzip.open("data/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz") as f:
    for line in f:
        if line.startswith("#CHROM"):
            header = line.strip().split()
            break
stream = StringIO(subprocess.check_output("tabix data/gtex/GTEx_Analysis_2017-06-05_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz chr17:44577102-47577102", shell=True))
genos = pd.read_csv(stream, sep="\t", names=header)

keep_head = list(genos.columns.values[:9]) + brain_ids
keep_head = [k for k in keep_head if k in list(genos.columns.values)]
brain_genos = genos[keep_head]

# Which allele corresponds to haplotype h1?
h1_codes = {}
with open("data/mapt/MAPT_dbsnp151_HapSpecificSNPs_fixed.txt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        data = line.strip().split()
        ID = tuple(data[:2])
        if data[10].split(":")[0] == "0|0":
            h1_codes[ID] = (data[3], data[4])   # ref is h1
        elif data[10].split(":")[0] == "1|1":
            h1_codes[ID] = (data[4], data[3])   # alt is h1
        else:
            print "Something weird happened", data[10]

# Map individuals to a specific haplotype (should hopefully be obvious based on just a few markers)
ind_haplos = {}
# For each individual
for ind in brain_genos.columns.values[9:]:
    ind_haplos[ind] = {"h2": 0, "h1": 0}
    subset = brain_genos[['#CHROM', "POS", "REF", "ALT", ind]]
    for i in range(subset.shape[0]):
        combo = (subset.iloc[i, 0], str(subset.iloc[i, 1]))
        haplo = subset.iloc[i, 4].split(":")[0].split("/")
        if combo in h1_codes:
            if subset.iloc[i, 2] == h1_codes[combo][0]:
                # Ref in GTEx is H1
                ind_haplos[ind]['h1'] += haplo.count("0")
            if subset.iloc[i, 3] == h1_codes[combo][1]:
                # Alt in GTEx is H2
                ind_haplos[ind]['h2'] += haplo.count("1")
            if subset.iloc[i, 2] == h1_codes[combo][1]:
                # Ref in GTEx is H2
                ind_haplos[ind]['h2'] += haplo.count("0")
            if subset.iloc[i, 3] == h1_codes[combo][0]:
                # Alt in GTEx is H1
                ind_haplos[ind]['h1'] += haplo.count("1")
    print ind_haplos[ind]
    h2_ratio = ind_haplos[ind]['h2'] * 1.0 / (ind_haplos[ind]['h2'] + ind_haplos[ind]['h1'])
    if ind_haplos[ind]['h2'] + ind_haplos[ind]['h1'] < 4000:
        print "Strange case: not enough determinable genotypes"
        continue
        ind_haplos[ind] = "weirdo"
    if h2_ratio < 0.01:
        print "H1/H1"
        ind_haplos[ind] = "H1/H1"
    elif h2_ratio > 0.99:
        print "H2/H2"
        ind_haplos[ind] = "H2/H2"
    elif h2_ratio > 0.45 and h2_ratio < 0.55:
        print "H2/H1"
        ind_haplos[ind] = "H2/H1"
    else:
        print "Strange case"

print ind_haplos
print "H1/H1:", ind_haplos.values().count("H1/H1")
print "H2/H1:", ind_haplos.values().count("H2/H1")
print "H2/H2:", ind_haplos.values().count("H2/H2")

# TODO: Cross-reference these genotypes with our actual VCFs to make sure we agree with our conclusion here
# just a few spot checks

saveable_genes = {} 
# Load location of each gene and
# get the set of genes that are in our region of interest
with gzip.open("data/gencode/gencode.v31.annotation.gtf.gz") as f:
    for i in range(5):
        f.readline()
    for line in f:
        data = line.strip().split()
        if data[2] != "gene":
            continue
        if data[0] != "chr17":
            continue
        if int(data[4]) < 44577102:
            continue
        if int(data[3]) > 47577102:
            continue

        # Swap "start" and "end" if gene is on reverse strand
        if data[6] == "+":
            saveable_genes[data[9].replace(";", "").replace('\"', "").split(".")[0]] = (data[3], data[4])
        elif data[6] == "-":
            saveable_genes[data[9].replace(";", "").replace('\"', "").split(".")[0]] = (data[4], data[3])

# Load expression matrix and get the genes we care about
# Output to "melted" data frame for R analysis: individual, haplotype, gene, expression
with gzip.open("data/gtex/GTEx_Analysis_2017-06-05_v8/rna_seq/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz") as f:
    with open("output/mapt/mapt_expression_matrix.txt", "w") as w: 
        w.write("gene\tsample\thaplotype\tindividual\texpression\tstart\tend\tbrain_region\n")
        f.readline()
        f.readline()
        header = f.readline().strip().split()
        for line in f:
            data = line.strip().split()
            gene = data[0].split(".")[0]
            # Do we care about this particular gene?
            if gene not in saveable_genes:
                continue
            print gene
            tissue_matches = [i for i in range(len(header)) if header[i] in brain_samps]

            # Save key information for this individual/gene combo
            for tm in tissue_matches:
                sample = header[tm]
                ind = "-".join(sample.split("-")[:2])
                if ind not in ind_haplos:
                    continue
                haplo = ind_haplos[ind]
                expression = data[tm]
                w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(gene, sample, haplo, ind, expression, saveable_genes[gene][0], saveable_genes[gene][1], brain_samps[sample]))




