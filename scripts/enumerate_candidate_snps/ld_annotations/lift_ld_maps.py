import gzip
import subprocess

# With the help of dbSNP, convert LD pairwise map from hg19 coordinates
# to hg38 coordinates

def main():
    db38 = load_hg38_rsid_keys()

    with gzip.open("data/ld/EUR_geno.txt.gz") as f:
        with open("data/ld/EUR_geno_hg38.presort.txt", "w") as w:
            for line in f:
                data = line.strip().split()
                if data[3] in db38:
                    snp1 = db38[data[3]]
                else:
                    continue
                if data[7] in db38:
                    snp2 = db38[data[7]]
                else:
                    continue

                w.write(line.strip() + "\t" + snp1[0] + "\t" + snp1[1] + "\t" + snp2[0] + "\t" + snp2[1] + "\n")
    subprocess.check_call('''echo "snp1_chrom_hg19\tsnp1_pos_hg19\tsnp1_pos2_hg19\tsnp1_rsid_hg19\tsnp2_chrom_hg19\tsnp2_pos_hg19\tsnp2_pos2_hg19\tsnp2_rsid_hg19\tld_r2\tsnp1_chrom_hg38\tsnp1_pos_hg38\tsnp2_chrom_hg38\tsnp2_pos_hg38" > data/ld/EUR_geno_hg38.txt''', shell=True)
    subprocess.check_call("sort -k10,10 -k11,11n data/ld/EUR_geno_hg38.presort.txt >> data/ld/EUR_geno_hg38.txt", shell=True)
    subprocess.check_call("bgzip -f data/ld/EUR_geno_hg38.txt", shell=True)
    subprocess.check_call("tabix -f -S 1 -s 10 -b 11 -e 11 data/ld/EUR_geno_hg38.txt.gz", shell=True)


def load_hg38_rsid_keys():
    return load_rsid_keys(rsid_to_pos_file="data/dbsnp/sorted_1kg_matched_hg38_snp150.txt.gz", \
            pos_to_rsid_file="data/dbsnp/sorted_1kg_matched_hg19_snp150.txt.gz")

def load_rsid_keys(rsid_to_pos_file, pos_to_rsid_file):

    rsid_to_pos = {}

    with gzip.open(rsid_to_pos_file) as f:
        line_no = 0
        for line in f:
            data = line.strip().split()
            chrom = data[0]
            rs_no = data[2]

            rsid_to_pos[rs_no] = (chrom, data[1])

            line_no += 1
            if line_no % 10000000 - 1 == -1:
                print line_no

    return rsid_to_pos

if __name__ == "__main__":
    main()
