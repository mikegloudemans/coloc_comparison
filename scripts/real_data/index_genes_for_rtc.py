import subprocess
import gzip
import os

# RTC requires that we input gene expression. Here we provide a convenient indexed directory structure
# so that the lookups only take a fraction of the time

def index_genes(rtc_index_dir = "/users/mgloud/projects/coloc_comparisons/data/rtc-gene-index"):

    rtc_index_dir = "/users/mgloud/projects/coloc_comparisons/data/rtc-gene-index"

    # Index genes if not already indexed.
    if "header.txt" in os.listdir(rtc_index_dir):
        return "Already indexed"


    # Load persistent data
    pheno_file = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz"

    with gzip.open(pheno_file) as f:
        f.readline()
        f.readline()
        header = f.readline()   # Load header; write it later to indicate we're finished
        for line in f:
            gene = line.strip().split()[0].replace("ENSG00000", "")
            g1 = gene[:2]
            g2 = gene[2:4]
            subprocess.check_call("mkdir -p {0}/{1}/{2}".format(rtc_index_dir, g1, g2), shell=True)

            with open("{0}/{1}/{2}/{3}.txt".format(rtc_index_dir, g1, g2, gene), "w") as w:
                w.write(line)

    # Since we check for the header to see if indexing is done,
    # don't write the header until the indexing is completed.
    # Otherwise if the program is interrupted, we'll never finish
    with open("{0}/header.txt".format(rtc_index_dir), "w") as w:
        w.write(header)


if __name__ == "__main__":
    index_genes()
