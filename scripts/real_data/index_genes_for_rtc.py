import subprocess
import gzip

# Load persistent data
pheno_file = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz"
rtc_index_dir = "/users/mgloud/projects/coloc_comparisons/data/rtc-gene-index"

with gzip.open(pheno_file) as f:
    f.readline()
    f.readline()
    with open("{0}/header.txt".format(rtc_index_dir), "w") as w:
        w.write(f.readline())

    for line in f:
        gene = line.strip().split()[0].replace("ENSG00000", "")
        g1 = gene[:2]
        g2 = gene[2:4]
        subprocess.check_call("mkdir -p {0}/{1}/{2}".format(rtc_index_dir, g1, g2), shell=True)

        with open("{0}/{1}/{2}/{3}.txt".format(rtc_index_dir, g1, g2, gene), "w") as w:
            w.write(line)

