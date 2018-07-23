# Author: Mike Gloudemans
# Date: 7/23/2018

import glob
import subprocess
import sys

for i in range(1250):
    print i
    eqtl_file = "/users/mgloud/projects/coloc_comparisons/output/simulations/eqtl/eqtl_sumstats{0}.txt".format(i)
    gwas_file = "/users/mgloud/projects/coloc_comparisons/output/simulations/gwas/gwas_sumstats{0}.txt".format(i)

    #
    # eQTL
    #

    # Write liftable file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/liftable.bed", "w") as w:
        with open(eqtl_file) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                w.write("{0}\t{1}\t{2}\t{3}\n".format("chr"+data[3], data[4], int(data[4])+1, data[1]))
    # Lift it
    subprocess.check_call("liftOver /users/mgloud/projects/coloc_comparisons/tmp/liftable.bed /users/mgloud/software/liftOver/chains/hg19ToHg38.over.chain.gz /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed /users/mgloud/projects/coloc_comparisons/tmp/unlifted.bed", shell=True)
    
    # Merge back together with the original version to save metadata; scrap unnecessary rows
    subprocess.check_call("sort -k4,4 /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed > /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed", shell=True)
    subprocess.check_call("sort -k2,2 {0} > /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed".format(eqtl_file), shell=True)
    subprocess.check_call("join -1 4 -2 2 /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed | sed s/\ /\\\t/g > /users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed", shell=True)

    # Reformat it one more time

    with open("/users/mgloud/projects/coloc_comparisons/output/simulations/hg38/eqtl/eqtl_sumstats{0}.txt".format(i), "w") as w:
        w.write("rsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref_allele\talt_allele\teffect_af\teffect_size\tse\tzscore\tpvalue\tN\n")
        with open("/users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed") as f:
            for line in f:
                data = line.strip().split()
                # temporary, since first pass screwed up the sample size by not recording it
                data = data + [-1]
                w.write("{0}\t{1}_{2}_{3}_{4}_hg38\thg38\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(data[4], data[1], data[2], data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15]))
                
    #
    # GWAS
    #

    # Write liftable file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/liftable.bed", "w") as w:
        with open(gwas_file) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                w.write("{0}\t{1}\t{2}\t{3}\n".format("chr"+data[3], data[4], int(data[4])+1, data[1]))
    # Lift it
    subprocess.check_call("liftOver /users/mgloud/projects/coloc_comparisons/tmp/liftable.bed /users/mgloud/software/liftOver/chains/hg19ToHg38.over.chain.gz /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed /users/mgloud/projects/coloc_comparisons/tmp/unlifted.bed", shell=True)
    
    # Merge back together with the original version to save metadata; scrap unnecessary rows
    subprocess.check_call("sort -k4,4 /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed > /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed", shell=True)
    subprocess.check_call("sort -k2,2 {0} > /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed".format(gwas_file), shell=True)
    subprocess.check_call("join -1 4 -2 2 /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed | sed s/\ /\\\t/g > /users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed", shell=True)

    # Reformat it one more time

    with open("/users/mgloud/projects/coloc_comparisons/output/simulations/hg38/gwas/gwas_sumstats{0}.txt".format(i), "w") as w:
        w.write("rsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref_allele\talt_allele\teffect_af\teffect_size\tse\tzscore\tpvalue\tn_cases\tn_controls\n")
        with open("/users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed") as f:
            for line in f:
                data = line.strip().split()
                w.write("{0}\t{1}_{2}_{3}_{4}_hg38\thg38\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(data[4], data[1], data[2], data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[16]))
