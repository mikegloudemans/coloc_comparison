#!/usr/bin/python
# Author: Mike Gloudemans
# 7/18/2018

# Write a genetic map file (showing recombination rates)
# including every 1000 Genomes varaint.
# Then bgzip and tabix it.

import gzip
import subprocess

# Genetic map
#for chrom in (range(1,23) + ["X"]):
for chrom in (range(7,23) + ["X"]):
    kg_file = "/mnt/lab_data/montgomery/shared/1KG/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz".format(chrom)

    filename = "/mnt/lab_data/montgomery/nicolerg/genetic-map/genetic_map_GRCh37_chr{0}.txt".format(chrom)
    output = "/users/mgloud/projects/coloc_comparisons/data/genetic-map/extended_" + filename.split("/")[-1]

    last_pos = 0
    next_pos = 0
    last_map = 0
    next_map = 0

    with open(output, "w") as w:
        w.write("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n")
        with gzip.open(kg_file) as f1k:
            with open(filename) as f:
                f.readline()
                for kg_var in f1k:
                    if kg_var.startswith("#"):
                        continue

                    kg_var = int(kg_var.split()[1])
                    if kg_var > next_pos:
                        map_var = f.readline()
                        if map_var == "":
                            break
                        map_var = map_var.strip().split()


                        last_pos = next_pos
                        next_pos = int(map_var[1])
                        last_map = next_map
                        next_map = float(map_var[3])
                        current_rate = map_var[2]

                    current_map = last_map + ((kg_var - last_pos) * 1.0 / (next_pos - last_pos)) * (next_map - last_map)

                    w.write("{3}\t{0}\t{1}\t{2}\n".format(kg_var, current_rate, current_map, chrom))

        subprocess.check_call("cat {0} | uniq | bgzip -f > {0}.gz".format(output), shell="True")
        subprocess.check_call("tabix -f -S 1 -s 1 -b 2 -e 2 " + output + ".gz", shell=True)


