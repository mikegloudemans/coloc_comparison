# Author: Mike Gloudemans
# Date: 1/25/2019
#
# Compare the results of COLOC when using data with the same p-values,
# but varying settings of the sample size parameter.
#

import glob

with open("coloc_sample_size_comparison_results.txt", "w") as w:

    # Write header
    with open(glob.glob("/users/mgloud/projects/brain_gwas/output/coloc-sample-size-test/eqtl_10/gwas_10/*/*h4*")[0]) as f:
        w.write(f.readline().strip())
        w.write("\teqtl_N\tgwas_N\n")

    size = [10, 100, 1000, 10000, 100000]

    for eqtl_ss in size:
        for gwas_ss in size:
            files = glob.glob("/users/mgloud/projects/brain_gwas/output/coloc-sample-size-test/eqtl_{0}/gwas_{1}/*/*h4*".format(eqtl_ss, gwas_ss))
            for filename in files:
                with open(filename) as f:
                    f.readline()
                    data = f.readline()
                    if data == "":
                        continue
                    w.write(data.strip())
                    w.write("\t{0}\t{1}\n".format(eqtl_ss, gwas_ss))


