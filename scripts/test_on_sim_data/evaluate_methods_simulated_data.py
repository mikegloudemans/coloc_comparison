#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 7/23/2018

import glob
import subprocess
import sys
import time
import json
from multiprocessing import Pool
import traceback

def main():
    
    subprocess.check_call("rm -rf /users/mgloud/projects/coloc_comparisons/tmp/simulations", shell=True)
    subprocess.call("mkdir -p /users/mgloud/projects/coloc_comparisons/tmp/simulations", shell=True)

    config_file = sys.argv[1]
    with open(config_file) as f:
        sf = f.read()
    settings = json.loads(sf)

    # For each simulated colocalization:

    base_dir = settings["sim_base_dir"]
    if base_dir[-1] == "/":
        base_dir = base_dir[:-1]
    base_last_dir = base_dir.strip().split("/")[-1]

    subprocess.call("rm -rf /users/mgloud/projects/brain_gwas/output/simulation-coloc-comparisons".format(base_last_dir), shell=True)

    # Load answer key
    # Note that it doesn't really matter that we're choosing the specific lead SNP,
    # because it's only used to select the wider region we're testing.
    answers = {}
    with open(base_dir + "/answer_key.txt") as f:
        f.readline()
        for line in f:
            data = line.strip().split("\t")
            answers[int(data[0])] = {}
            answers[int(data[0])]["seed_chrom"] = int(data[1])
            answers[int(data[0])]["seed_pos"] = int(data[2])
            answers[int(data[0])]["n_cases"] = int(data[4])
            answers[int(data[0])]["n_controls"] = int(data[5])
            answers[int(data[0])]["n_eqtl"] = int(data[6])

    # Dispatch workers
    pool = Pool(int(settings["max_threads"]))
    for i in range(int(settings["num_tests"])):
        pool.apply_async(test_wrapper, args=(i, answers, config_file, base_dir, base_last_dir))
    pool.close()
    pool.join()

    # Clean up tmp files after we're done running
    subprocess.call("rm -rf /users/mgloud/projects/coloc_comparisons/tmp/simulations*", shell=True)


def test_wrapper(i, answers, config_file, base_dir, base_last_dir):
    try:
        test(i, answers, config_file, base_dir, base_last_dir)
    except Exception:
        # Exceptions should already be logged by the coloc pipeline itself
        traceback.print_exc(file=sys.stdout)

def test(i, answers, config_file, base_dir, base_last_dir):


    with open(config_file) as f:
        sf = f.read()
    settings = json.loads(sf)["coloc_settings"]

    n_gwas = answers[i]["n_cases"] + answers[i]["n_controls"]
    n_eqtl = answers[i]["n_eqtl"]
    cc_ratio_gwas = answers[i]["n_cases"] * 1.0 / (answers[i]["n_cases"] + answers[i]["n_controls"])
    
    settings["out_dir_group"] = settings["out_dir_group"] + "/" + base_last_dir
    settings["snp_list_file"] = "/users/mgloud/projects/coloc_comparisons/tmp/simulations/snp_list{0}.txt".format(i)

    settings["gwas_experiments"]["{1}/hg19/gwas/gwas_sumstats{0}.txt.gz".format(i, base_dir)] = \
            {"ref": "1kgenomes", "gwas_format": "effect_size", "type": "cc", "N": str(n_gwas), "s": str(cc_ratio_gwas)}
    settings["eqtl_experiments"]["{1}/hg19/eqtl/eqtl_sumstats{0}.txt.gz".format(i, base_dir)] = \
            {"ref": "1kgenomes", "N": str(n_eqtl), "ref": "ref-file", "rtc_genos": "{0}/hg19/eqtl/eqtl_genotypes{1}.vcf.gz".format(base_dir, i), "rtc_phenos": "{0}/hg19/eqtl/eqtl_phenotypes{1}.bed.gz".format(base_dir, i), "phenos": "{0}/hg19/eqtl/eqtl_phenotypes{1}.bed.gz".format(base_dir, i)}
    settings["ref_genomes"]["ref-file"]["file"] = "{1}/hg19/eqtl/eqtl_genotypes{0}.vcf.gz".format(i, base_dir)
    
    # Dump the config template to a file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/simulations/config{0}.config".format(i), "w") as w:
        json.dump(settings, w)

    with open("/users/mgloud/projects/coloc_comparisons/tmp/simulations/snp_list{0}.txt".format(i), "w") as w:
        w.write("{0}\t{1}\n".format(answers[i]["seed_chrom"], answers[i]["seed_pos"]))

    # Get it going
    subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/coloc_comparisons/tmp/simulations/config{0}.config 1".format(i), shell=True)

if __name__ == "__main__":
    main()
