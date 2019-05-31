#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 7/23/2018

import glob
import subprocess
import sys
from multiprocessing import Pool


def main():

    config_template = '''
    {{
            "rsid_index_file": "/users/mgloud/projects/index-dbsnp/data/hg19/common_all_20170710.vcf.gz",

            "out_dir_group": "smr-comparisons/{2}",

            "selection_basis": "snps_from_list",

            "snp_list_file":
                    "/users/mgloud/projects/coloc_comparisons/tmp/smr_snp_list{0}.txt",

            "selection_thresholds":
            {{
                "eqtl": 1,
                "gwas": 1
            }},

            "gwas_experiments": 
            {{
                "{1}/hg19/gwas/gwas_sumstats{0}.txt.gz": {{"ref": "1kgenomes", "gwas_format": "effect_size"}}
            }},
            
            "eqtl_experiments":	
            {{
                "{1}/hg19/eqtl/eqtl_sumstats{0}.txt.gz": {{"ref": "ref-file"}}
            }},

            "methods": 
            {{
                    "smr":{{}}
            }},

            "ref_genomes": 
            {{
                    "1kgenomes": 
                    {{
                            "file": 
                                    "/mnt/lab_data/montgomery/shared/1KG/ALL.chr{{0}}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",

                            "af_attribute": 
                                    "AF",

                            "N": 
                                    2504
                    }},
                    "ref-file":
                    {{
                        "file":
                            "{1}/hg19/eqtl/eqtl_genotypes{0}.vcf.gz", 

                        "phenos": "{1}/hg19/eqtl/eqtl_phenotypes{0}.bed.gz",

                        "af_attribute": "",

                        "N":
                            1000
                    }}
            }}

    }}
    '''

    base_dir = sys.argv[1]
    num_tests = int(sys.argv[2])
    max_threads = int(sys.argv[3])
    
    # For each simulated colocalization:

    if base_dir[-1] == "/":
        base_dir = base_dir[:-1]
    base_last_dir = base_dir.strip().split("/")[-1]

    subprocess.call("rm -rf /users/mgloud/projects/brain_gwas/output/smr-comparisons/{0}".format(base_last_dir), shell=True)

    answers = {}
    with open(base_dir + "/answer_key.txt") as f:
        f.readline()
        for line in f:
            data = line.strip().split("\t")
            answers[int(data[0])] = {}
            answers[int(data[0])]["seed_chrom"] = int(data[1])
            answers[int(data[0])]["seed_pos"] = int(data[2])

    pool = Pool(max_threads)
    for i in range(num_tests):        
        pool.apply_async(test_wrapper, args=(i, answers, config_template, base_dir, base_last_dir))
    pool.close()
    pool.join()
    

def test_wrapper(i, answers, config_template, base_dir, base_last_dir):
    try:
        run_test(i, answers, config_template, base_dir, base_last_dir)
    except Exception:
        traceback.print_exc(file=sys.stdout)


def run_test(i, answers, config_template, base_dir, base_last_dir):

        with open("/users/mgloud/projects/coloc_comparisons/tmp/smr_snp_list{0}.txt".format(i), "w") as w:
            w.write("{0}\t{1}\n".format(answers[i]["seed_chrom"], answers[i]["seed_pos"]))

        # Dump the config template to a file
        with open("/users/mgloud/projects/coloc_comparisons/tmp/smr{0}.config".format(i), "w") as w:
            w.write(config_template.format(i, base_dir, base_last_dir))

        # Get it going
        subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/coloc_comparisons/tmp/smr{0}.config 1".format(i), shell=True)


if __name__ == "__main__":
    main()
