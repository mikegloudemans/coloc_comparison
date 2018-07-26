#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 7/23/2018

import glob
import subprocess
import sys

base_dir = sys.argv[1]
num_tests = int(sys.argv[2])

config_template = '''
{{
        "out_dir_group": "coloc-comparisons/{2}",

        "selection_basis": "snps_from_list",

	"snp_list_file":
                "/users/mgloud/projects/coloc_comparisons/tmp/snp_list.txt",

        "selection_thresholds":
        {{
            "eqtl": 1,
            "gwas": 1
        }},

        "gwas_experiments": 
	{{
            "{1}/hg19/gwas/gwas_sumstats{0}.txt.gz": {{"ref": "1kgenomes", "gwas_format": "effect_size", "type": "cc", "N": {3}, "s": "{4}"}}
	}},
	
	"eqtl_experiments":	
	{{
            "{1}/hg19/eqtl/eqtl_sumstats{0}.txt.gz": {{"ref": "1kgenomes", "N": {5}}}
	}},

	"methods": 
	{{
		"coloc":{{}}
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
	        }}
        }}
}}
'''

# For each simulated colocalization:

if base_dir[-1] == "/":
    base_dir = base_dir[:-1]
base_last_dir = base_dir.strip().split("/")[-1]

subprocess.call("rm -rf /users/mgloud/projects/brain_gwas/output/coloc-comparisons/{0}".format(base_last_dir), shell=True)

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

for i in range(num_tests):
    
    with open("/users/mgloud/projects/coloc_comparisons/tmp/snp_list.txt", "w") as w:
        w.write("{0}\t{1}\n".format(answers[i]["seed_chrom"], answers[i]["seed_pos"]))

    n_gwas = answers[i]["n_cases"] + answers[i]["n_controls"]
    cc_ratio_gwas = answers[i]["n_cases"] * 1.0 / (answers[i]["n_cases"] + answers[i]["n_controls"])
    n_eqtl = answers[i]["n_eqtl"]

    # Dump the config template to a file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/coloc.config", "w") as w:
        w.write(config_template.format(i, base_dir, base_last_dir, n_gwas, cc_ratio_gwas, n_eqtl))

    # Get it going
    subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/coloc_comparisons/tmp/coloc.config", shell=True)

# (Later: Separate script)

# Parse eCAVIAR results file to see what we find

# Compile all results into a sort of ROC curve;
# probably want to do this in R

# Can also stratify results based on the magnitude of
# effect sizes and stuff, see how much we're sensitive to
# parameters like that
