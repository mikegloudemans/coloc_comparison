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
        "out_dir_group": "ecaviar-comparisons/{2}",

        "selection_basis": "gwas",

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
            "{1}/hg19/eqtl/eqtl_sumstats{0}.txt.gz": {{"ref": "1kgenomes"}}
	}},

	"methods": 
	{{
		"finemap":{{}}
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

for i in range(num_tests):
#for i in [9]:
    
    # Dump the config template to a file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/ecaviar.config", "w") as w:
        w.write(config_template.format(i, base_dir, base_last_dir))

    # Get it going
    subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/coloc_comparisons/tmp/ecaviar.config", shell=True)

# (Later: Separate script)

# Parse eCAVIAR results file to see what we find

# Compile all results into a sort of ROC curve;
# probably want to do this in R

# Can also stratify results based on the magnitude of
# effect sizes and stuff, see how much we're sensitive to
# parameters like that
