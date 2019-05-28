#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 7/23/2018

import glob
import subprocess
import sys
import time

base_dir = sys.argv[1]
num_tests = int(sys.argv[2])

config_template = '''
{{
        "out_dir_group": "baseline-comparisons/{2}",

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
		"baseline":{{}}
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

subprocess.call("rm -rf /users/mgloud/projects/brain_gwas/output/baseline-comparisons/{0}".format(base_last_dir), shell=True)

for i in range(num_tests):
   
    # Dump the config template to a file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/baseline{0}.config".format(i), "w") as w:
        w.write(config_template.format(i, base_dir, base_last_dir))

    # Get it going
    subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/coloc_comparisons/tmp/baseline{0}.config 1 &".format(i), shell=True)
    
    while int(subprocess.check_output('''ps -ef | grep "python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/coloc_comparisons/tmp/baseline" | wc -l''', shell=True)) > 30:
        time.sleep(0.01)

