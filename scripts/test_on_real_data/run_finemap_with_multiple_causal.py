import operator
import subprocess
import json
import sys
import time

# Define the subset of traits we actually want to run...
# Or maybe, just to make things easy for now...let's run all from the other test, but with blood only.

def main():
    # Reset things fresh on each run, so we're not mixing results
    subprocess.call("rm -rf /users/mgloud/projects/brain_gwas/output/finemap-multiple-causal/*", shell=True)
    subprocess.call("rm -rf /users/mgloud/projects/coloc_comparisons/tmp/*", shell=True)

    kept_data = []
    with open("/users//mgloud/projects/gwas/output/snps_to_test.txt") as f:
        all_data = []
        f.readline()
        for line in f:
            if "Adipose_Sub" not in line:
                continue
            data = line.strip().split()
            kept_data.append(data)

    kept_data = sorted(kept_data, key=operator.itemgetter(2))

    # Then for every locus in the "kept data"...
    for i in range(len(kept_data)):

        test = kept_data[i]
        print test
        
        if test[2].split("/")[-1] in ["GWAS_Adiponectin_Dastani_2012.txt.gz"]:
            continue

        temp = json.loads(template)
        temp["snp_list_file"] = "/users/mgloud/projects/gwas/tmp/snp_list{0}.txt".format(i)

        # Add locus to SNP list...but only once for each gene
        with open("/users/mgloud/projects/gwas/tmp/snp_list{0}.txt".format(i), "w") as w:
            w.write("{0}\t{1}\t{2}\n".format(test[0], test[1], test[7]))
               
        # Add corresponding gwas experiment to the list, if not already present
        temp["gwas_experiments"][test[2]] = {"ref": "1kgenomes", "gwas_format": "pval_only", "N": "50000", "type":"quant"}
        if test[2] != test[4]:
            temp["gwas_experiments"][test[2]]["traits"] = [test[4]]

        # Add corresponding eQTL tissue to the list
        temp["eqtl_experiments"][test[3]] = {"ref": "1kgenomes", "gwas_format": "effect_size", "N": "500"}

        # Write config file to the appropriate directory
        with open("/users/mgloud/projects/gwas/tmp/lc_config{0}.config".format(i), "w") as w:
            json.dump(temp, w)

        # Run the test
        subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/gwas/tmp/lc_config{0}.config 1 &".format(i), shell=True)

        while int(subprocess.check_output('''ps -ef | grep "python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/gwas/tmp/lc_config" | grep -v grep | wc -l''', shell=True)) > 8:
            time.sleep(1)

template = '''
{
        "out_dir_group": "finemap-multiple-causal",

       "gwas_experiments": 
	{
	},
	
	"eqtl_experiments":	
	{
	},

	"eqtl_threshold": 
		1,

	"selection_basis": 
		"snps_from_list",

	"snp_list_file":
                "/users/mgloud/projects/gwas/tmp/snp_list.txt",

	"plot_all":
                "True",

	"methods": 
	{
                "finemap": {"max_causal": "2", 
                            "save_all_finemap": "True"}
	},

        "ref_genomes": 
	{
		"1kgenomes": 
		{
			"file": 
				"/mnt/lab_data/montgomery/shared/1KG/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",

                	"af_attribute": 
				"AF",

                        "N": 
				2504
	        }
        }
}
'''

if __name__ == "__main__":
    main()
