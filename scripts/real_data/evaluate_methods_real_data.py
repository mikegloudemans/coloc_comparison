import operator
import subprocess
import json
import sys
import time
import datetime
import pandas as pd

def main():

    config_file = sys.argv[1]
    with open(config_file) as f:
        template = f.read()

    # Purge tmp files
    subprocess.call("rm -rf /users/mgloud/projects/coloc_comparisons/tmp/*", shell=True)
    
    # Create new subdirectory with time-stamp for this run
    now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S.%f')
    #out_dir = "/users/mgloud/projects/brain_gwas/output/coloc_comparisons/{0}".format(now)
    out_dir = "coloc_comparisons/{0}".format(now)
    subprocess.check_call("mkdir -p {0}".format(out_dir), shell=True)

    temp = json.loads(template)

    kept_data = []
    with open(temp["snp_list_file"]) as f:
        all_data = []
        f.readline()
        for line in f:
            if "Adipose_Sub" not in line:
                continue
            data = line.strip().split()
            kept_data.append(data)

    kept_data = sorted(kept_data, key=operator.itemgetter(2))

    # Special loading phase for RTC...
    if "rtc" in temp["coloc_settings"]["methods"]:

        # Load persistent data
        pheno_file = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz"
        id_file = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt"
        
        # Load in map of IDs 
        id_map = pd.read_csv(id_file, sep="\t")
        id_map['SMTSD'] = id_map['SMTSD'].apply(lambda x: x.replace(" - ", "_")).apply(lambda x: x.replace(" ", "_")).apply(lambda x: x.replace("(", "")).apply(lambda x: x.replace(")", ""))
            
        pheno_map = pd.read_csv(pheno_file, sep="\t", skiprows=2)

    # Then for every locus in the "kept data"...
    for i in range(len(kept_data)):

        test = kept_data[i]
        print test
        
        temp = json.loads(template)["coloc_settings"]
        temp["snp_list_file"] = "/users/mgloud/projects/coloc_comparisons/tmp/snp_list{0}.txt".format(i)
        temp["out_dir_group"] = out_dir

        if "smr" in temp["methods"]:
            # Needed for SMR
            # NOTE: We're doing this for now because the new munged data have beta/se
            # for some GWAS, but the old do not for any. Eventually, change this to make
            # it run on the exact same set of files
            test[2] = test[2].replace("munged", "munged/new")
            test[4] = test[4].replace("munged", "munged/new")     
   
        # This file has issues I guess, figure this out later
        if test[2].split("/")[-1] in ["GWAS_Adiponectin_Dastani_2012.txt.gz"]:
            continue

        if "Anthro" not in test[2]:
            continue

 
        # Add locus to SNP list...but only once for each gene
        with open("/users/mgloud/projects/coloc_comparisons/tmp/snp_list{0}.txt".format(i), "w") as w:
            # Chrom / SNP / gene
            w.write("{0}\t{1}\t{2}\n".format(test[0], test[1], test[7]))
               
        # Add corresponding gwas experiment to the list, if not already present
        temp["gwas_experiments"][test[2]] = {"ref": "1kgenomes", "gwas_format": "pval_only", "N": "50000", "type":"quant"}
        if test[2] != test[4]:
            temp["gwas_experiments"][test[2]]["traits"] = [test[4]]

        # Add corresponding eQTL tissue to the list
        temp["eqtl_experiments"][test[3]] = {"ref": "1kgenomes", "gwas_format": "effect_size", "N": "500"}

        # Prep for certain special cases for certain methods

        if "twas" in temp["methods"]:
            # Add corresponding eQTL tissue model
            temp["eqtl_experiments"][test[3]]["model"] = {
                    "pos":"/users/mgloud/projects/brain_gwas/scripts/data_prep/twas/Adipose_Visceral_Omentum.P01.pos", 
                    "gene":"/users/mgloud/projects/brain_gwas/scripts/data_prep/twas/Adipose_Visceral_Omentum/Adipose_Visceral_Omentum.{0}.wgt.RDat", 
                    "gene_dir_name": "Adipose_Visceral_Omentum"
            }

        if "rtc" in temp["methods"]:
            temp["eqtl_experiments"][test[3]]["rtc_genos"] = "/users/mgloud/projects/coloc_comparisons/tmp/genos{0}.vcf.gz".format(i)
            temp["eqtl_experiments"][test[3]]["rtc_phenos"] = "/users/mgloud/projects/coloc_comparisons/tmp/phenos{0}.bed.gz".format(i)

            prep_rtc_input(id_map, pheno_map, test, i)

        # Write config file to the appropriate directory
        with open("/users/mgloud/projects/coloc_comparisons/tmp/coloc_comparisons{0}.config".format(i), "w") as w:
            json.dump(temp, w)

        # TODO: Turn this into not a BG but a foreground task, then have it called by a multi-threading pool
        # Run the test
        subprocess.call("python /users/mgloud/projects/brain_gwas/scripts/dispatch.py /users/mgloud/projects/coloc_comparisons/tmp/coloc_comparisons{0}.config 1".format(i), shell=True)

    # Purge tmp files
    subprocess.call("rm -rf /users/mgloud/projects/coloc_comparisons/tmp/*", shell=True)

def prep_rtc_input(id_map, pheno_map, test, i):
        # Do preprocessing to prepare the additional input for RTC! (Note: eventually
        # we should make this into a separate optional script, where the required input
        # files are specified and then the script pulls out the necessary subset of info)
       
        geno_file = "/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz"

        # Figure out which tissue we're going for
        tissue = test[3].split("/")[-1].split(".")[0]

       
        # Get all samples corresponding to our tissue of interest
        id_map = id_map[(id_map['SMTSD'] == tissue).values]
        tissue_samples = list(id_map['SAMPID'])

        # Subset genotype file down to region of interest
        # Header
        subprocess.check_call("zcat {0} 2> /dev/null | head -n 10000 | grep \#CHROM > /users/mgloud/projects/coloc_comparisons/tmp/genotype_vcf{1}.tmp".format(geno_file, i), shell=True)
        subprocess.check_call("zcat {0} 2> /dev/null | head -n 10000 | grep \## > /users/mgloud/projects/coloc_comparisons/tmp/genos{1}.vcf".format(geno_file, i), shell=True)
        # Everything else we need
        subprocess.check_call('''tabix {0} {1}:{2}-{3} | awk '{{split($3, a, "_");  $3=a[1] "_" a[2]; print $0}}' | sed s/\\ /\\\\t/g >> /users/mgloud/projects/coloc_comparisons/tmp/genotype_vcf{4}.tmp'''.format(geno_file, test[0], int(test[1]) - 1000000, int(test[1]) + 1000000, i), shell=True)


        # Load subsetted genotype (No filtering is actually needed for the
        # individuals in the genotype matrix, because QTLtools performs this filtering)
        genotypes = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/genotype_vcf{0}.tmp".format(i), sep="\t")
        genotypes["#CHROM"] = genotypes["#CHROM"].apply(lambda x: "chr" + str(x))

        # Filter expression matrix down to a single gene of interest
        # and to samples we want
        pheno_sub = pheno_map[pheno_map['Name'] == test[7]]
        tissue_subset = ["Name"] + [e for e in list(pheno_sub.columns.values) if e in tissue_samples]
        expression = pheno_sub[tissue_subset]

        # NOTE: The expression matrix from GTEx does not contain values matched to a normal distribution. However,
        # since QTLtools does this automatically, we don't have to worry about formatting this ourselves.

        # Change column headers in the expression matrix to include only the donor ID, not the full sample ID
        expression = expression.rename(lambda x: "-".join(x.split("-")[:2]), axis='columns')
        # Add additional necessary headers
        # NOTE: "chr" trick here is kind of hackish and may cause bugs, try to fix this
        expression["#chr"] = "chr" + str(test[0])
        # TODO: This might not be appropriate, but for now we're doing it to make things work
        expression["start"] = int(test[1]) - 1000000
        expression["end"] = int(test[1]) + 1000000
        expression = expression.rename(index=str, columns={"Name": "gene"})
        expression["length"] = 2000000
        expression["strand"] = "+"
        column_ids = ["#chr", "start", "end", "gene", "length", "strand"]
        column_ids += [s for s in list(expression.columns.values) if "GTEX" in s]
        expression = expression[column_ids]

        # Rewrite genotype matrix, bgzip and tabix
        with open("/users/mgloud/projects/coloc_comparisons/tmp/genos{0}.vcf".format(i), "a") as a:
            genotypes.to_csv(a, sep="\t", index=False)
        subprocess.check_call("bgzip -f /users/mgloud/projects/coloc_comparisons/tmp/genos{0}.vcf".format(i), shell=True)
        subprocess.check_call("tabix -f -S 2 -s 1 -b 2 -e 2 /users/mgloud/projects/coloc_comparisons/tmp/genos{0}.vcf.gz".format(i), shell=True)

        # Rewrite expression matrix and bgzip and tabix
        expression.to_csv("/users/mgloud/projects/coloc_comparisons/tmp/phenos{0}.bed".format(i), sep="\t", index=False)
        subprocess.check_call("bgzip -f /users/mgloud/projects/coloc_comparisons/tmp/phenos{0}.bed".format(i), shell=True)
        subprocess.check_call("tabix -f -S 1 -s 1 -b 2 -e 3 /users/mgloud/projects/coloc_comparisons/tmp/phenos{0}.bed.gz".format(i), shell=True)


if __name__ == "__main__":
    main()
