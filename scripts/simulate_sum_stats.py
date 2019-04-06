#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 7/16/2018

import math
import config
import random
import sys
import numpy as np
import gzip
import operator
import pandas as pd
import subprocess
if sys.version_info[0] < 3: 
   from StringIO import StringIO
else:
   from io import StringIO
from scipy import stats
from shutil import copyfile
import datetime

config_file = sys.argv[1]
settings = config.load_config(config_file)

now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
if "out_dir_group" in settings:
    base_output_dir = "/users/mgloud/projects/coloc_comparisons/output/simulations/{0}/{1}".format(settings["out_dir_group"], now)
else:
    base_output_dir = "/users/mgloud/projects/coloc_comparisons/output/simulations/{0}".format(now)
subprocess.call("mkdir -p {0}/hg19/gwas".format(base_output_dir), shell=True)
subprocess.call("mkdir -p {0}/hg19/eqtl".format(base_output_dir), shell=True)
subprocess.call("mkdir -p {0}/hg38/gwas".format(base_output_dir), shell=True)
subprocess.call("mkdir -p {0}/hg38/eqtl".format(base_output_dir), shell=True)

copyfile(config_file, "{0}/settings_used.config".format(base_output_dir))

def main():
    with open("{0}/answer_key.txt".format(base_output_dir), "w") as w:
        w.write("test_case\tseed_chrom\tseed_pos\tcausal_variants\tcases_n\tcontrols_n\teqtl_n\tsnps_n\n")

    # Get possible GWAS loci upfront, so we don't
    # waste any more time on this
    locus_bank = get_possible_loci(settings)

    for i in range(settings["total_test_sites"]):

        random.seed(i)
        
        # Loop until we get a variant that works
        while True:

            # Store settings for current run, for easy reference
            settings["current_run"] = {}

            print "Starting new locus..."
            # Choose a locus 
            locus = random.choice(locus_bank)
            print locus

            # Figure out if there's going to be a causal GWAS variant
            gwas_effect_size = get_gwas_effect_size(settings)

            # Get all haplotypes using HAPGEN2
            gwas_effect_sizes = run_hapgen2(settings, locus, gwas_effect_size)
            settings["current_run"]["gwas_effect_sizes"] = gwas_effect_sizes
            if gwas_effect_sizes == "Bad variant" or gwas_effect_sizes == "Fail":
                # Make sure HAPGEN2 actually succeeded, i.e. it was a valid site
                continue

            # Get effect sizes, using LD pairings at restrictions
            eqtl_effect_sizes = get_eqtl_effect_sizes(settings, gwas_effect_sizes)
            settings["current_run"]["eqtl_effect_sizes"] = eqtl_effect_sizes
            if eqtl_effect_sizes == "Impossible LD matrix":
                continue

            eqtl_phenotypes = get_expression_phenotypes(eqtl_effect_sizes)
            if "simulate_peer" in settings and settings["simulate_peer"] == "True":
                peer_factors = simulate_peer_factors(eqtl_phenotypes, i)

            # Finally, create summary statistics
            sum_stats = make_sum_stats(eqtl_phenotypes)
            (eqtls, gwas) = sum_stats

            if "simulate_peer" in settings and settings["simulate_peer"] == "True":
                peer_eqtls = call_eqtls_with_peer(eqtl_phenotypes, peer_factors)
                write_peer_sumstats(peer_eqtls, i)

            write_sumstats(eqtls, gwas, i)
            write_answer_key(gwas_effect_sizes, eqtl_effect_sizes, locus, i)
            write_expression_phenotypes(eqtl_phenotypes, i, locus)
            if not "no_genotypes" in settings or settings["no_genotypes"] == "False":
                write_genotypes_as_vcf(i, locus)

            break

    # liftOver to hg38
    if not "liftover" in settings or settings["liftover"] == "True":
        run_liftover(settings)

def get_eqtl_effect_sizes(settings, gwas_effect_sizes):

    eqtl_effect_sizes = [0] * len(gwas_effect_sizes)
    potential_eqtl_effect = random.uniform(settings["eqtl_min_effect_size"], settings["eqtl_max_effect_size"])
    if max([abs(ges) for ges in gwas_effect_sizes]) > 0:
        # GWAS is causal
        gwas_causal_indices = [i for i,ges in enumerate(gwas_effect_sizes) if abs(ges) > 0]

        if random.random() < settings["p_eqtl_causal_given_gwas_causal"]:
            # eQTL is causal

            if random.random() < settings["p_same_causal_given_both_causal"]:
                # Same causal variant for both
                eqtl_effect_sizes[gwas_causal_indices[0]] = potential_eqtl_effect
            else:
                # Different eQTL causal variant
                causal_eqtl_index = random.randint(len(eqtl_effect_sizes) * 1 / 4, len(eqtl_effect_sizes) * 3 / 4)
                eqtl_effect_sizes[causal_eqtl_index] = potential_eqtl_effect
    else:
        # GWAS is not causal
        if random.random() < settings["p_eqtl_causal_given_no_gwas_causal"]:
            # eQTL is causal
            # Choose causal eQTL; don't let it be right at the edge of the measured region though
            causal_eqtl_index = random.randint(len(gwas_effect_sizes) * 1 / 4, len(gwas_effect_sizes) * 3 / 4)
            eqtl_effect_sizes[causal_eqtl_index] = potential_eqtl_effect

    return eqtl_effect_sizes

# Get the set of loci that are hits in the given GWAS files
def get_possible_loci(settings):
    possible_loci = []
    for gwas_file in settings["gwas_reference_files"]:

        with gzip.open(gwas_file) as f:
            header = f.readline().strip().split()

            pval_index = header.index("pvalue")
            chr_index = header.index("chr")
            snp_pos_index = header.index("snp_pos")

            all_snps = []

            for line in f:
                data = line.strip().split("\t")
                try:
                    pvalue = float(data[pval_index])
                except:
                    continue
                chr = data[chr_index]
                pos = int(data[snp_pos_index])
                if pvalue > settings["gwas_threshold"]:
                    continue

                all_snps.append((chr, pos, pvalue))
       
        # For now, include only autosomal SNPs.
        filtered = []
        for s in all_snps:
            if "chr" in str(s[0]):
                filtered.append((s[0][3:], s[1], s[2]))
            else:
                filtered.append((s[0], s[1], s[2]))

        all_snps = sorted(filtered, key=operator.itemgetter(2)) 
        
        # Go through the list of SNPs in order, adding the ones
        # passing our criteria.
        for snp in all_snps:

            # For now, ignore a SNP if it's in the MHC region -- this
            # would require alternative methods.
            if (snp[0] == "6") and snp[1] > 25000000 and snp[1] < 35000000:
                    continue

            # Before adding a SNP, make sure it's not right next
            # to another SNP that we've already selected.
            skip = False
            for kept_snp in possible_loci:
                    if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < settings["window"]:
                            skip = True
                            break
            if not skip:
                possible_loci.append(snp)

    return possible_loci

def run_hapgen2(settings, locus, gwas_effect_size):

    # Get the data from 1000 Genomes VCF near locus
    filename = "/mnt/lab_data/montgomery/shared/1KG/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz".format(locus[0])
    header = subprocess.check_output("zcat {0} 2> /dev/null | head -n 500 | grep \\#CHROM".format(filename), shell=True).strip().split()
    stream = StringIO(subprocess.check_output("tabix {0} {1}:{2}-{3}".format(filename, locus[0], locus[1] - settings["window"], locus[1] + settings["window"]), shell=True))
    vcf = pd.read_csv(stream, sep="\t", names=header)

    # Remove variants with position appearing multiple times
    dup_counts = {}
    for pos in vcf["POS"]:
            dup_counts[pos] = dup_counts.get(pos, 0) + 1
    vcf["dup_counts"] = [dup_counts[pos] for pos in vcf['POS']]
    vcf = vcf[vcf["dup_counts"] == 1]

    # Filter down to a min MAF of 0.05 for now
    def maf(x):
        af = x.split("AF=")[1].split(";")[0]
        if "," in af:
            return 0
        af = float(af)
        maf = min(abs(1-af), abs(af))
        return maf
    vcf['MAF'] = vcf['INFO'].apply(maf)
    vcf = vcf[vcf['MAF'] > 0.05]

    # Make sure the variant of interest hasn't been filtered out;
    # if it has, then get a new one.
    if int(locus[1]) not in set(list(vcf['POS'])):
        return "Bad variant"
    
    # Format it for HAPGEN
    # Write .haps file
    vcf_genos = vcf.iloc[:,10:(vcf.shape[1]-2)]
    with open("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2.haps", "w") as w:
        for i, row in vcf_genos.iterrows():
            w.write("|".join(list(row)).replace("|", " ") + "\n")

    # Write .leg file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2.leg", "w") as w:
        w.write("rs position X0 X1\n")
        for i, row in vcf.iterrows():
            w.write(" ".join([str(r) for r in [row['ID'], row['POS'], row['REF'], row['ALT']]]) + "\n")

    # Write .map centimorgans recombination rate file
    map_file = "/users/mgloud/projects/coloc_comparisons/data/genetic-map/extended_genetic_map_GRCh37_chr{0}.txt.sorted.gz".format(locus[0])
    head = subprocess.check_output("zcat {0} | head -n 1".format(map_file), shell=True)
    data = subprocess.check_output("tabix {0} {1}:{2}-{3} | sed s/\\t/\ /g | cut -f2-4 -d ' '".format(map_file, locus[0], locus[1] - settings["window"], locus[1] + settings["window"]), shell=True).split("\n")
    with open("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2.map", "w") as w:
        w.write(head)
        for line in data:
            if line == "":
                continue
            info = line.split()
            if int(info[1]) in set(list(vcf['POS'])):
                w.write(line + "\n")

    gwas_control_sample_size = int(math.floor(10**random.uniform(settings["gwas_min_control_log_sample_size"], settings["gwas_max_control_log_sample_size"])))
    gwas_case_sample_size = int(math.floor(10**random.uniform(settings["gwas_min_case_log_sample_size"], settings["gwas_max_case_log_sample_size"])))
    eqtl_sample_size = int(math.floor(10**random.uniform(settings["eqtl_min_log_sample_size"], settings["eqtl_max_log_sample_size"])))
    settings["current_run"]["gwas_control_sample_size"] = gwas_control_sample_size
    settings["current_run"]["gwas_case_sample_size"] = gwas_case_sample_size
    settings["current_run"]["eqtl_sample_size"] = eqtl_sample_size

    # Then run HAPGEN2 to get genotypes
    # NOTE: This could be dangerous because it could lead to infinite loops
    try:
        subprocess.check_call("hapgen2 -m /users/mgloud/projects/coloc_comparisons/tmp/hapgen2.map -l /users/mgloud/projects/coloc_comparisons/tmp/hapgen2.leg -h /users/mgloud/projects/coloc_comparisons/tmp/hapgen2.haps -o /users/mgloud/projects/coloc_comparisons/tmp/hapgen2_gwas.out -dl {0} 1 {1} {2} -n {3} {4}".format(locus[1], 10**gwas_effect_size, 10**(2*gwas_effect_size), gwas_control_sample_size, gwas_case_sample_size), shell=True)
        subprocess.check_call("hapgen2 -m /users/mgloud/projects/coloc_comparisons/tmp/hapgen2.map -l /users/mgloud/projects/coloc_comparisons/tmp/hapgen2.leg -h /users/mgloud/projects/coloc_comparisons/tmp/hapgen2.haps -o /users/mgloud/projects/coloc_comparisons/tmp/hapgen2_eqtl.out -dl {0} 1 1 1 -n {1} 1".format(locus[1], eqtl_sample_size), shell=True)
    except:
        return "Fail"

    # Write VCF to tmp file in case we need to compute LD
    vcf.to_csv('/users/mgloud/projects/coloc_comparisons/tmp/tmp.vcf', sep="\t", index=False, header=True)
    settings["current_run"]["rsids"] = vcf['ID']

    gwas_effect_sizes = [0] * vcf.shape[0]
    gwas_effect_sizes[list(vcf['POS']).index(locus[1])] = gwas_effect_size

    return gwas_effect_sizes
    
def get_gwas_effect_size(settings):
  
    if random.random() < settings["p_gwas_causal"]:
        return random.uniform(settings["gwas_min_effect_size"], settings["gwas_max_effect_size"])
    else:
        return 0

def get_expression_phenotypes(eqtl_effect_sizes):

    # Read haplotypes for each individual
    haps = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
    genos = np.add(haps.iloc[:, range(0,haps.shape[1]-1,2)].as_matrix(), haps.iloc[:, range(1,haps.shape[1]-1,2)].as_matrix())

    means = np.matmul(genos.T, eqtl_effect_sizes)
    phenos = np.random.normal(means, 1, len(means))

    return phenos

def simulate_peer_factors(phenos, index):

    all_peer_factors = []

    for frac in settings["peer_fraction_var_explained"]:
        pheno_var = np.var(phenos)
        total_var = pheno_var / frac
        peer_factors = np.random.normal(phenos, total_var-pheno_var, len(phenos))

        with open("{0}/hg19/eqtl/peer_factors{1}_{2}.txt".format(base_output_dir, index, frac), "w") as w:
            header = "\t".join(["ID" + str(i) for i in range(len(phenos))]) + "\n"
            w.write(header)
            w.write("\t".join([str(pf) for pf in peer_factors]))
            w.write("\n")

        all_peer_factors.append(peer_factors)

    return all_peer_factors

def write_expression_phenotypes(phenos, index, locus):
    # We do this because RTC needs to re-call eQTLs with various SNPs regressed out
    with open("{0}/hg19/eqtl/eqtl_phenotypes{1}.bed".format(base_output_dir, index), "w") as w:

        header = "#chr\tstart\tend\tgene\tlength\tstrand\t" + "\t".join(["ID" + str(i) for i in range(len(phenos))]) + "\n"
        w.write(header)

        # Now put in the phenotypes for our fake gene
        w.write("chr{0}\t{1}\t{2}\tnameless_gene\t500\t+\t".format(locus[0], locus[1], int(locus[1]) + 500))
        w.write("\t".join([str(p) for p in phenos]))
        w.write("\n")

    # bgzip and tabix it
    subprocess.check_call("bgzip -f {0}/hg19/eqtl/eqtl_phenotypes{1}.bed".format(base_output_dir, index), shell=True)
    subprocess.check_call("tabix -f -S 1 -s 1 -b 2 -e 3 {0}/hg19/eqtl/eqtl_phenotypes{1}.bed.gz".format(base_output_dir, index), shell=True)

def write_genotypes_as_vcf(index, locus):
    
    metadata = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2_eqtl.out.controls.gen", sep=" ", header=None)
    haps = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
    
    # We do this because RTC needs to re-call eQTLs with various SNPs regressed out
    with open("{0}/hg19/eqtl/eqtl_genotypes{1}.vcf".format(base_output_dir, index), "w") as w:

        w.write("##fileformat=VCFv4.1\n") # This is essential; otherwise QTLtools will not work
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(["ID" + str(i) for i in range(haps.shape[1] / 2)]) + "\n"
        w.write(header)

        # Now put in our simulated genotypes

        # NOTE: For large sample sizes this can be incredibly slow.
        # Should find a faster way, or limit this part to eQTLs w/ smaller sample sizes.
        for i in range(metadata.shape[0]):
            w.write("chr{0}\t{1}\t{0}_{1}\t{2}\t{3}\t100\tPASS\tNONE\tGT".format(locus[0], metadata.iloc[i,2], metadata.iloc[i,3], metadata.iloc[i,4]))
            for j in range(haps.shape[1] / 2):
                w.write("\t" + "|".join([str(int(s)) for s in list(haps.iloc[i,2*j:(2*(j+1))])]))
            w.write("\n")

    # bgzip and tabix it
    subprocess.check_call("bgzip -f {0}/hg19/eqtl/eqtl_genotypes{1}.vcf".format(base_output_dir, index), shell=True)
    subprocess.check_call("tabix -f -S 1 -s 1 -b 2 -e 2 {0}/hg19/eqtl/eqtl_genotypes{1}.vcf.gz".format(base_output_dir, index), shell=True)

def make_sum_stats(eqtl_phenotypes):
    eqtls = call_eqtls(eqtl_phenotypes)
    gwas = run_gwas()
    return (eqtls, gwas)

def call_eqtls(eqtl_phenotypes):
    # Load genotypes
    # Read haplotypes for each individual
    haps = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
    genos = np.add(haps.iloc[:, range(0,haps.shape[1]-1,2)].as_matrix(), haps.iloc[:, range(1,haps.shape[1]-1,2)].as_matrix())

    # Call eQTLs
    eqtls = []
    for i in range(genos.shape[0]):
        predictors = genos[i,:]
        slope, intercept, r_value, p_value, std_err = stats.linregress(predictors, eqtl_phenotypes)
        if np.isnan(slope):
            slope = 0
            std_err = 1
        eqtls.append((slope, std_err, p_value))
    return eqtls

def call_eqtls_with_peer(eqtl_phenotypes, peer):

    # Load genotypes
    # Read haplotypes for each individual
    haps = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
    genos = np.add(haps.iloc[:, range(0,haps.shape[1]-1,2)].as_matrix(), haps.iloc[:, range(1,haps.shape[1]-1,2)].as_matrix())

    all_eqtls = []

    for p in peer:
        # Get residuals for gene expression
        slope, intercept, r_value, p_value, std_err = stats.linregress(p, eqtl_phenotypes)
        preds = p * slope + intercept
        residuals = eqtl_phenotypes - preds

        # Call eQTLs
        eqtls = []
        for i in range(genos.shape[0]):
            predictors = genos[i,:]
            slope, intercept, r_value, p_value, std_err = stats.linregress(predictors, residuals)
            if np.isnan(slope):
                slope = 0
                std_err = 1
            eqtls.append((slope, std_err, p_value))

        all_eqtls.append(eqtls)

    return all_eqtls



def run_gwas():
    # Load genotypes
    # Read haplotypes for each individual
    case_haps = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2_gwas.out.cases.haps", sep=" ", header=None)
    case_genos = np.add(case_haps.iloc[:, range(0,case_haps.shape[1]-1,2)].as_matrix(), case_haps.iloc[:, range(1,case_haps.shape[1]-1,2)].as_matrix())
    control_haps = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/hapgen2_gwas.out.controls.haps", sep=" ", header=None)
    control_genos = np.add(control_haps.iloc[:, range(0,control_haps.shape[1]-1,2)].as_matrix(), control_haps.iloc[:, range(1,control_haps.shape[1]-1,2)].as_matrix())

    # Run GWAS (for now, maybe just do it as a linear regression too? Makes it easily adaptable for continous phenotypes)
    gwas = []
    for i in range(case_genos.shape[0]):
        predictors = list(case_genos[i,:]) + list(control_genos[i,:])
        gwas_response = [1] * case_genos.shape[1] + [0] * control_genos.shape[1]
        slope, intercept, r_value, p_value, std_err = stats.linregress(predictors, gwas_response)
        # Catch issues in cases where one allele never occurred
        if np.isnan(slope):
            slope = 0
            std_err = 1
        gwas.append((slope, std_err, p_value))
    return gwas

def write_sumstats(eqtls, gwas, index):
    #header = subprocess.check_output("cat /users/mgloud/projects/coloc_comparisons/tmp/tmp.vcf 2> /dev/null | grep \\#CHROM", shell=True).strip().split()
    vcf = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/tmp.vcf", sep="\t")

    # Filter down to a min MAF of 0.05 for now
    def alt_af(x):
        af = x.split("AF=")[1].split(";")[0]
        if "," in af:
            return 0
        af = float(af)
        return(af)
    vcf['AF'] = vcf['INFO'].apply(alt_af)
 
    with open("{0}/hg19/gwas/gwas_sumstats{1}.txt".format(base_output_dir, index), "w") as w:
        w.write("rsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\tbeta\tse\tzscore\tpvalue\tn_cases\tn_controls\n")
        for i in range(vcf.shape[0]):
            var = vcf.iloc[i, :]
            g = gwas[i]
            variant_id = "{0}_{1}_{2}_{3}_hg19".format(var['#CHROM'], var["POS"], var["REF"], var["ALT"])
            w.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\n".format(var['ID'], variant_id, "hg19", var['#CHROM'], var["POS"], var["REF"], var["ALT"], var["AF"], g[0], g[1], g[0] / g[1], g[2], settings["current_run"]["gwas_case_sample_size"], settings["current_run"]["gwas_control_sample_size"]))
    
    with open("{0}/hg19/eqtl/eqtl_sumstats{1}.txt".format(base_output_dir, index), "w") as w:
        w.write("feature\trsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\tbeta\tse\tzscore\tpvalue\tN\n")
        for i in range(vcf.shape[0]):
            var = vcf.iloc[i, :]
            e = eqtls[i]
            variant_id = "{0}_{1}_{2}_{3}_hg19".format(var['#CHROM'], var["POS"], var["REF"], var["ALT"])
            w.write("nameless_gene\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(var['ID'], variant_id, "hg19", var['#CHROM'], var["POS"], var["REF"], var["ALT"], var["AF"], e[0], e[1], e[0] / e[1], e[2], settings["current_run"]["eqtl_sample_size"]))

    # bgzip and tabix
    subprocess.check_call("bgzip -f {0}/hg19/gwas/gwas_sumstats{1}.txt".format(base_output_dir, index), shell=True)
    subprocess.check_call("bgzip -f {0}/hg19/eqtl/eqtl_sumstats{1}.txt".format(base_output_dir, index), shell=True)
    subprocess.check_call("tabix -f -S 1 -s 4 -b 5 -e 5 {0}/hg19/gwas/gwas_sumstats{1}.txt.gz".format(base_output_dir, index), shell=True)
    subprocess.check_call("tabix -f -S 1 -s 5 -b 6 -e 6 {0}/hg19/eqtl/eqtl_sumstats{1}.txt.gz".format(base_output_dir, index), shell=True)

def write_peer_sumstats(eqtls, index):
    vcf = pd.read_csv("/users/mgloud/projects/coloc_comparisons/tmp/tmp.vcf", sep="\t")

    # Filter down to a min MAF of 0.05 for now
    def alt_af(x):
        af = x.split("AF=")[1].split(";")[0]
        if "," in af:
            return 0
        af = float(af)
        return(af)
    vcf['AF'] = vcf['INFO'].apply(alt_af)
 
    for eqtl_index in range(len(eqtls)): 

        with open("{0}/hg19/eqtl/eqtl_peer_sumstats{1}_{2}.txt".format(base_output_dir, index, settings["peer_fraction_var_explained"][eqtl_index]), "w") as w:
            w.write("feature\trsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\tbeta\tse\tzscore\tpvalue\tN\n")
            for i in range(vcf.shape[0]):
                var = vcf.iloc[i, :]
                e = eqtls[eqtl_index][i]
                variant_id = "{0}_{1}_{2}_{3}_hg19".format(var['#CHROM'], var["POS"], var["REF"], var["ALT"])
                w.write("nameless_gene\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\n".format(var['ID'], variant_id, "hg19", var['#CHROM'], var["POS"], var["REF"], var["ALT"], var["AF"], e[0], e[1], e[0] / e[1], e[2], settings["current_run"]["eqtl_sample_size"]))

        # bgzip and tabix
        subprocess.check_call("bgzip -f {0}/hg19/eqtl/eqtl_peer_sumstats{1}_{2}.txt".format(base_output_dir, index, settings["peer_fraction_var_explained"][eqtl_index]), shell=True)
        subprocess.check_call("tabix -f -S 1 -s 5 -b 6 -e 6 {0}/hg19/eqtl/eqtl_peer_sumstats{1}_{2}.txt.gz".format(base_output_dir, index, settings["peer_fraction_var_explained"][eqtl_index]), shell=True)


def write_answer_key(gwas_effect_sizes, eqtl_effect_sizes, locus, index):
    with open("{0}/answer_key.txt".format(base_output_dir), "a") as a:
        # Write GWAS/eQTL variants, effect sizes...might also add more info later
        info = ""
        for i in range(len(gwas_effect_sizes)):
            if gwas_effect_sizes[i] != 0:
                info += "gwas:" + settings["current_run"]["rsids"].iloc[i] + ":" + str(gwas_effect_sizes[i]) + ","
        for i in range(len(eqtl_effect_sizes)):
            if eqtl_effect_sizes[i] != 0:
                info += "eqtl:" + settings["current_run"]["rsids"].iloc[i] + ":" + str(eqtl_effect_sizes[i]) + ","
        a.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(index, locus[0], locus[1], info, settings["current_run"]["gwas_case_sample_size"], settings["current_run"]["gwas_control_sample_size"], settings["current_run"]["eqtl_sample_size"], len(gwas_effect_sizes)))

def run_liftover(settings):
    subprocess.check_call("python /users/mgloud/projects/coloc_comparisons/scripts/liftover_sumstats_hg19_to_hg38.py {0} {1}".format(base_output_dir, settings["total_test_sites"]), shell=True)


if __name__ == "__main__":
    main()
