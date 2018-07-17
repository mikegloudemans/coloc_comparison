#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 7/16/2018

import config
import random
import sys
import numpy as np

config_file = sys.argv[1]
settings = config.load_config(config_file)

def main():
    for i in range(settings["total_test_sites"]):

        # Choose a locus (perhaps based on where an existing GWAS hit
        # is known to be)
        #locus = choose_locus(settings)

        # Get genotypes using HAPGEN2
        #genotypes = run_hapgen2(settings)

        # Compute LD for those genotypes
        #ld = compute_ld(genotypes)
        ld = np.random.random(size=(100,100))
        ld = (ld + ld.T) / 2

        # Get effect sizes, using LD pairings at restrictions
        effect_sizes = get_effect_sizes(ld, settings)
        print zip(effect_sizes[0], effect_sizes[1])

        # Finally, create summary statistics
        #sum_stats = make_sum_stats(effect_sizes, genotypes, settings)

        # Output summary statistics to a file
        out_file = "TODO"
        #write_sum_stats(genotypes, out_file)


# Function: Get effect sizes
# Input: 
#   - A matrix showing LD between all positions
#   - Settings to specify probabilities of different configurations
# Output:
#   - A tuple (gwas_effect_sizes, eqtl_effect_sizes)
#     - Two vectors with effect sizes for the 
#
# The variants with non-zero effects, if any, are the
# causal ones. If eQTL and GWAS have same variants, then
# 

def get_effect_sizes(ld_matrix, settings):

    # Make sure LD matrix is square
    assert ld_matrix.shape[0] == ld_matrix.shape[1]
    assert np.array_equal(ld_matrix, ld_matrix.T)

    gwas_effect_sizes = [0] * ld_matrix.shape[0]
    eqtl_effect_sizes = [0] * ld_matrix.shape[0]

    potential_gwas_effect = random.uniform(settings["gwas_min_fraction_variance_explained"], settings["gwas_max_fraction_variance_explained"])
    potential_eqtl_effect = random.uniform(settings["eqtl_min_fraction_variance_explained"], settings["eqtl_max_fraction_variance_explained"])

    if random.random() < settings["p_gwas_causal"]:
        # GWAS is causal

        causal_gwas_index = ld_matrix.shape[0] / 2
        gwas_effect_sizes[causal_gwas_index] = potential_gwas_effect

        if random.random() < settings["p_eqtl_causal_given_gwas_causal"]:
            # eQTL is causal

            if random.random() < settings["p_same_causal_given_both_causal"]:
                # Same causal variant for both
                eqtl_effect_sizes[causal_gwas_index] = potential_eqtl_effect
            else:
                # Different eQTL causal variant

                done = False
                for j in range(20):
                    causal_eqtl_index = random.randint(ld_matrix.shape[0] * 1 / 4, ld_matrix.shape[0] * 3 / 4)
                    if ld_matrix[causal_gwas_index][causal_eqtl_index] > settings["max_ld_for_different_causal"] \
                            or causal_gwas_index == causal_eqtl_index:
                        continue
                    else:
                        eqtl_effect_sizes[causal_eqtl_index] = potential_eqtl_effect
                        done = True
                        break

                if not done:
                    # If we can't pick a reasonable eQTL for this after 20 tries, then give up
                    # TODO: If we let this happen to often, it will skew the probabilities,
                    # so find a good way to fix this if it really is happening, or change
                    # the settings to make it easier on the generator.
                    return "Impossible LD matrix"

    else:
        # GWAS is not causal
        if random.random() < settings["p_eqtl_causal_given_no_gwas_causal"]:
            # eQTL is causal

            # Choose causal eQTL; don't let it be right at the edge of the measured region though
            causal_eqtl_index = random.randint(ld_matrix.shape[0] * 1 / 4, ld_matrix.shape[0] * 3 / 4)


    return(gwas_effect_sizes, eqtl_effect_sizes)

if __name__ == "__main__":
    main()
