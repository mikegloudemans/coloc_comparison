import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

rankless_colocalization_matrix = pd.read_csv("/users/j29tien/colocalization_ML/______.tsv", sep="\t")
true_matrix = pd.DataFrame(columns=['coloc_h4', 'rtc_score', 'finemap_clpp', 'finemap_clppmod', 'caviarbf_clpp', 'baseline_neglogpvalue', 'smartbaseline_neglogpvalue', 'smr_neglogpvalue', 'smrheidiadjusted_neglogpvalue', 'gsmr_neglogpvalue'])
false_matrix = pd.DataFrame(columns=['coloc_h4', 'rtc_score', 'finemap_clpp', 'finemap_clppmod', 'caviarbf_clpp', 'baseline_neglogpvalue', 'smartbaseline_neglogpvalue', 'smr_neglogpvalue', 'smrheidiadjusted_neglogpvalue', 'gsmr_neglogpvalue'])

true_matrix = colocalization_matrix[colocalization_matrix['colocalization_status'] == True]
false_matrix = colocalization_matrix[colocalization_matrix['colocalization_status'] == False]

true_matrix.to_csv('/users/j29tien/colocalization_ML/true_colocalizations.tsv', sep='\t', index=False)
false_matrix.to_csv('/users/j29tien/colocalization_ML/false_colocalizations.tsv', sep='\t', index=False)

### PLAN: train the random forest classifier on current training dataset, but without the ranks
### PLAN: then, evaluate it on a training dataset that has a higher proportion of TRUE colocalizations

# create set with higher proportion of TRUEs


