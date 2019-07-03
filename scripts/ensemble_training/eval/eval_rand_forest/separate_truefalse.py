import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv

rankless_colocalization_matrix = pd.read_csv("/users/j29tien/colocalization_ML/colocalization_matrix_UNRANKED.tsv", sep="\t")
true_matrix = pd.DataFrame(columns=['coloc_h4', 'rtc_score', 'finemap_clpp', 'finemap_clppmod', 'caviarbf_clpp', 'baseline_neglogpvalue', 'smartbaseline_neglogpvalue', 'smr_neglogpvalue', 'smrheidiadjusted_neglogpvalue', 'gsmr_neglogpvalue'])
false_matrix = pd.DataFrame(columns=['coloc_h4', 'rtc_score', 'finemap_clpp', 'finemap_clppmod', 'caviarbf_clpp', 'baseline_neglogpvalue', 'smartbaseline_neglogpvalue', 'smr_neglogpvalue', 'smrheidiadjusted_neglogpvalue', 'gsmr_neglogpvalue'])

true_matrix = rankless_colocalization_matrix[rankless_colocalization_matrix['colocalization_status'] == True]
true_matrix.reset_index(drop=True, inplace=True)

false_matrix = rankless_colocalization_matrix[rankless_colocalization_matrix['colocalization_status'] == False]
false_matrix.reset_index(drop=True, inplace=True)

true_matrix.to_csv('/users/j29tien/colocalization_ML/true_colocalizations.tsv', sep='\t', index=False)
false_matrix.to_csv('/users/j29tien/colocalization_ML/false_colocalizations.tsv', sep='\t', index=False)

### PLAN: train the random forest classifier on current training dataset, but without the ranks
### PLAN: then, evaluate it on a training dataset that has a higher proportion of TRUE colocalizations

# create set with higher proportion of TRUEs
high_trues = true_matrix.loc[0:int(len(true_matrix.index)*0.667), :]
low_falses = false_matrix.loc[0:int(len(false_matrix.index)*0.1), :]

low_trues = true_matrix.loc[int(len(true_matrix.index)*0.667)+1 :len(true_matrix.index)-1, :]
high_falses = false_matrix.loc[int(len(false_matrix.index)*0.1)+1 :len(false_matrix.index)-1, :]

heavily_true = pd.concat([high_trues,low_falses], ignore_index=True)
heavily_false = pd.concat([low_trues,high_falses], ignore_index=True)

heavily_true.to_csv('/users/j29tien/colocalization_ML/highprop_true.tsv', sep='\t', index=False)
heavily_false.to_csv('/users/j29tien/colocalization_ML/highprop_false.tsv', sep='\t', index=False)
