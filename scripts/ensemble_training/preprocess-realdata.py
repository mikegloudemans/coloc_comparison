## PLEASE READ:
# This script will not run on its own, it is merely a collection of key commands that can be run in Python shell
##

import sys

import numpy as np
import pandas as pd
import csv
from itertools import cycle

in_matrix_path = "/users/mgloud/projects/coloc_comparisons/jeremy/colocalization_matrix_real_data.tsv"
in_matrix = pd.read_csv(in_matrix_path, sep='\t')

# remove the columns with metadata
in_matrix.drop(in_matrix.columns[[0, 1, 2, 3, 4]], axis=1, inplace=True)

# remove test sites (rows) with NA -- sites where method wasn't able to run
in_matrix.dropna(inplace=True)

# rename the column names to vv
#in_matrix = pd.read_csv("/users/j29tien/colocalization_ML/real_matrix.tsv", sep='\t')

in_matrix = in_matrix[['COLOC', 'RTC', 'FINEMAP_STANDARD', 'FINEMAP_MODIFIED', 'BASELINE', 'SMARTBASELINE', 'SMR', 'GSMR', 'TWAS']]
# COLOC, RTC, FINEMAP_STANDARD, FINEMAP_MODIFIED, CAVIARBF, BASELINE, SMARTBASELINE, SMR, SMR_HEIDIADJUSTED, GSMR, TWAS

in_matrix.to_csv("/users/j29tien/colocalization_ML/real_matrix.tsv", sep='\t', index=False)


