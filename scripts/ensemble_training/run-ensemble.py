## This script runs the ensemble, given the submethods that are to be included.
## Command line format: python eval_rand_forest.py METHOD_TO_REMOVE1 METHOD_TO_REMOVE2

import sys
## Runs on Python version >= 3.6.8, sklearn version = 0.21.2

from sklearn.ensemble import RandomForestClassifier

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from itertools import cycle

from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from scipy import interp
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

###################
## PREPROCESSING ##
###################

# Load matrix of simulated data. Default is rankless, because real data will not have ranks
matrix_save_path = "/users/j29tien/colocalization_ML/colocalization_matrix_NEW_UNRANKED.tsv"
#matrix_save_path = "/users/j29tien/colocalization_ML/colocalization_matrix_NEW.tsv"
matrix = pd.read_csv(matrix_save_path, sep="\t")

# remove any submethods based on input
print(sys.argv)
if len(sys.argv) > 1:
    matrix.drop(sys.argv[1:len(sys.argv)], axis=1, inplace=True)
print(matrix)
# split matrix into test and train
trainX, testX, trainY, testY = train_test_split(matrix.drop('colocalization_status', 1), matrix['colocalization_status'], test_size=0.25)

#convert training and testing data to numpy arrays
trainX = trainX.to_numpy()  # access the numpy array containing values
trainY = trainY.to_numpy()  # access the numpy array containing values

## WARNING: will not be using testX and testY
testX = testX.to_numpy()  # access the numpy array containing values
testY = testY.to_numpy()  # access the numpy array containing values

# Binarize the output for training
#trainY = label_binarize(trainY, classes=[0, 1])
#n_classes = trainY.shape[1]

# Load REAL data to be predicted on
real_data_path = "/users/j29tien/colocalization_ML/real_matrix.tsv"
realX = pd.read_csv(real_data_path, sep='\t').to_numpy()


####################
## TRAIN ENSEMBLE ##
####################

# initiate and fit the random forest classifier -- need to tweak parameters --> n_estimators = 150, max_depth = 7, min_samples_leaf=1
#Python 3 version
randfor = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini', max_depth=7, max_features='auto', max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, min_samples_leaf=1, min_samples_split=2, min_weight_fraction_leaf=0.0, n_estimators=150, n_jobs=None, oob_score=False, random_state=None, verbose=0, warm_start=False)

#Python 2 version
#randfor = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini', max_depth=7, max_features='auto', max_leaf_nodes=None, min_impurity_split=0, min_samples_leaf=1, min_samples_split=2, min_weight_fraction_leaf=0.0, n_estimators=150, oob_score=False, random_state=None, verbose=0, warm_start=False)

randfor.fit(trainX, trainY)

###############################
## COLOCALIZATION PREDICTION ##
###############################

#predict using real data
predY = randfor.predict(realX)
print(predY)

# This is only for evaluating model
#print("Accuracy score for model:")
#print(accuracy_score(testY, predY))

predProbY = randfor.predict_proba(realX)
print("COLOCALIZATION SCORES:")
print(predProbY[:,1])
pd.DataFrame(predProbY[:,1]).to_csv("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/realdata_scores.tsv", sep='\t', header=False, index=False)


############################
## PERFORMANCE EVALUATION ##
############################

##load submethod names
#method_names = np.delete(matrix.columns.values, len(matrix.columns.values)-1) #delete the last column header ('colocalization_status')
#
## Compute ROC curve and ROC area for each method
#fpr = dict()
#tpr = dict()
#roc_auc = dict()
#
#n_submethods = int(len(method_names)) #coloc, rtc, finemap-clpp, finemap-clpp_mod, caviarbf, baseline, smart_baseline, smr, smr-heidi, gsmr, twas
#for i in range(n_submethods):
#    fpr[i], tpr[i], _ = roc_curve(testY, testX[:, i]) #pulling out rows of raw scores from each method (row 1, 2, 3, ... 11)
#    roc_auc[i] = auc(fpr[i], tpr[i])
## store ensemble performance at n_submethods'th index
#fpr[n_submethods],tpr[n_submethods],_ = roc_curve(testY, predProbY[:,1]) # index gives probability for 0/FALSE, then 1/TRUE
#roc_auc[n_submethods] = auc(fpr[n_submethods], tpr[n_submethods])
#
## plot ROC for the random forest ensemble method
#plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=28)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=28)    # fontsize of the tick labels
#plt.rc('legend', fontsize=22)    # legend fontsize
#plt.rc('axes', titlesize=28)
#
#plt.figure(figsize=(20, 15))
#lw = 6
#colors = cycle(['aqua', 'dimgray', 'cornflowerblue', 'orangered', 'gold', 'chartreuse', 'forestgreen', 'darkviolet', 'deeppink', 'crimson', 'rosybrown']) #keeping same number of colors as n_submethods
#for i, color in zip(range(n_submethods), colors):
#    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
#             label='{0} (area = {1:0.2f})'
#             ''.format(method_names[i], roc_auc[i]))
#
#plt.plot(fpr[n_submethods], tpr[n_submethods], color='darkorange', linewidth=15, linestyle = ':', label='ENSEMBLE, random forest (area = %0.2f)' % roc_auc[n_submethods])
#plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
#
#
#plt.xlim([0.0, 1.0])
#plt.ylim([0.0, 1.05])
##plt.xlabel('False Positive Rate')
##plt.ylabel('True Positive Rate')
#plt.title('Comparison of colocalization identification methods')
#plt.legend(loc="lower right")
#plt.savefig("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/test_ROC.png")
#
#
##feature importances for the random forest
##print(randfor.feature_importances_)
#importances = randfor.feature_importances_
#std = np.std([tree.feature_importances_ for tree in randfor.estimators_], axis=0)
#indices = np.argsort(importances)[::-1]
## Print the feature ranking
#print("Feature ranking:")
#for f in range(testX.shape[1]):
#    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

