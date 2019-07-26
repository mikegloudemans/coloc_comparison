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

matrix_save_path = "/users/j29tien/colocalization_ML/colocalization_matrix_NEW.tsv"
matrix = pd.read_csv(matrix_save_path, sep="\t")
trainX, testX, trainY, testY = train_test_split(matrix.drop('colocalization_status', 1), matrix['colocalization_status'], test_size=0.25)

#convert training and testing data to numpy arrays
trainX = trainX.to_numpy()  # access the numpy array containing values
trainY = trainY.to_numpy()  # access the numpy array containing values
testX = testX.to_numpy()  # access the numpy array containing values
testY = testY.to_numpy()  # access the numpy array containing values

# Binarize the output for training
trainY = label_binarize(trainY, classes=[0, 1])
n_classes = trainY.shape[1]


# initiate and fit the random forest classifier -- need to tweak parameters --> n_estimators = 150, max_depth = 7, min_samples_leaf=1
randfor = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini', max_depth=7, max_features='auto', max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, min_samples_leaf=1, min_samples_split=2, min_weight_fraction_leaf=0.0, n_estimators=150, n_jobs=None, oob_score=False, random_state=None, verbose=0, warm_start=False)
randfor.fit(trainX, trainY)


#predict and score using test data
predY = randfor.predict(testX)
print(predY)
print("Accuracy score for model:")
print(accuracy_score(testY, predY))
predProbY = randfor.predict_proba(testX)
print(predProbY)

#load submethod names
method_names = matrix.columns.values

# Compute ROC curve and ROC area for each method
fpr = dict()
tpr = dict()
roc_auc = dict()

n_submethods = int(len(matrix.columns)/2) #coloc, rtc, finemap-clpp, finemap-clpp_mod, caviarbf, baseline, smart_baseline, smr, smr-heidi, gsmr, twas
for i in range(n_submethods):
    fpr[i], tpr[i], _ = roc_curve(testY, testX[:, i*2]) #pulling out rows of raw scores from each method (row 0, 2, 4, ... 18)
    roc_auc[i] = auc(fpr[i], tpr[i])
# store ensemble performance at n_submethods'th index
fpr[n_submethods],tpr[n_submethods],_ = roc_curve(testY, predProbY[:,1]) # index gives probability for 0/FALSE, then 1/TRUE
roc_auc[n_submethods] = auc(fpr[n_submethods], tpr[n_submethods])

# plot ROC for the random forest ensemble method
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=28)    # fontsize of the tick labels
plt.rc('ytick', labelsize=28)    # fontsize of the tick labels
plt.rc('legend', fontsize=22)    # legend fontsize
plt.rc('axes', titlesize=28)

plt.figure(figsize=(20, 15))
lw = 6
colors = cycle(['aqua', 'dimgray', 'cornflowerblue', 'orangered', 'gold', 'chartreuse', 'forestgreen', 'darkviolet', 'deeppink', 'crimson', 'rosybrown']) #keeping same number of colors as n_submethods
for i, color in zip(range(n_submethods), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='{0} (area = {1:0.2f})'
             ''.format(method_names[i*2], roc_auc[i]))

plt.plot(fpr[n_submethods], tpr[n_submethods], color='darkorange', linewidth=15, linestyle = ':', label='ENSEMBLE, random forest (area = %0.2f)' % roc_auc[n_submethods])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')


plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
#plt.xlabel('False Positive Rate')
#plt.ylabel('True Positive Rate')
plt.title('Comparison of colocalization identification methods')
plt.legend(loc="lower right")
plt.savefig("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/eval/eval_rand_forest/comp_ROC_NEW.png")





#feature importances for the random forest
#print(randfor.feature_importances_)
importances = randfor.feature_importances_
std = np.std([tree.feature_importances_ for tree in randfor.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
# Print the feature ranking
print("Feature ranking:")
for f in range(testX.shape[1]):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

