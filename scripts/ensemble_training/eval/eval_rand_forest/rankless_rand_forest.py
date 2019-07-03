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


training_save_path = "/users/j29tien/colocalization_ML/highprop_false.tsv" # train on dataset that has a larger proportion of falses. (this dataset happens to be larger; more suitable for training)

#load trainX, the training dataset
df = pd.read_csv(training_save_path, sep="\t")   # read .tsv file into memory
trainX = df.drop('colocalization_status', 1).to_numpy()  # access the numpy array containing values

#load trainY, the answer key for training
trainY = df['colocalization_status'].to_numpy()  # access the numpy array containing values
trainY = np.ravel(trainY)

# Binarize the output
trainY = label_binarize(trainY, classes=[0, 1])
n_classes = trainY.shape[1]

# initiate and fit the random forest classifier -- need to tweak parameters --> n_estimators = 150, max_depth = 7, min_samples_leaf=1
randfor = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini', max_depth=7, max_features='auto', max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, min_samples_leaf=1, min_samples_split=2, min_weight_fraction_leaf=0.0, n_estimators=150, n_jobs=None, oob_score=False, random_state=None, verbose=0, warm_start=False)
randfor.fit(trainX, trainY)


#predict using test dataset
#load test data and KEY
test_save_path = "/users/j29tien/colocalization_ML/highprop_true.tsv"
df = pd.read_csv(test_save_path, sep="\t")   # read .tsv file into memory
testX = df.drop('colocalization_status', 1).to_numpy()  # access the numpy array containing values

testY = df['colocalization_status'].to_numpy()  # access the numpy array containing values

#predict and score using test data
predY = randfor.predict(testX)
print(predY)
print("Accuracy score for model:")
print(accuracy_score(testY, predY))
predProbY = randfor.predict_proba(testX)
print(predProbY)

#load submethod names
method_names_path = "/users/j29tien/colocalization_ML/feature_names.tsv"
df = pd.read_csv(method_names_path, sep="\t")
method_names = df.columns.values

# Compute ROC curve and ROC area for each method
fpr = dict()
tpr = dict()
roc_auc = dict()

n_submethods = 10 #coloc, rtc, finemap-clpp, finemap-clpp_mod, caviarbf, baseline, smart_baseline, smr, smr-heidi, gsmr
for i in range(n_submethods):
    fpr[i], tpr[i], _ = roc_curve(testY, testX[:, i]) #pulling out rows of raw scores from each method (row 0, 1, 2, ... 9)
    roc_auc[i] = auc(fpr[i], tpr[i])
# store ensemble performance at n_submethods'th index
fpr[n_submethods],tpr[n_submethods],_ = roc_curve(testY, predProbY[:,1]) # index gives probability for 0/FALSE, then 1/TRUE
roc_auc[n_submethods] = auc(fpr[n_submethods], tpr[n_submethods])

# plot ROC for the random forest ensemble method
plt.figure(figsize=(20, 15))
lw = 2
plt.plot(fpr[n_submethods], tpr[n_submethods], color='darkorange', linewidth=4, linestyle = ':', label='ensemble, random forest (area = %0.2f)' % roc_auc[n_submethods])
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')

colors = cycle(['aqua', 'dimgray', 'cornflowerblue', 'orangered', 'gold', 'chartreuse', 'forestgreen', 'darkviolet', 'deeppink', 'crimson'])
for i, color in zip(range(n_submethods), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=lw,
             label='{0} (area = {1:0.2f})'
             ''.format(method_names[i*2], roc_auc[i]))

plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Comparison of colocalization identification methods')
plt.legend(loc="lower right")
plt.savefig("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/eval/eval_rand_forest/rankless_comp_ROC.png")


#feature importances for the random forest
#print(randfor.feature_importances_)
importances = randfor.feature_importances_
std = np.std([tree.feature_importances_ for tree in randfor.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
# Print the feature ranking
print("Feature ranking:")
for f in range(testX.shape[1]):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

