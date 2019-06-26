import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from itertools import cycle

from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from scipy import interp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score


training_save_path = "/users/j29tien/colocalization_ML/training_set.tsv"
trainingKEY_save_path = "/users/j29tien/colocalization_ML/training_set_KEY.tsv"

#load trainX, the training dataset
df = pd.read_csv(training_save_path, sep="\t")   # read .tsv file into memory
trainX = df.to_numpy()  # access the numpy array containing values

#load trainY, the answer key for training
df = pd.read_csv(trainingKEY_save_path, sep="\t")   # read .tsv file into memory
trainY = df.to_numpy()  # access the numpy array containing values
trainY = np.ravel(trainY)

# Binarize the output
trainY = label_binarize(trainY, classes=[0, 1])
n_classes = trainY.shape[1]

#create and fit multiclass logistic regression model
logreg = LogisticRegression(multi_class='ovr') #specify 'one versus rest' option 
logreg = logreg.fit(trainX, trainY)

#predict using test dataset
#load test data and KEY
test_save_path = "/users/j29tien/colocalization_ML/test_set.tsv"
testKEY_save_path = "/users/j29tien/colocalization_ML/test_set_KEY.tsv"
df = pd.read_csv(test_save_path, sep="\t")   # read .tsv file into memory
testX = df.to_numpy()  # access the numpy array containing values
df = pd.read_csv(testKEY_save_path, sep="\t")   # read .tsv file into memory
testY = df.to_numpy()  # access the numpy array containing values

predY = logreg.predict(testX)
print(predY)
print(accuracy_score(testY, predY))
predProbY = logreg.predict_proba(testX)
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
    fpr[i], tpr[i], _ = roc_curve(testY, testX[:, i*2]) #pulling out rows of raw scores from each method (row 0, 2, 4, ... 18)
    roc_auc[i] = auc(fpr[i], tpr[i])
# store ensemble performance at n_submethods'th index
fpr[n_submethods],tpr[n_submethods],_ = roc_curve(testY, predProbY[:,1]) # index gives probability for 0/FALSE, then 1/TRUE
roc_auc[n_submethods] = auc(fpr[n_submethods], tpr[n_submethods])

#weights = np.ravel(logreg.coef_) # returns a matrix of weights (coefficients)
#print(weights)

# plot ROC for the multilogreg ensemble method
plt.figure(figsize=(20, 15))
lw = 2
plt.plot(fpr[n_submethods], tpr[n_submethods], color='darkorange', linewidth=4, linestyle = ':', label='ensemble, logisitic regression (area = %0.2f)' % roc_auc[n_submethods])
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
plt.savefig("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/eval/comp_ROC.png")

