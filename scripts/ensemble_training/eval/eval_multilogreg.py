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

# Compute ROC curve and ROC area for each class
fpr = dict()
tpr = dict()
roc_auc = dict()
#for i in range(n_classes):
    #fpr[i], tpr[i], _ = roc_curve(testY[:, i], predProbY[:, i])
    #roc_auc[i] = auc(fpr[i], tpr[i])

fpr,tpr,_ = roc_curve(testY, predProbY[:,1])
roc_auc = auc(fpr, tpr)

# Compute micro-average ROC curve and ROC area --> NOT NECESSARY, not multiclass (only TRUE or FALSE)
#fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
#roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])


weights = np.ravel(logreg.coef_) # returns a matrix of weights (coefficients)
print(weights)

# plot ROC for the multilogreg ensemble method
plt.figure()
lw = 2
plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic example')
plt.legend(loc="lower right")
plt.savefig("/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/eval/multilogreg_ROC.png")

