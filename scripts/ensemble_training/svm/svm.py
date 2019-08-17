from sklearn import svm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
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

# initiate and fit the support vector machine classifier -- ("one-vs-rest" multi-class strategy)
suppvec = svm.SVC(gamma='scale', decision_function_shape='ovo')
# suppvec = svm.LinearSVC() --> DOES NOT CONVERGE, but allows us to use coef_ to view weights
suppvec.fit(trainX, trainY)

#predict using cross-validation dataset
#load validation data and KEY
validation_save_path = "/users/j29tien/colocalization_ML/validation_set.tsv"
validationKEY_save_path = "/users/j29tien/colocalization_ML/validation_set_KEY.tsv"
df = pd.read_csv(validation_save_path, sep="\t")   # read .tsv file into memory
validX = df.to_numpy()  # access the numpy array containing values
df = pd.read_csv(validationKEY_save_path, sep="\t")   # read .tsv file into memory
validY = df.to_numpy()  # access the numpy array containing values


#predict and score using validation data
predY = suppvec.predict(validX)
print(predY)
print("Accuracy score for model:")
print(accuracy_score(validY, predY))

# save support vector properties and weights
properties = suppvec.support_vectors_
#weights = suppvec.coef_ ONLY viewable with LinearSVC
with open('/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/svm/support_vectors.csv', 'w') as sv_file:
    writer = csv.writer(sv_file, delimiter = ',', quotechar   = '"', quoting = csv.QUOTE_NONE, escapechar = ' ')
    for row in properties:
        writer.writerow(row)

predProbY = suppvec.decision_function(validX) #predict_proba and predict_log_proba are only enabled when probability=True in constructor. not recommended; expensive operation. Returns confidence scores per (sample, class) combination. In the binary case, confidence score for self.classes_[1] where >0 means this class would be predicted.
print(predProbY)
print(predProbY.shape)
# save the prediction probabilities into a file
with open('/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/svm/prediction_probabilities.csv', 'w') as probability_file:
    writer = csv.writer(probability_file, delimiter = ',', quotechar   = '"', quoting = csv.QUOTE_NONE, escapechar = ' ')
    writer.writerow(predProbY) # for SVM with two classes (binary), returns a single score for each sample


