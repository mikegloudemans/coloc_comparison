import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn import datasets
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

#create and fit multiclass logistic regression model
logreg = LogisticRegression(multi_class='ovr') #specify 'one versus rest' option 
logreg.fit(trainX, trainY)

#predict using cross-validation dataset
#load validation data and KEY
validation_save_path = "/users/j29tien/colocalization_ML/validation_set.tsv"
validationKEY_save_path = "/users/j29tien/colocalization_ML/validation_set_KEY.tsv"
df = pd.read_csv(validation_save_path, sep="\t")   # read .tsv file into memory
validX = df.to_numpy()  # access the numpy array containing values
df = pd.read_csv(validationKEY_save_path, sep="\t")   # read .tsv file into memory
validY = df.to_numpy()  # access the numpy array containing values

predY = logreg.predict(validX)
print(predY)
print(accuracy_score(validY, predY))
# worth saving the predictions in a file?

weights = logreg.coef_ # returns a matrix of weights (coefficients)
print(weights)
#save weights in a file

