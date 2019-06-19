from sklearn.ensemble import RandomForestClassifier
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

