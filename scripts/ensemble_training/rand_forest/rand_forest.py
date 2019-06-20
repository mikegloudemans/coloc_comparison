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

# initiate and fit the random forest classifier -- need to tweak parameters --> n_estimators = 150, max_depth = 7, min_samples_leaf=1
randfor = RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini', max_depth=7, max_features='auto', max_leaf_nodes=None, min_impurity_decrease=0.0, min_impurity_split=None, min_samples_leaf=1, min_samples_split=2, min_weight_fraction_leaf=0.0, n_estimators=150, n_jobs=None, oob_score=False, random_state=None, verbose=0, warm_start=False)
randfor.fit(trainX, trainY)

#predict using cross-validation dataset
#load validation data and KEY
validation_save_path = "/users/j29tien/colocalization_ML/validation_set.tsv"
validationKEY_save_path = "/users/j29tien/colocalization_ML/validation_set_KEY.tsv"
df = pd.read_csv(validation_save_path, sep="\t")   # read .tsv file into memory
validX = df.to_numpy()  # access the numpy array containing values
df = pd.read_csv(validationKEY_save_path, sep="\t")   # read .tsv file into memory
validY = df.to_numpy()  # access the numpy array containing values


#predict and score using validation data
predY = randfor.predict(validX)
print(predY)
print("Accuracy score for model:")
print(accuracy_score(validY, predY))


predProbY = randfor.predict_proba(validX)
print(predProbY)
# save the prediction probabilities into a file
with open('/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/rand_forest/prediction_probabilities.csv', 'w') as probability_file:
    writer = csv.writer(probability_file, delimiter = ',', quotechar   = '"', quoting = csv.QUOTE_NONE, escapechar = ' ')
    for row in predProbY:
        writer.writerow(row)
#feature importances for the random forest
#print(randfor.feature_importances_)
importances = randfor.feature_importances_
std = np.std([tree.feature_importances_ for tree in randfor.estimators_], axis=0)
indices = np.argsort(importances)[::-1]
# Print the feature ranking
print("Feature ranking:")
with open('/users/j29tien/colocalization_ML/coloc_comparison/scripts/ensemble_training/rand_forest/feature_importances.csv', 'w') as feature_file:
	for f in range(validX.shape[1]):
	    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
	    feature_file.write("%d,\"feature %d\",%f" % (f + 1, indices[f], importances[indices[f]]))
	    feature_file.write('\n')

