import argparse
import os.path as op
import sys

import numpy as np
import scipy.io as sio
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import r2_score
from sklearn.metrics import precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix

from mat2py import mat2py


def runCrossValidation(train, RFfile):

	train_tracks = []
	for feature in train:
		if feature[0] != 0.:
			train_tracks.append(feature)
	train_tracks = np.array(train_tracks)
	# Gets parameter values for training data
	trainArr = train_tracks[:,1:]
	# Gets class label of all training data
	trainRes = train_tracks[:,0]

	# Convert all NaNs to 0 for RF to work properly
	trainArr = np.nan_to_num(trainArr)
	trainRes = np.nan_to_num(trainRes)

	# Load the classifier
	rf = joblib.load(RFfile)

	# Stratified KFolds cross validation
	cv = StratifiedKFold(n_splits = 5)

	precision   = []
	accuracy    = []
	sensitivity = []
	matthews    = []
	r2          = []
	f1          = []
	auroc       = []
	cm          = [[0, 0], [0, 0]]

	for train_index, test_index in cv.split(trainArr, trainRes):
	    probas_     = rf.fit(trainArr[train_index], trainRes[train_index]).predict_proba(trainArr[test_index])
	    classes     = rf.fit(trainArr[train_index], trainRes[train_index]).predict(trainArr[test_index])
	   # r2          = np.append(r2, (r2_score(trainRes[test_index], probas_[:, 1])))
	    precision   = np.append(precision, (precision_score(trainRes[test_index], classes)))
	   # auroc       = np.append(auroc, (roc_auc_score(trainRes[test_index], classes)))
	    accuracy    = np.append(accuracy, (accuracy_score(trainRes[test_index], classes)))
	    sensitivity = np.append(sensitivity, (recall_score(trainRes[test_index], classes)))
	    f1          = np.append(f1, (f1_score(trainRes[test_index], classes)))
	   # matthews    = np.append(matthews, (matthews_corrcoef(trainRes[test_index], classes)))
	    #cma         = np.add(cma, (confusion_matrix(trainRes[test_index], classes)))

	# cma         = np.array(cma)
	# r2          = np.array(r2)
	precision   = np.array(precision)
	accuracy    = np.array(accuracy)
	sensitivity = np.array(sensitivity)
	f1          = np.array(f1)
	# auroc       = np.array(auroc)
	# matthews    = np.array(matthews)

	return accuracy, precision, sensitivity, f1

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Random Forest classification of tracks')
	parser.add_argument('RFfile', help="Classifier to load or name of file to save classifier in")
	parser.add_argument('training', help="Labeled data to be used for crossvalidation in npy format")
	args = parser.parse_args()

	savedir = op.dirname(args.RFfile)

	train = np.load(args.training)
	fields = list(train.dtype.names)

	train = np.array(train.tolist())

	accuracy, precision, sensitivity, f1 = runCrossValidation(train, args.RFfile)

	#n = open(op.join(savedir, 'notes.txt'), 'a')
	print('\n KF Accuracy: %0.2f (+/- %0.2f)' % (accuracy.mean(), accuracy.std() * 2))
	print('\n KF Precision: %0.2f (+/- %0.2f)' % (precision.mean(), precision.std() * 2))
	print('\n KF Sensitivity: %0.2f (+/- %0.2f)' % (sensitivity.mean(), sensitivity.std() * 2))
	print('\n KF F1: %0.2f (+/- %0.2f)' % (f1.mean(), f1.std() * 2))
	#n.close()
