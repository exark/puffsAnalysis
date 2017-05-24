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


def runCrossValidation(train, test, RFfile):

	train_tracks = []
	for feature in train:
		if feature[0] != 0.:
			train_tracks.append(feature)
	train_tracks = np.array(train_tracks)
	test_tracks = np.array(test)
	# Gets parameter values for training data
	trainArr = train_tracks[:,1:]

	# Gets class label of all training data
	trainRes = train_tracks[:,0]

	# Gets parameter values for test data
	testArr = test_tracks[:,1:]

	# Convert all NaNs to 0 for RF to work properly
	trainArr = np.nan_to_num(trainArr)
	trainRes = np.nan_to_num(trainRes)
	testArr = np.nan_to_num(testArr)

	if not op.exists(RFfile):
		# Train the classifier 
		rf = RandomForestClassifier(n_estimators = 500, oob_score = True, class_weight="balanced")
		rf.fit(trainArr, trainRes)

		# Save the classifier as RFfile
		joblib.dump(rf, RFfile, compress=True)

	else:
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
	   # sensitivity = np.append(sensitivity, (recall_score(trainRes[test_index], classes)))
	   # f1          = np.append(f1, (f1_score(trainRes[test_index], classes)))
	   # matthews    = np.append(matthews, (matthews_corrcoef(trainRes[test_index], classes)))
	    #cma         = np.add(cma, (confusion_matrix(trainRes[test_index], classes)))

	# cma         = np.array(cma)
	# r2          = np.array(r2)
	precision   = np.array(precision)
	accuracy    = np.array(accuracy)
	# sensitivity = np.array(sensitivity)
	# f1          = np.array(f1)
	# auroc       = np.array(auroc)
	# matthews    = np.array(matthews)

	return accuracy, precision 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Random Forest classification of tracks')
	parser.add_argument('RFfile', help="Classifier to load or name of file to save classifier in")
	parser.add_argument('training', help="Numpy array or HDF5-encoded MATLAB file")
	parser.add_argument('--fields', default = [], nargs="+", help="Field(s) to extract from .mat file")
	parser.add_argument('--testing', dest='testing',default=[])
	args = parser.parse_args()

	traindir = op.join(op.dirname(op.dirname(args.training)),'Classification')
	savedir = traindir

	if (args.training).endswith('.npy'):
		train = np.load(args.training)
		fields = list(train.dtype.names)
	else:
		if not args.fields or len(args.fields) <= 1:
			raise ValueError('Must import at least two parameters (including labeled variable) for RF classification')
			quit()
		else:
			fields = args.fields
		train = mat2py(args.training, fields, traindir)

	if args.testing:
		testdir = op.join(op.dirname(op.dirname(args.testing)),'Classification')
		if (args.testing).endswith('.npy'):
			test = np.load(args.testing)
		else:
			test = mat2py(args.testing, fields, testdir)
		savedir = testdir
	else:
		test = np.array(train)

	train = np.array(train.tolist())
	test = np.array(test.tolist())
	
	accuracy, precision = runCrossValidation(train, test, args.RFfile)

	#Get feature importances from classifier
	importances = rf.feature_importances_
	std = np.std([tree.feature_importances_ for tree in rf.estimators_], axis=0)
	indices = np.argsort(importances)[::-1]

	n = open(op.join(savedir, 'notes.txt'), 'a')
	# n.write('\n Classifier built from: ' + args.training)
	# n.write('\n Params used: ' 	  + ','.join(fields[1:]))
	# n.write('\n Puffs/Total: ' 	  + str(ntracks[2]) + '/' + str(ntracks[0]) + ' (' + str((ntracks[2]/ntracks[0]) *100) + '%)')
	# n.write('\n Nonpuffs/Total: ' + str(ntracks[1]) + '/' + str(ntracks[0]) + ' (' + str((ntracks[1]/ntracks[0]) *100) + '%)')
	# n.write('\n Feature Importances: ')
	# for f in range(len(fields)-1):
	# 	n.write("\n\t %d. %s (%f)" % (f + 1, fields[indices[f]+1], importances[indices[f]]))
	n.write('\n OOB Error: ' 	  + str(rf.oob_score_))
	n.write('\n KF Accuracy: %0.2f (+/- %0.2f)' % (accuracy.mean(), accuracy.std() * 2))
	n.write('\n KF Precision: %0.2f (+/- %0.2f)' % (precision.mean(), precision.std() * 2))
	n.close()