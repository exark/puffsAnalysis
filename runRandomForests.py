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

# runRandomForests builds a classifier based on train and runs it on test (generated from mat2py.py)
#	Classifier is saved as RFfile. Classifier results are saved as RFresults.mat in savedir.
# 	Parameter values for tracks in each class label are returned.
# runRandomForests assumes first parameter passed is always the class variable.

def runRandomForests(train, test, RFfile, savedir):

	# Add the track to training data if it has been scored
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

	# Use classifier to predict class of all test data
	testRes = rf.predict_proba(testArr)

	# Create 2D array for all class labels where for each one,
	# class[1] is a list of track indices predicted to be that class
	# class[2],[3],[4] (...) are that track's 1st, 2nd, 3rd (...) parameter values
	nparams = np.zeros((testArr.shape[1]))
	puffs = [[] for _ in range(len(nparams)+1)]
	nonpuffs = [[] for _ in range(len(nparams)+1)]

	for i, res in enumerate(testRes):
		if res[0] < 0.314917127072:
			nonpuffs[0].append(i)
			for j, param in enumerate(nparams):
				nonpuffs[j+1].append(testArr[i,j])
		else:
			puffs[0].append(i)
			for j, param in enumerate(nparams):
				puffs[j+1].append(testArr[i,j])

	# Create a dictionary of track indices for all class labels (adjust to 1 indexing for MATLAB)
	nonp= [x+1 for x in nonpuffs[0]]
	p= [x+1 for x in puffs[0]]

	labels = ['nonpuffs', 'puffs']
	ID = [nonp,p]
	idx = dict(zip(labels, ID))

	# Save dictionary as RF.mat
	sio.savemat(op.join(savedir,'RFresults'), idx)
	ntracks = [len(nonpuffs[0]) + len(puffs[0])]
	for x in [len(nonpuffs[0]), len(puffs[0])]:
		ntracks.append(x)

	return nonpuffs[1:], puffs[1:], ntracks, rf,

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Random Forest classification of tracks')
	parser.add_argument('RFfile', help="Classifier to load or name of file to save classifier in")
	parser.add_argument('training', help="Numpy array or HDF5-encoded MATLAB file")
	parser.add_argument('--fields', default = [], nargs="+", help="Field(s) to extract from .mat file")
	parser.add_argument('--testing', dest='testing',default=[])
	args = parser.parse_args()

	testdir = op.join(op.dirname(op.dirname(args.testing)),'Classification')
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
		if (args.testing).endswith('.npy'):
			test = np.load(args.testing)
		else:
			test = mat2py(args.testing, fields, testdir)
		savedir = testdir
	else:
		test = np.array(train)

	train = np.array(train.tolist())
	test = np.array(test.tolist())

	nonpuffs, puffs, ntracks, rf = runRandomForests(train, test, args.RFfile, savedir)

	#Get feature importances from classifier
	importances = rf.feature_importances_
	std = np.std([tree.feature_importances_ for tree in rf.estimators_], axis=0)
	indices = np.argsort(importances)[::-1]

	n = open(op.join(savedir, 'notes.txt'), 'a')
	n.write('\n Classifier built from: ' + args.training)
	n.write('\n Params used: ' 	  + ','.join(fields[1:]))
	n.write('\n Puffs/Total: ' 	  + str(ntracks[2]) + '/' + str(ntracks[0]) + ' (' + str((ntracks[2]/ntracks[0]) *100) + '%)')
	n.write('\n Nonpuffs/Total: ' + str(ntracks[1]) + '/' + str(ntracks[0]) + ' (' + str((ntracks[1]/ntracks[0]) *100) + '%)')
	n.write('\n OOB Error: ' 	  + str(rf.oob_score_))
	n.write('\n Feature Importances: ')
	for f in range(len(fields)-1):
		n.write("\n\t %d. %s (%f)" % (f + 1, fields[indices[f]+1], importances[indices[f]]))
	n.close()
