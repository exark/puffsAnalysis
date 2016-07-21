import argparse
import os.path as op
import numpy as np

from mat2py import mat2py
from runRandomForests import runRandomForests
from plotRandomForests import plotRandomForests

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

	nonpuffs, puffs, maybe, ntracks = runRandomForests(train, test, args.RFfile, savedir)

	n = open(op.join(savedir, 'notes.txt'), 'a')
	n.write('\n Classifier built from: ' + args.training)
	n.write('\n Params used: ' + ', '.join(fields[1:]))
	n.write('\n Puffs/Total: ' + str(ntracks[2]) + '/' + str(ntracks[0]) + ' (' + str((ntracks[2]/ntracks[0]) *100) + '%)')
	n.write('\n Nonpuffs/Total: ' + str(ntracks[1]) + '/' + str(ntracks[0]) + ' (' + str((ntracks[1]/ntracks[0]) *100) + '%)')
	n.write('\n Maybe/Total: ' + str(ntracks[3]) + '/' + str(ntracks[0]) + ' (' + str((ntracks[3]/ntracks[0]) *100) + '%)\n')
	n.close()

	plotRandomForests(nonpuffs, puffs, maybe, fields[1:], savedir, p2D=[2,0])

	

