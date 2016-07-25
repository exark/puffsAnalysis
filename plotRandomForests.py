import os.path as op

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# plotRandomForests plots 2D and 3D scatters for params specified by p2D and p3D respectively. 
# 	nonpuffs, puffs, maybe are lists of param values for all tracks in that class, as returned by runRandomForests.py
#	Scatters are saved as .jpg in savedir. 
#   p2D and p3D are either names of params to plot or indices of params as entered on command line (excluding class variable)

def plotRandomForests(nonpuffs, puffs, maybe, pnames, savedir, p2D = [0,1], save2D = '2Dfig.jpg', p3D = [0,1,2], save3D = '3Dfig.jpg'):

	if len(pnames) < 2:
		raise ValueError('Cannot plot with less than 2 parameters')
		return
	elif ((all(isinstance(x, (int)) for x in p2D)
		or all(isinstance(x, (str)) for x in p2D)
		or all(isinstance(x, (int)) for x in p3D) 
		or all(isinstance(x, (str)) for x in p3D)) is False): 
		raise TypeError('Parameters must be specified either all as indices or names')
		return
	else:
		# Gets indices of parameter names
		if isinstance(p2D[0], str): 
			for i, pname in enumerate(p2D): 
				if not pname in pnames: 
					raise ValueError('Parameter' + '\'pname\'' + 'does not exist for plotting')
				else: 
					p2D[i] = pnames.index(p2D[i])
	
		if isinstance(p3D[0], str):  
			for i,pname in enumerate(p3D): 
				if not pname in pnames: 
					raise ValueError('Parameter' + '\'pname\'' + 'does not exist for plotting')
				else: 
					p3D[i] = pnames.index(p3D[i])

		# Plot 2D scatter
		fig = plt.figure()
		plt.plot(nonpuffs[p2D[0]], nonpuffs[p2D[1]], 'c.', label = 'Nonpuffs')
		plt.plot(maybe[p2D[0]], maybe[p2D[1]], 'g.', label = 'Maybe')
		plt.plot(puffs[p2D[0]], puffs[p2D[1]],'r.', label = 'Puffs')
		plt.xlabel(pnames[p2D[0]])
		plt.ylabel(pnames[p2D[1]])
		plt.title('2D Scatter Plot for RF Results')
		plt.legend(bbox_to_anchor = (1.13, 1.1), numpoints = 1, fontsize = 'small')
		plt.savefig(op.join(savedir,save2D))

		if len(pnames) > 2:
			#Plot 3D scatter
			fig = plt.figure()
			ax = fig.add_subplot(111,projection='3d')

			for c, m, l, xs, ys, zs in [ ('c', 'o', 'Nonpuffs', nonpuffs[p3D[0]], nonpuffs[p3D[1]],nonpuffs[p3D[2]]),
			('g', 'o', 'Maybe', maybe[p3D[0]], maybe[p3D[1]],maybe[p3D[2]]),
			('r','o', 'Puffs', puffs[p3D[0]],puffs[p3D[1]],puffs[p3D[2]]),]:
				ax.scatter(xs,ys,zs,c=c, marker = m, label = l)

			ax.set_xlabel(pnames[p3D[0]], labelpad = 10)
			ax.set_ylabel(pnames[p3D[1]], labelpad = 12)
			ax.set_zlabel(pnames[p3D[2]], labelpad = 6)
			ax.set_title('3D Scatter Plot for RF Results')
			plt.legend(bbox_to_anchor = (0.1, 1.05), scatterpoints = 1, fontsize = 'small')
			plt.savefig(op.join(savedir,save3D))

		plt.show()


