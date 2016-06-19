import math
import numpy as np
import copy
import operator
import os
import argparse
import psutil
import time
import pdb

import matplotlib.pyplot as plt

try:
	import libcarma as libcarma
except ImportError:
	print 'libcarma is not setup. Setup libcarma by sourcing bin/setup.sh'
	sys.exit(1)
import util.mcmcviz as mcmcviz
from util.mpl_settings import set_plot_params
import util.triangle as triangle

def _readCRTS(filename, path):
	filePath = os.path.join(path, filename)
	with open(filePath) as fileIn:
		allLines = fileIn.readlines()
	num = len(allLines) - 1
	tList = list()
	magList = list()
	magerrList = list()
	yList = list()
	yerrList = list()
	for i in xrange(1, num+1):
		words = allLines[i].rstrip('\n').split(',')
		f, ferr = libcarma.pogsonFlux(float(words[1]), float(words[2]))
		tList.append(float(words[5]))
		yList.append(f)
		#yerrList.append(ferr)
		yerrList.append(ferr/5.2496267109062256) # CRTS errors may be overestimated. This scaling factor was empiracally detrmined by looking at the SF of the light curve and making np.mean(PG1302_102.yerr)/math.sqrt(sfE[10]) == 1.0
		#yerrList.append(ferr/10.931733864782307) # CRTS errors may be overestimated. This scaling factor was empiracally detrmined by looking at the SF of the light curve and making np.mean(PG1302_102.yerr)/math.sqrt(np.min(sfE[np.where(sfE != 0.0)])) == 1.0
		magList.append(float(words[1]))
		magerrList.append(float(words[2]))
	zipped = zip(tList, magList, magerrList, yList, yerrList)
	zipped.sort()
	tList, magList, magErrList, yList, yerrList = zip(*zipped)

	tIn = np.array(tList)
	magIn = np.array(magList)
	magerrIn = np.array(magerrList)
	yIn = np.array(yList)
	yerrIn = np.array(yerrList)
	maskIn = np.array(num*[1.0])

	return tIn, magIn, magerrIn, yIn, yerrIn, maskIn

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-b', '--band', type = str, default = 'r', help = r'SDSS band')
	parser.add_argument('-nsteps', '--nsteps', type = int, default = 1000, help = r'Number of steps per walker')
	parser.add_argument('-nwalkers', '--nwalkers', type = int, default = 25*psutil.cpu_count(logical = True), help = r'Number of walkers')
	parser.add_argument('-pMax', '--pMax', type = int, default = 1, help = r'Maximum C-AR order')
	parser.add_argument('-pMin', '--pMin', type = int, default = 1, help = r'Minimum C-AR order')
	parser.add_argument('-qMax', '--qMax', type = int, default = -1, help = r'Maximum C-MA order')
	parser.add_argument('-qMin', '--qMin', type = int, default = -1, help = r'Minimum C-MA order')
	parser.add_argument('--plot', dest = 'plot', action = 'store_true', help = r'Show plot?')
	parser.add_argument('--no-plot', dest = 'plot', action = 'store_false', help = r'Do not show plot?')
	parser.set_defaults(plot = True)
	parser.add_argument('-minT', '--minTimescale', type = float, default = 2.0, help = r'Minimum allowed timescale = minTimescale*lc.dt')
	parser.add_argument('-maxT', '--maxTimescale', type = float, default = 0.5, help = r'Maximum allowed timescale = maxTimescale*lc.T')
	parser.add_argument('-maxS', '--maxSigma', type = float, default = 2.0, help = r'Maximum allowed sigma = maxSigma*var(lc)')
	parser.add_argument('-xTol', '--xTol', type = float, default = 0.001, help = r'Relative tolerance on parameters during optimization phase')
	parser.add_argument('-maxE', '--maxEvals', type = int, default = 10000, help = r'Maximum number of evaluations per walker during optimization phase')
	parser.add_argument('--stop', dest = 'stop', action = 'store_true', help = r'Stop at end?')
	parser.add_argument('--no-stop', dest = 'stop', action = 'store_false', help = r'Do not stop at end?')
	parser.set_defaults(stop = False)
	parser.add_argument('--save', dest = 'save', action = 'store_true', help = r'Save files?')
	parser.add_argument('--no-save', dest = 'save', action = 'store_false', help = r'Do not save files?')
	parser.set_defaults(save = False)
	parser.add_argument('--fit', dest = 'fit', action = 'store_true', help = r'Fit CARMA model')
	parser.add_argument('--no-fit', dest = 'fit', action = 'store_false', help = r'Do not fit CARMA model')
	parser.set_defaults(fit = True)
	parser.add_argument('--viewer', dest = 'viewer', action = 'store_true', help = r'Visualize MCMC walkers')
	parser.add_argument('--no-viewer', dest = 'viewer', action = 'store_false', help = r'Do not visualize MCMC walkers')
	parser.set_defaults(viewer = True)
	args = parser.parse_args()

	if (args.qMax >= args.pMax):
		raise ValueError('pMax must be greater than qMax')
	if (args.qMax == -1):
		args.qMax = args.pMax - 1
	if (args.qMin == -1):
		args.qMin = 0
	if (args.pMin < 1):
		raise ValueError('pMin must be greater than or equal to 1')
	if (args.qMin < 0):
		raise ValueError('qMin must be greater than or equal to 0')

	tIn, magIn, magerrIn, yIn, yerrIn, maskIn = _readCRTS('PG 1302-102.dat', '/home/vish/code/trunk/crts/data')
	PG1302_102 = libcarma.externalLC(name = 'PG 1302-102', band = 'V', t = tIn, y = yIn, yerr = yerrIn, mask = maskIn)
	PG1302_102.yunit = r'$F$ (Jy)'
	lagsE, sfE, sferrE = PG1302_102.sf(newdt = PG1302_102.dt)

	#PG1302_102.plot()
	#PG1302_102.plotsf(newdt = PG1302_102.dt)
	#plt.show()

	PG1302_102.minTimescale = args.minTimescale
	PG1302_102.maxTimescale = args.maxTimescale
	PG1302_102.maxSigma = args.maxSigma

	taskDict = dict()
	DICDict= dict()

	for p in xrange(args.pMin, args.pMax + 1):
		for q in xrange(args.qMin, p):
			nt = libcarma.basicTask(p, q, nwalkers = args.nwalkers, nsteps = args.nsteps)

			print 'Starting libcarma fitting for p = %d and q = %d...'%(p, q)
			startLCARMA = time.time()
			nt.fit(PG1302_102)
			stopLCARMA = time.time()
			timeLCARMA = stopLCARMA - startLCARMA
			print 'libcarma took %4.3f s = %4.3f min = %4.3f hrs'%(timeLCARMA, timeLCARMA/60.0, timeLCARMA/3600.0)

			Deviances = copy.copy(nt.LnPosterior[:,args.nsteps/2:]).reshape((-1))
			DIC = 0.5*math.pow(np.std(-2.0*Deviances),2.0) + np.mean(-2.0*Deviances)
			print 'C-ARMA(%d,%d) DIC: %+4.3e'%(p, q, DIC)
			DICDict['%d %d'%(p, q)] = DIC
			taskDict['%d %d'%(p, q)] = nt

	sortedDICVals = sorted(DICDict.items(), key = operator.itemgetter(1))
	pBest = int(sortedDICVals[0][0].split()[0])
	qBest = int(sortedDICVals[0][0].split()[1])
	print 'Best model is C-ARMA(%d,%d)'%(pBest, qBest)

	bestTask = taskDict['%d %d'%(pBest, qBest)]

	loc0 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[0][0]
	loc1 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[1][0]
	
	lbls = list()
	for i in xrange(pBest):
		lbls.append(r'$\tau_{AR, %d}$ ($d$)'%(i + 1))
	for i in xrange(qBest):
		lbls.append(r'$\tau_{MA, %d} ($d$)$'%(i))
	lbls.append(r'Amp. ($F$ $d^{%2.1f}$)'%(qBest + 0.5 - pBest))
	mcmcviz.vizTriangle(pBest, qBest, bestTask.timescaleChain, labelList = lbls, figTitle = r'PG1302-102')
	
	Theta = bestTask.Chain[:, loc0, loc1]
	nt = libcarma.basicTask(pBest, qBest)
	nt.set(PG1302_102.dt, Theta)
	PG1302_102.dtSmooth = PG1302_102.dt
	nt.smooth(PG1302_102)
	PG1302_102.plot()
	bestTask.plotsf(PG1302_102, newdt = PG1302_102.dt)
	bestTask.plotacf(PG1302_102, newdt = PG1302_102.dt)

	plt.show()

	pdb.set_trace()