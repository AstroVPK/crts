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

def _writeResult(filename, path, p, q, DIC, Chain, timescaleChain, LnPosterior):
	ndims = p + q + 1
	nwalkers = Chain.shape[1]
	nsteps = Chain.shape[2]
	bestWalker = np.where(LnPosterior == np.max(LnPosterior))[0][0]
	bestStep = np.where(LnPosterior == np.max(LnPosterior))[1][0]

	fullFilename = filename + '_%d_%d_Chain.dat'%(p, q)
	filePath = os.path.join(path, fullFilename)
	with open(filePath, 'w') as fileOut:
		line = 'p: %d q: %d ndims: %d nwalkers: %d nsteps: %d DIC: %+16.15e\n'%(p, q, ndims, nwalkers, nsteps, DIC)
		fileOut.write(line)
		line = 'MLE\n'
		fileOut.write(line)
		line = ''
		for dimNum in xrange(ndims):
			line += '%+16.15e '%(Chain[dimNum, bestWalker, bestStep])
		line += '%+16.15e\n'%(LnPosterior[bestWalker, bestStep])
		fileOut.write(line)
		line = 'Draws\n'
		fileOut.write(line)
		for stepNum in xrange(nsteps):
			for walkerNum in xrange(nwalkers):
				line = ''
				for dimNum in xrange(ndims):
					line += '%+16.15e '%(Chain[dimNum, walkerNum, stepNum])
				line += '%+16.15e\n'%(LnPosterior[walkerNum, stepNum])
				fileOut.write(line)

	fullFilename = filename + '_%d_%d_timescaleChain.dat'%(p, q)
	filePath = os.path.join(path, fullFilename)
	with open(filePath, 'w') as fileOut:
		line = 'p: %d q: %d ndims: %d nwalkers: %d nsteps: %d DIC: %+16.15e\n'%(p, q, ndims, nwalkers, nsteps, DIC)
		fileOut.write(line)
		line = 'MLE\n'
		fileOut.write(line)
		line = ''
		for dimNum in xrange(ndims):
			line += '%+16.15e '%(timescaleChain[dimNum, bestWalker, bestStep])
		line += '%+16.15e\n'%(LnPosterior[bestWalker, bestStep])
		fileOut.write(line)
		line = 'Draws\n'
		fileOut.write(line)
		for stepNum in xrange(nsteps):
			for walkerNum in xrange(nwalkers):
				line = ''
				for dimNum in xrange(ndims):
					line += '%+16.15e '%(timescaleChain[dimNum, walkerNum, stepNum])
				line += '%+16.15e\n'%(LnPosterior[walkerNum, stepNum])
				fileOut.write(line)

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
	parser.add_argument('-n', '--name', type = str, default = 'PG1302-102', help = r'CRTS Object name')
	parser.add_argument('-b', '--band', type = str, default = 'V', help = r'CRTS band')
	parser.add_argument('-pwd', '--pwd', type = str, default = os.environ['CRTSDATADIR'], help = r'CRTS data directory')
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
	parser.add_argument('--save', dest = 'save', action = 'store_true', help = r'Save MCMC files?')
	parser.add_argument('--no-save', dest = 'save', action = 'store_false', help = r'Do not save MCMC files?')
	parser.set_defaults(save = True)
	parser.add_argument('--fit', dest = 'fit', action = 'store_true', help = r'Fit CARMA model')
	parser.add_argument('--no-fit', dest = 'fit', action = 'store_false', help = r'Do not fit CARMA model')
	parser.set_defaults(fit = True)
	parser.add_argument('--viewer', dest = 'viewer', action = 'store_true', help = r'Visualize MCMC walkers')
	parser.add_argument('--no-viewer', dest = 'viewer', action = 'store_false', help = r'Do not visualize MCMC walkers')
	parser.set_defaults(viewer = True)
	parser.add_argument('--show', dest = 'show', action = 'store_true', help = r'Show figures?')
	parser.add_argument('--no-show', dest = 'show', action = 'store_false', help = r'Do not show figures')
	parser.set_defaults(show = False)
	parser.add_argument('--savefig', dest = 'savefig', action = 'store_true', help = r'Save figures?')
	parser.add_argument('--no-savefig', dest = 'savefig', action = 'store_false', help = r'Do not save figures')
	parser.set_defaults(savefig = True)
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

	tIn, magIn, magerrIn, yIn, yerrIn, maskIn = _readCRTS(args.name + '.dat', args.pwd)
	PG1302_102 = libcarma.externalLC(name = args.name, band = args.band, t = tIn, y = yIn, yerr = yerrIn, mask = maskIn)
	PG1302_102.yunit = r'$F$ (Jy)'
	PG1302_102.minTimescale = args.minTimescale
	PG1302_102.maxTimescale = args.maxTimescale
	PG1302_102.maxSigma = args.maxSigma

	if args.savefig or args.show:
		print 'Plotting light curve for %s'%(PG1302_102.name)
		figName = os.path.join(args.pwd,'%s_LC.jpg'%(PG1302_102.name))
		if not os.path.isfile(figName):
			PG1302_102.plot()
			if args.savefig:
				plt.savefig(figName, dpi = 1000)
			if args.show:
				plt.show()
			plt.clf()

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

			if args.save:
				_writeResult(args.name, args.pwd, p, q, DIC, nt.Chain, nt.timescaleChain, nt.LnPosterior)

	sortedDICVals = sorted(DICDict.items(), key = operator.itemgetter(1))
	pBest = int(sortedDICVals[0][0].split()[0])
	qBest = int(sortedDICVals[0][0].split()[1])
	print 'Best model is C-ARMA(%d,%d)'%(pBest, qBest)

	bestTask = taskDict['%d %d'%(pBest, qBest)]

	loc0 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[0][0]
	loc1 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[1][0]

	if args.savefig or args.show:
		lblsTau = list()
		for i in xrange(pBest):
			lblsTau.append(r'$\tau_{AR, %d}$ ($d$)'%(i + 1))
		for i in xrange(qBest):
			lblsTau.append(r'$\tau_{MA, %d}$ ($d$)'%(i))
		lblsTau.append(r'Amp. ($Jy$ $d^{%2.1f}$)'%(qBest + 0.5 - pBest))
		try:
			mcmcviz.vizTriangle(pBest, qBest, bestTask.timescaleChain, labelList = lblsTau, figTitle = r'PG1302-102 Timescales Chain')
		except ValueError as err:
			print str(err)
		else:
			print 'Plotting triangle plot of timescales for the %s-band light curve of %s'%(PG1302_102.band, PG1302_102.name)
			figName = os.path.join(args.pwd,'%s_Tau.jpg'%(PG1302_102.name))
			if not os.path.isfile(figName):
				if args.savefig:
					plt.savefig(figName, dpi = 1000)
				if args.show:
					plt.show()
				plt.clf()

		lblsTheta = list()
		for i in xrange(pBest):
			lblsTheta.append(r'$\alpha_{%d}$'%(i + 1))
		for i in xrange(qBest + 1):
			lblsTheta.append(r'$\beta_{%d}$'%(i))
		try:
			mcmcviz.vizTriangle(pBest, qBest, bestTask.Chain, labelList = lblsTheta, figTitle = r'PG1302-102 Coefficients Chain')
		except ValueError as err:
			print str(err)
		else:
			print 'Plotting triangle plot of parameters for the %s-band light curve of %s'%(PG1302_102.band, PG1302_102.name)
			figName = os.path.join(args.pwd,'%s_Theta.jpg'%(PG1302_102.name))
			if not os.path.isfile(figName):
				if args.savefig:
					plt.savefig(figName, dpi = 1000)
				if args.show:
					plt.show()
				plt.clf()

	if args.savefig or args.show:
		print 'Plotting the %s-band light curve of %s'%(PG1302_102.band, PG1302_102.name)
		figName = os.path.join(args.pwd,'%s_smoothLC.jpg'%(PG1302_102.name))
		if not os.path.isfile(figName):
			Theta = bestTask.Chain[:, loc0, loc1]
			nt = libcarma.basicTask(pBest, qBest)
			nt.set(PG1302_102.dt, Theta)
			PG1302_102.dtSmooth = PG1302_102.dt
			nt.smooth(PG1302_102)
			PG1302_102.plot()
			if args.savefig:
				plt.savefig(figName, dpi = 1000)
			if args.show:
				plt.show()
			plt.clf()


		print 'Plotting the structure function of the %s-band light curve of %s'%(PG1302_102.band, PG1302_102.name)
		figName = os.path.join(args.pwd,'%s_SF.jpg'%(PG1302_102.name))
		if not os.path.isfile(figName):
			bestTask.plotsf(PG1302_102, newdt = PG1302_102.dt)
			if args.savefig:
				plt.savefig(figName, dpi = 1000)
			if args.show:
				plt.show()
			plt.clf()

	if args.stop:
		pdb.set_trace()