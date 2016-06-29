import math
import numpy as np
import copy
import operator
import os
import argparse
import psutil
import time
import sys
import pdb

from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV

try:
	os.environ['DISPLAY']
except KeyError:
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

try:
	import libcarma as libcarma
except ImportError:
	print 'libcarma is not setup. Setup libcarma by sourcing bin/setup.sh'
	sys.exit(1)
import util.mcmcviz as mcmcviz
from util.mpl_settings import set_plot_params
import util.triangle as triangle

fhgt = 10
fwid = 16
set_plot_params(useTex = True)

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

def plotkde(p, q, Chain, kernel = 'epanechnikov', bwidth = None, dim1 = 0, dim2 = 1, nbins = 100, pointsPerDim = 100):
	ndims = Chain.shape[0]
	nwalkers = Chain.shape[1]
	nsteps = Chain.shape[2]
	divisors = np.zeros(ndims)
	flatChain = np.swapaxes(copy.copy(Chain[:,:,nsteps/2:]).reshape((ndims,-1), order = 'F'), axis1 = 0, axis2 = 1)
	for dim in xrange(ndims):
		divisors[dim] = np.mean(flatChain[:,dim])
		flatChain[:,dim] /= divisors[dim]
	if bwidth is None:
		print 'Finding optimal bandwidth...'
		startOpt = time.time()
		grid = GridSearchCV(KernelDensity(kernel = kernel), {'bandwidth': np.linspace(0.01, 1.0, 10)}, cv=10) # 20-fold cross-validation
		grid.fit(flatChain)
		stopOpt = time.time()
		bwidth = grid.best_params_['bandwidth']
		print 'Optimal bandwidth is %f'%(bwidth)
		totTime = stopOpt - startOpt
		print 'Optimal bandwidth estimation by cross-validation took %3.2e sec = %3.2e min = %3.2e hrs'%(totTime, totTime/60.0, totTime/3600.0)
	else:
		bwidth = 0.1
	kde = KernelDensity(kernel = kernel, bandwidth = bwidth).fit(flatChain)

	xVec, yVec, Grid = _makeGrid(flatChain, dim1 = dim1, dim2 = dim2, pointsPerDim = pointsPerDim)
	Density = kde.score_samples(Grid)
	Density2D = np.reshape(Density, (-1, pointsPerDim))

	plt.figure(-5, figsize = (fwid, fwid))
	levels = MaxNLocator(nbins = nbins).tick_values(np.nanmin(Density[np.where(np.isinf(Density) == False)]), np.nanmax(Density))
	cmap = plt.get_cmap('Blues')
	norm = BoundaryNorm(levels, ncolors = cmap.N, clip = True)
	plt.contourf(xVec, yVec, Density2D, levels = levels, cmap = cmap)
	plt.colorbar()
	plt.scatter(flatChain[:,0], flatChain[:,1], marker = '.', color = '#000000')

	if dim1 < p + q:
		plt.xlabel(r'$\frac{\tau_{%d}}{%+4.3e}$'%(dim1+1, divisors[dim1]))
	elif dim1 == p + q:
		plt.xlabel(r'$\frac{\mathrm{Amp.}}{%+4.3e}$'%(divisors[dim1]))
	if dim2 < p + q:
		plt.ylabel(r'$\frac{\tau_{%d}}{%+4.3e}$'%(dim2+1, divisors[dim2]))
	elif dim2 == p + q:
		plt.ylabel(r'$\frac{\mathrm{Amp.}}{%+4.3e}$'%(divisors[dim2]))

	plt.tight_layout()
	pdb.set_trace()

def _makeGrid(flatChain, dim1 = 0, dim2 = 1, pointsPerDim = 100):
	ndims = flatChain.shape[1]
	Grid = np.zeros((int(math.pow(pointsPerDim, ndims)), ndims))
	minVal1 = 0.95*np.min(flatChain[:,dim1])
	maxVal1 = 1.05*np.max(flatChain[:,dim1])
	minVal2 = 0.95*np.min(flatChain[:,dim2])
	maxVal2 = 1.05*np.max(flatChain[:,dim2])
	xVec = np.linspace(start = minVal1, stop = maxVal1, num = pointsPerDim)
	yVec = np.linspace(start = minVal2, stop = maxVal2, num = pointsPerDim)
	for i in xrange(pointsPerDim):
		for j in xrange(pointsPerDim):
			Grid[j+i*pointsPerDim,dim1] = minVal1 + j*((maxVal1 - minVal1)/pointsPerDim)
			Grid[j+i*pointsPerDim,dim2] = minVal2 + i*((maxVal2 - minVal2)/pointsPerDim)
	return xVec, yVec, Grid

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
	parser.add_argument('--kde', dest = 'kde', action = 'store_true', help = r'Compute KDE estimate?')
	parser.add_argument('--no-kde', dest = 'kde', action = 'store_false', help = r'Do not compute KDE estimate?')
	parser.set_defaults(kde = False)
	parser.add_argument('-k', '--kernel', type = str, default = 'epanechnikov', help = r'KDE kernel to use')
	parser.add_argument('-bw', '--bandwidth', type = float, default = 0.1, help = r'KDE bandwidth to use')
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
	PG1302_102_Mag = libcarma.externalLC(name = args.name, band = args.band, t = tIn, y = magIn, yerr = magerrIn, mask = maskIn)
	PG1302_102 = libcarma.externalLC(name = args.name, band = args.band, t = tIn, y = yIn, yerr = yerrIn, mask = maskIn)
	PG1302_102_Mag.yunit = r'$m$ (V-band)'
	PG1302_102.yunit = r'$F$ (Jy)'
	PG1302_102.minTimescale = args.minTimescale
	PG1302_102.maxTimescale = args.maxTimescale
	PG1302_102.maxSigma = args.maxSigma

	if args.savefig or args.show:
		print 'Plotting %s-band magnitude light curve for %s'%(PG1302_102.band, PG1302_102.name)
		figName = os.path.join(args.pwd,'%s_magLC.jpg'%(PG1302_102_Mag.name))
		if not os.path.isfile(figName):
			PG1302_102_Mag.plot()
			if args.savefig:
				plt.savefig(figName, dpi = 1000)
			if args.show:
				plt.show()
			plt.clf()

		print 'Plotting %s-band flux light curve for %s'%(PG1302_102.band, PG1302_102.name)
		figName = os.path.join(args.pwd,'%s_LC.jpg'%(PG1302_102.name))
		if not os.path.isfile(figName):
			PG1302_102.plot()
			if args.savefig:
				plt.savefig(figName, dpi = 1000)
			if args.show:
				plt.show()
			plt.clf()

################################################## CARMA fitting #############################################

	taskDict = dict()
	DICDict= dict()
	for p in xrange(args.pMin, args.pMax + 1):
		for q in xrange(args.qMin, min(p, args.qMax)+1):
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

			if args.kde:
				plotkde(p, q, nt.timescaleChain, kernel = args.kernel, bwidth = args.bandwidth)
				#plotkde(p, q, nt.timescaleChain, kernel = args.kernel) ## Takes forever!!!!!
				if args.show:
					plt.show()
					pdb.set_trace()

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
		print 'Plotting the %s-band smoothed light curve of %s'%(PG1302_102.band, PG1302_102.name)
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

########################################## done CARMA fitting ################################################

	if args.stop:
		pdb.set_trace()