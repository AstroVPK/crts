import math as math
import numpy as np
import urllib, urllib2
import os as os
import sys as sys
import subprocess
import argparse
import matplotlib.pyplot as plt
import pdb

try:
	import libcarma as libcarma
except ImportError:
	print 'libcarma is not setup. Setup libcarma by sourcing bin/setup.sh'
	sys.exit(1)

class crtsLC(libcarma.basicLC):

	def read(self, name, band, path, **kwargs):
	
		datasplit = []
		IDs = []
		rmags = []
		runc = []
		rmjd = []
		uncf = []
		mjdf = []
		uncp = []
		magsf = []
		rblend = []
		blend = []
		t = []
		yarray = []
		yerrarray = []
		x = []
		mask = []
		
		### CODE here to open the data file ####
		#path = input('file here --> ')
		txt = open(path, 'r')
		data = txt.read()
		datasplit = data.split()
		if (datasplit[0] == "MasterID"):
			del datasplit[0:7]
			
		### CODE HERE to construct t, x, y, yerr, & mask + dt, T, startT + other properties you want to track.
		for a in range (0, len(datasplit), 7):
			IDs.append(datasplit[a])
		for b in range (1, len(datasplit), 7):	
			rmags.append(datasplit[b])
		for c in range (2, len(datasplit), 7):
			runc.append(datasplit[c])
		for d in range (5, len(datasplit), 7):
			rmjd.append(datasplit[d])
		for e in range (6, len(datasplit), 7):
			rblend.append(datasplit[e])
			
		for f in range (len(rblend)):
			x.append(0)
			
		for l in range (len(rmags)):
			a = float(rmags[l])
			magsf.append(a)
		for m in range (len(runc)):
			b = float(runc[m])
			uncf.append(b)
		for n in range (len(rmjd)):
			c = float(rmjd[n])
			mjdf.append(c)
		for o in range (len(rblend)):
			d = int(rblend[o])
			blend.append(d)
			
		for u in range (len(mjdf)):
			t.append(mjdf[u]-mjdf[0])
		for v in range (len(magsf)):
			y, yerr = libcarma.pogsonFlux(magsf[v], uncf[v])
			yarray.append(y)
			yerrarray.append(yerr)
		
		for w in range (len(blend)):
			if blend[w] == 0:
				mask.append(1)
			else:
				mask.append(0)
		
		self.mask = np.require(np.array(mask))
		self.t = np.require(np.array(t))
		self.yarray = np.require(np.array(yarray))
		self.yerrarray = np.require(np.array(yerrarray))
		self.x = np.require(np.array(x))
		dt = t[1] - t[0]
		
		'''print mask[:]
		print '-'*20
		print t[:]
		print '-'*20
		print yarray[:]
		print '-'*20
		print yerrarray[:]
		print '-'*20
		print x[:]'''
	

		### Boilerplate follows - you don't have to mess with it
		self._computedCadenceNum = -1
		self._tolIR = 1.0e-3
		self._fracIntrinsicVar = 0.0
		self._fracNoiseToSignal = 0.0
		self._maxSigma = 2.0
		self._minTimescale = 2.0
		self._maxTimescale = 0.5
		self._pSim = 0
		self._qSim = 0
		self._pComp = 0
		self._qComp = 0
		self._isSmoothed = False ## Has the LC been smoothed?
		self._dtSmooth = 0.0
		self._isRegular = True
		self.XSim = np.require(np.zeros(self.pSim), requirements=['F', 'A', 'W', 'O', 'E']) ## State of light curve at last timestamp
		self.PSim = np.require(np.zeros(self.pSim*self.pSim), requirements=['F', 'A', 'W', 'O', 'E']) ## Uncertainty in state of light curve at last timestamp.
		self.XComp = np.require(np.zeros(self.pComp), requirements=['F', 'A', 'W', 'O', 'E']) ## State of light curve at last timestamp
		self.PComp = np.require(np.zeros(self.pComp*self.pComp), requirements=['F', 'A', 'W', 'O', 'E']) ## Uncertainty in state of light curve at last timestamp.
		self._name = str(name) ## The name of the light curve (usually the object's name).
		self._band = str(r'V') ## The name of the photometric band (eg. HSC-I or SDSS-g etc..).
		self._xunit = r'$d$' ## Unit in which time is measured (eg. s, sec, seconds etc...).
		#self._yunit = r'who the f*** knows?' ## Unit in which the flux is measured (eg Wm^{-2} etc...).
		self._yunit = r'$F$ (Jy)' ## Unit in which the flux is measured (eg Wm^{-2} etc...).

		count = int(np.sum(self.mask))
		y_meanSum = 0.0
		yerr_meanSum = 0.0
		for i in xrange(self.numCadences):
			y_meanSum += self.mask[i]*self.y[i]
			yerr_meanSum += self.mask[i]*self.yerr[i]
		if count > 0.0:
			self._mean = y_meanSum/count
			self._meanerr = yerr_meanSum/count
		else:
			self._mean = 0.0
			self._meanerr = 0.0
		y_stdSum = 0.0
		yerr_stdSum = 0.0
		for i in xrange(self.numCadences):
			y_stdSum += math.pow(self.mask[i]*self.y[i] - self._mean, 2.0)
			yerr_stdSum += math.pow(self.mask[i]*self.yerr[i] - self._meanerr, 2.0)
		if count > 0.0:
			self._std = math.sqrt(y_stdSum/count)
			self._stderr = math.sqrt(yerr_stdSum/count)
		else:
			self._std = 0.0
			self._stderr = 0.0

	def write(self, name, path = None, **kwrags):
		pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', '--name', type = str, default = 'PG1302102', help = r'Object name')
	args = parser.parse_args()

	LC = crtsLC(name = args.name, band = 'V')

	LC.plot()
	LC.plotacf()
	LC.plotsf()
	plt.show(False)