import numpy as np
import matplotlib.pyplot as plt
import argparse
import operator
import copy
import os
import sys
import psutil
import time
import pdb
import cPickle as pickle
import platform

from matplotlib import cm

try:
    import kali.crts
    import kali.carma
    import kali.mbhbcarma
    import kali.util.triangle
except ImportError:
    print 'kali is not setup. Setup kali by sourcing bin/setup.sh'
    sys.exit(1)

plt.ion()

parser = argparse.ArgumentParser()
parser.add_argument('-n', '--name', type=str, default='PG1302-102', help=r'Object name')
parser.add_argument('-z', '--z', type=float, default=0.278400, help=r'Redshift')
parser.add_argument('-maxEvals', '--maxEvals', type=int, default=10000,
                    help=r'Number of steps in initial maximization')
parser.add_argument('-nthreads', '--nthreads', type=int, default=psutil.cpu_count(logical=False),
                    help=r'Number of threads')
parser.add_argument('-nwalkers', '--nwalkers', type=int, default=50*psutil.cpu_count(logical=False),
                    help=r'Number of walkers')
parser.add_argument('-nsteps', '--nsteps', type=int, default=10000, help=r'Number of steps per walker')
parser.add_argument('-pMax', '--pMax', type=int, default=1, help=r'Maximum C-AR order')
parser.add_argument('-pMin', '--pMin', type=int, default=1, help=r'Minimum C-AR order')
parser.add_argument('-qMax', '--qMax', type=int, default=-1, help=r'Maximum C-MA order')
parser.add_argument('-qMin', '--qMin', type=int, default=-1, help=r'Minimum C-MA order')
parser.add_argument('--plot', dest='plot', action='store_true', help=r'Show plot?')
parser.add_argument('--no-plot', dest='plot', action='store_false', help=r'Do not show plot?')
parser.set_defaults(plot=True)
parser.add_argument('-minT', '--minTimescale', type=float, default=2.0,
                    help=r'Minimum allowed timescale = minTimescale*lc.dt')
parser.add_argument('-maxT', '--maxTimescale', type=float, default=0.5,
                    help=r'Maximum allowed timescale = maxTimescale*lc.T')
parser.add_argument('-maxS', '--maxSigma', type=float, default=2.0,
                    help=r'Maximum allowed sigma = maxSigma*var(lc)')
parser.add_argument('--stop', dest='stop', action='store_true', help=r'Stop at end?')
parser.add_argument('--no-stop', dest='stop', action='store_false', help=r'Do not stop at end?')
parser.set_defaults(stop=False)
parser.add_argument('--save', dest='save', action='store_true', help=r'Save files?')
parser.add_argument('--no-save', dest='save', action='store_false', help=r'Do not save files?')
parser.set_defaults(save=False)
parser.add_argument('--pdf', dest='pdfYN', action='store_true', help=r'Save PDF figures?')
parser.add_argument('--no-pdf', dest='pdfYN', action='store_false', help=r'Do not save PDF figures?')
parser.set_defaults(pdfYN=False)
parser.add_argument('--rerun', dest='rerun', action='store_true', help=r'Re-run fit?')
parser.add_argument('--no-rerun', dest='rerun', action='store_false', help=r'Do not re-run fit?')
parser.set_defaults(rerun=False)
args = parser.parse_args()

if not args.pdfYN:
    System = platform.system()
    if System == 'Linux':
        ext = '.jpg'
    elif System == 'Darwin':
        ext = '.png'
else:
    ext = '.pdf'

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

Obj = kali.crts.crtsLC(name=args.name, band='V', z=args.z)
Obj.minTimescale = args.minTimescale
Obj.maxTimescale = args.maxTimescale
Obj.maxSigma = args.maxSigma
Obj.dtSmooth = Obj.dt
taskDict = dict()
DICDict = dict()

outDir = os.path.join(os.environ['CRTSDATADIR'], args.name)
try:
    os.stat(outDir)
except:
    os.mkdir(outDir)


def fitCARMA(pVal, qVal, Obj, args):
    carmaTask = kali.carma.CARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                     maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                     nsteps=args.nsteps)
    print 'Starting kali.carma fitting for p = %d and q = %d...'%(pVal, qVal)
    startCARMATask = time.time()
    carmaTask.fit(Obj)
    stopCARMATask = time.time()
    timeCARMATask = stopCARMATask - startCARMATask
    print 'kali.carma took %4.3f s = %4.3f min = %4.3f hrs'%(timeCARMATask, timeCARMATask/60.0,
                                                             timeCARMATask/3600.0)
    pickle.dump(carmaTask, open(os.path.join(outDir, 'kali.carma.%d.%d.pkl'%(pVal, qVal)), 'wb'))
    return carmaTask


def fitMBHBCARMA(pVal, qVal, Obj, args):
    mbhbcarmaTask = kali.mbhbcarma.MBHBCARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                                 maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                                 nsteps=args.nsteps)
    print 'Starting kali.mbhbcarma fitting for p = %d and q = %d...'%(pVal, qVal)
    startMBHBCARMATask = time.time()
    mbhbcarmaTask.fit(Obj)
    stopMBHBCARMATask = time.time()
    timeMBHBCARMATask = stopMBHBCARMATask - startMBHBCARMATask
    print 'kali.mbhbcarma took %4.3f s = %4.3f min = %4.3f hrs'%(timeMBHBCARMATask,
                                                                 timeMBHBCARMATask/60.0,
                                                                 timeMBHBCARMATask/3600.0)
    pickle.dump(mbhbcarmaTask, open(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.pkl'%(pVal, qVal)), 'wb'))
    return mbhbcarmaTask


for pVal in xrange(args.pMin, args.pMax + 1):
    for qVal in xrange(args.qMin, args.qMax + 1):
        if args.rerun:
            carmaTask = fitCARMA(pVal, qVal, Obj, args)
        else:
            if os.path.isfile(os.path.join(outDir, 'kali.carma.%d.%d.pkl'%(pVal, qVal))):
                print 'Restoring kali.carma task with p = %d and q = %d...'%(pVal, qVal)
                carmaTask = pickle.load(open(os.path.join(outDir,
                                                          'kali.carma.%d.%d.pkl'%(pVal, qVal)), 'rb'))
            else:
                carmaTask = fitCARMA(pVal, qVal, Obj, args)
        print 'kali.carma.%d.%d DIC: %+4.3e'%(pVal, qVal, carmaTask.dic)
        taskDict['kali.carma %d %d'%(pVal, qVal)] = carmaTask
        DICDict['kali.carma %d %d'%(pVal, qVal)] = carmaTask.dic
        res = carmaTask.plottriangle()
        res[0][0].savefig(os.path.join(outDir, 'kali.carma.%d.%d.sto%s'%(pVal, qVal, ext)))
        theta_carma = carmaTask.bestTheta
        bestCarmaTask = kali.carma.CARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                             maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                             nsteps=args.nsteps)
        bestCarmaTask.set(Obj.dt, theta_carma)
        bestCarmaTask.smooth(Obj, stopT=Obj.t[-1]*1.5)
        res = Obj.plot()
        res.savefig(os.path.join(outDir, 'kali.carma.%d.%d.lc%s'%(pVal, qVal, ext)))

        if args.rerun:
            mbhbcarmaTask = fitMBHBCARMA(pVal, qVal, Obj, args)
        else:
            if os.path.isfile(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.pkl'%(pVal, qVal))):
                print 'Restoring kali.mbhbcarma task with p = %d and q = %d...'%(pVal, qVal)
                mbhbcarmaTask = pickle.load(open(os.path.join(outDir,
                                                              'kali.mbhbcarma.%d.%d.pkl'%(pVal, qVal)), 'rb'))
            else:
                mbhbcarmaTask = fitMBHBCARMA(pVal, qVal, Obj, args)
        print 'kali.mbhbcarma (%d,%d) DIC: %+4.3e'%(pVal, qVal, mbhbcarmaTask.dic)
        taskDict['kali.mbhbcarma %d %d'%(pVal, qVal)] = mbhbcarmaTask
        DICDict['kali.mbhbcarma %d %d'%(pVal, qVal)] = mbhbcarmaTask.dic
        res = mbhbcarmaTask.plottriangle()
        res[0][0].savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.sto%s'%(pVal, qVal, ext)))
        res[1][0].savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.orb%s'%(pVal, qVal, ext)))
        res[2][0].savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.aux%s'%(pVal, qVal, ext)))
        theta_mbhbcarma = mbhbcarmaTask.bestTheta
        bestMBHBCarmaTask = kali.mbhbcarma.MBHBCARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                                         maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                                         nsteps=args.nsteps)
        bestMBHBCarmaTask.set(Obj.dt, theta_mbhbcarma)
        bestMBHBCarmaTask.smooth(Obj, stopT=Obj.t[-1]*1.5)
        res = Obj.plot()
        res.savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.lc%s'%(pVal, qVal, ext)))

sortedDICVals = sorted(DICDict.items(), key=operator.itemgetter(1))
modelBest = str(sortedDICVals[0][0].split()[0])
pBest = int(sortedDICVals[0][0].split()[1])
qBest = int(sortedDICVals[0][0].split()[2])
print 'Best model is %s (%d,%d)'%(modelBest, pBest, qBest)
bestTask = taskDict['%s %d %d'%(modelBest, pBest, qBest)]
loc0 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[0][0]
loc1 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[1][0]
res = bestTask.plottriangle()
res[0][0].savefig(os.path.join(outDir, 'best_%s.%d.%d.sto%s'%(modelBest, pBest, qBest, ext)))
if modelBest == 'kali.mbhbcarma':
    res[1][0].savefig(os.path.join(outDir, 'best_%s.%d.%d.orb%s'%(modelBest, pBest, qBest, ext)))
    res[2][0].savefig(os.path.join(outDir, 'best_%s.%d.%d.aux%s'%(modelBest, pBest, qBest, ext)))
theta_best = bestTask.bestTheta
if modelBest == 'kali.carma':
    optTask = kali.carma.CARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                   maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                   nsteps=args.nsteps)
elif modelBest == 'kali.mbhbcarma':
    optTask = kali.mbhbcarma.MBHBCARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                           maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                           nsteps=args.nsteps)
optTask.set(Obj.dt, theta_best)
optTask.smooth(Obj, stopT=Obj.t[-1]*1.5)
res = Obj.plot()
res.savefig(os.path.join(outDir, 'best_%s.%d.%d.lc%s'%(modelBest, pVal, qVal, ext)))
if args.stop:
    pdb.set_trace()
