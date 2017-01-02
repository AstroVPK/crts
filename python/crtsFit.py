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
    import kali.mbhb
    import kali.mbhbcarma
    import kali.util.triangle
except ImportError:
    print 'kali is not setup. Setup kali by sourcing bin/setup.sh'
    sys.exit(1)

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
parser.add_argument('-minT', '--minTimescale', type=float, default=2.0,
                    help=r'Minimum allowed timescale = minTimescale*lc.dt')
parser.add_argument('-maxT', '--maxTimescale', type=float, default=0.5,
                    help=r'Maximum allowed timescale = maxTimescale*lc.T')
parser.add_argument('-maxS', '--maxSigma', type=float, default=2.0,
                    help=r'Maximum allowed sigma = maxSigma*var(lc)')
parser.add_argument('-dpi', '--dpi', type=int, default=300,
                    help=r'Dots Per inch to save raster figures at?')
parser.add_argument('--plot', dest='plot', action='store_true', help=r'Show plot?')
parser.add_argument('--no-plot', dest='plot', action='store_false', help=r'Do not show plot?')
parser.set_defaults(plot=True)
parser.add_argument('--stop', dest='stop', action='store_true', help=r'Stop at end?')
parser.add_argument('--no-stop', dest='stop', action='store_false', help=r'Do not stop at end?')
parser.set_defaults(stop=False)
parser.add_argument('--save', dest='save', action='store_true', help=r'Save files?')
parser.add_argument('--no-save', dest='save', action='store_false', help=r'Do not save files?')
parser.set_defaults(save=True)
parser.add_argument('--pdf', dest='pdfYN', action='store_true', help=r'Save PDF figures?')
parser.add_argument('--no-pdf', dest='pdfYN', action='store_false', help=r'Do not save PDF figures?')
parser.set_defaults(pdfYN=False)
parser.add_argument('--rerun', dest='rerun', action='store_true', help=r'Re-run fit?')
parser.add_argument('--no-rerun', dest='rerun', action='store_false', help=r'Do not re-run fit?')
parser.set_defaults(rerun=False)
parser.add_argument('--ion', dest='ion', action='store_true', help=r'Display plots interactively?')
parser.add_argument('--no-ion', dest='ion', action='store_false', help=r'Do not Display plots interactively?')
parser.set_defaults(ion=False)
args = parser.parse_args()

if args.ion:
    plt.ion()

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

outDir = os.path.join(os.environ['CRTSDATADIR'], args.name)
try:
    os.stat(outDir)
except:
    os.mkdir(outDir)

Obj = kali.crts.crtsLC(name=args.name, band='V', z=args.z)
Obj.minTimescale = args.minTimescale
Obj.maxTimescale = args.maxTimescale
Obj.maxSigma = args.maxSigma
Obj.dtSmooth = 0.5
taskDict = dict()
DICDict = dict()

if args.plot:
    basic = Obj.plot(fig=100, colory=r'#000000')
    if args.save:
        basic.savefig(os.path.join(outDir, 'kali.lc%s'%(ext)), dpi=args.dpi)


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


def fitMBHB(Obj, args):
    mbhbTask = kali.mbhb.MBHBTask(nthreads=args.nthreads,
                                  maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                  nsteps=args.nsteps)
    print 'Starting kali.mbhb fitting for p = %d and q = %d...'%(0, 0)
    startMBHBTask = time.time()
    mbhbTask.fit(Obj)
    stopMBHBTask = time.time()
    timeMBHBTask = stopMBHBTask - startMBHBTask
    print 'kali.mbhb took %4.3f s = %4.3f min = %4.3f hrs'%(timeMBHBTask,
                                                            timeMBHBTask/60.0,
                                                            timeMBHBTask/3600.0)
    pickle.dump(mbhbTask, open(os.path.join(outDir, 'kali.mbhb.%d.%d.pkl'%(0, 0)), 'wb'))
    return mbhbTask


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


if args.rerun:
    mbhbTask = fitMBHB(Obj, args)
else:
    if os.path.isfile(os.path.join(outDir, 'kali.mbhb.%d.%d.pkl'%(0, 0))):
        print 'Restoring kali.mbhb task with p = %d and q = %d...'%(0, 0)
        mbhbTask = pickle.load(open(os.path.join(outDir,
                                                 'kali.mbhb.%d.%d.pkl'%(0, 0)), 'rb'))
    else:
        mbhbTask = fitMBHB(Obj, args)
print 'kali.mbhb (%d,%d) DIC: %+4.3e'%(0, 0, mbhbTask.dic)
if args.plot:
    res = mbhbTask.plottriangle()
    if args.save:
        res[0][0].savefig(os.path.join(outDir, 'kali.mbhb.%d.%d.orb%s'%(0, 0, ext)),
                          dpi=args.dpi)
        res[1][0].savefig(os.path.join(outDir, 'kali.mbhb.%d.%d.aux%s'%(0, 0, ext)),
                          dpi=args.dpi)
taskDict['kali.mbhb %d %d'%(0, 0)] = mbhbTask
DICDict['kali.mbhb %d %d'%(0, 0)] = mbhbTask.dic
theta_mbhb = mbhbTask.bestTheta
bestMBHBTask = kali.mbhb.MBHBTask(nthreads=args.nthreads,
                                  maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                  nsteps=args.nsteps)
bestMBHBTask.set(theta_mbhb)
bestMBHBTask.smooth(Obj, stopT=(Obj.t[-1] + Obj.T*0.5))
if args.plot:
    res = Obj.plot(colory=r'#000000', colors=[r'#7b3294', r'#c2a5cf'])
    if args.save:
        res.savefig(os.path.join(outDir, 'kali.mbhb.%d.%d.lc%s'%(0, 0, ext)), dpi=args.dpi)

if args.plot:
    comp = Obj.plot(fig=100, colory=r'#000000', colors=[r'#7b3294', r'#c2a5cf'])

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
        print 'kali.carma (%d,%d) DIC: %+4.3e'%(pVal, qVal, carmaTask.dic)
        if args.plot:
            res = carmaTask.plottriangle()
            if args.save:
                res[0][0].savefig(os.path.join(outDir, 'kali.carma.%d.%d.sto%s'%(pVal, qVal, ext)),
                                  dpi=args.dpi)
        taskDict['kali.carma %d %d'%(pVal, qVal)] = carmaTask
        DICDict['kali.carma %d %d'%(pVal, qVal)] = carmaTask.dic
        theta_carma = carmaTask.bestTheta
        bestCarmaTask = kali.carma.CARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                             maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                             nsteps=args.nsteps)
        bestCarmaTask.set(Obj.dt, theta_carma)
        bestCarmaTask.smooth(Obj, stopT=(Obj.t[-1] + Obj.T*0.5))
        if args.plot:
            res = Obj.plot(colory=r'#000000', colors=[r'#a6611a', r'#dfc27d'])
            if args.save:
                res.savefig(os.path.join(outDir, 'kali.carma.%d.%d.lc%s'%(pVal, qVal, ext)), dpi=args.dpi)

        if args.plot:
            comp = Obj.plot(fig=100, clearFig=False, colory=r'#000000', colors=[r'#a6611a', r'#dfc27d'])

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
        if args.plot:
            res = mbhbcarmaTask.plottriangle()
            if args.save:
                res[0][0].savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.sto%s'%(pVal, qVal, ext)),
                                  dpi=args.dpi)
                res[1][0].savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.orb%s'%(pVal, qVal, ext)),
                                  dpi=args.dpi)
                res[2][0].savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.aux%s'%(pVal, qVal, ext)),
                                  dpi=args.dpi)
        taskDict['kali.mbhbcarma %d %d'%(pVal, qVal)] = mbhbcarmaTask
        DICDict['kali.mbhbcarma %d %d'%(pVal, qVal)] = mbhbcarmaTask.dic
        theta_mbhbcarma = mbhbcarmaTask.bestTheta
        bestMBHBCarmaTask = kali.mbhbcarma.MBHBCARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                                         maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                                         nsteps=args.nsteps)
        bestMBHBCarmaTask.set(Obj.dt, theta_mbhbcarma)
        bestMBHBCarmaTask.smooth(Obj, stopT=(Obj.t[-1] + Obj.T*0.5))
        if args.plot:
            res = Obj.plot(colory=r'#000000', colors=[r'#018571', r'#80cdc1'])
            if args.save:
                res.savefig(os.path.join(outDir, 'kali.mbhbcarma.%d.%d.lc%s'%(pVal, qVal, ext)), dpi=args.dpi)

        if args.plot:
            comp = Obj.plot(fig=100, clearFig=False, colory=r'#000000', colors=[r'#018571', r'#80cdc1'])
        # beamObj = bestMBHBCarmaTask.beam(duration=1.5*Obj.T, startT=Obj.startT)
        # if args.plot:
            # comp = beamObj.plot(fig=100, clearFig=False, colorx=r'#ca0020')
        if args.plot:
            if args.save:
                comp.savefig(os.path.join(outDir, 'comp_kali.%d.%d.lc%s'%(pVal, qVal, ext)), dpi=args.dpi)

sortedDICVals = sorted(DICDict.items(), key=operator.itemgetter(1))
modelBest = str(sortedDICVals[0][0].split()[0])
pBest = int(sortedDICVals[0][0].split()[1])
qBest = int(sortedDICVals[0][0].split()[2])
print 'Best model is %s (%d,%d)'%(modelBest, pBest, qBest)
bestTask = taskDict['%s %d %d'%(modelBest, pBest, qBest)]
loc0 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[0][0]
loc1 = np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[1][0]
if args.plot:
    res = bestTask.plottriangle()
    if args.save:
        res[0][0].savefig(os.path.join(outDir, 'best_%s.%d.%d.sto%s'%(modelBest, pBest, qBest, ext)),
                          dpi=args.dpi)
        if modelBest == 'kali.mbhbcarma':
            res[1][0].savefig(os.path.join(outDir, 'best_%s.%d.%d.orb%s'%(modelBest, pBest, qBest, ext)),
                              dpi=args.dpi)
            res[2][0].savefig(os.path.join(outDir, 'best_%s.%d.%d.aux%s'%(modelBest, pBest, qBest, ext)),
                              dpi=args.dpi)
theta_best = bestTask.bestTheta
if modelBest == 'kali.mbhb':
    optTask = kali.mbhb.MBHBTask(nthreads=args.nthreads,
                                 maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                 nsteps=args.nsteps)
    optTask.set(theta_best)
elif modelBest == 'kali.carma':
    optTask = kali.carma.CARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                   maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                   nsteps=args.nsteps)
    optTask.set(Obj.dt, theta_best)
elif modelBest == 'kali.mbhbcarma':
    optTask = kali.mbhbcarma.MBHBCARMATask(p=pVal, q=qVal, nthreads=args.nthreads,
                                           maxEvals=args.maxEvals, nwalkers=args.nwalkers,
                                           nsteps=args.nsteps)
    optTask.set(Obj.dt, theta_best)
optTask.smooth(Obj, stopT=(Obj.t[-1] + Obj.T*0.5))
if args.plot:
    res = Obj.plot()
    if args.save:
        res.savefig(os.path.join(outDir, 'best_%s.%d.%d.lc%s'%(modelBest, pVal, qVal, ext)), dpi=args.dpi)
if args.stop:
    pdb.set_trace()
