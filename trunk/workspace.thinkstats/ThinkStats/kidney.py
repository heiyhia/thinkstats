"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import numpy
import random
import sys

import Cdf
import correlation
import myplot
import matplotlib.pyplot as pyplot

INTERVAL = 245/365.0

def DoTheMath():
    interval = 3291.0
    dt = 811.0 * 3

    doublings = interval / dt

    d0 = 15.5
    d1 = d0 / math.pow(2.0, doublings)

    print 'interval', interval
    print 'dt', dt
    print 'doublings', doublings
    print 'd0', d0
    print 'd1', d1

    def log2(x):
        return math.log(x) / math.log(2)

    d0 = 0.1
    d1 = 15.5

    doublings = log2(15.5 / 0.1)
    dt = interval/doublings

    print 'doublings', doublings
    print 'dt', dt

    vdt = dt / 3
    rdt = 365/vdt

    print 'vdt', vdt
    print 'rdt', rdt


def MakeCdf():
    """Use the data from Zhang et al. to construct a CDF."""
    n = 53.0
    freqs = [0, 2, 31, 42, 48, 51, 52, 53]
    ps = [freq/n for freq in freqs]
    xs = numpy.arange(-1.5, 6.5, 1.0)

    cdf = Cdf.Cdf(xs, ps)
    return cdf


def PlotCdf(cdf, fit):
    """Plots the actual and fitted distributions."""
    xs, ps = cdf.xs, cdf.ps
    cps = [1-p for p in ps]

    # CCDF on logy scale: shows exponential behavior
    myplot.Plot(xs, cps, 'bo-',
                root = 'kidney1',
                xlabel='RDT',
                ylabel='log CCDF',
                yscale='log',
                )

    # CDF, model and data
    myplot.Cdf(fit,
               axis=[-2, 7, 0, 1])

    myplot.Plot(xs, ps, 'gs',
                clf=False,
                root = 'kidney2',
                xlabel='RDT',
                ylabel='CDF',
                axis=[-2, 7, 0, 1])


def QQPlot(cdf, fit):
    """Makes a QQPlot of the values from actual and fitted distributions."""
    xs = [-1.5, 5.5]
    myplot.Plot(xs, xs, 'b-')

    xs, ps = cdf.xs, cdf.ps
    fs = [fit.Value(p) for p in ps]

    myplot.Plot(xs, fs, 'gs',
                clf=False,
                root = 'kidney3',
                xlabel='Actual',
                ylabel='Model')
    

def FitCdf(cdf):
    """Fits a line to the log CCDF and returns the slope."""
    xs, ps = cdf.xs, cdf.ps
    cps = [1-p for p in ps]

    xs = xs[1:-1]
    lcps = [math.log(p) for p in cps[1:-1]]
    
    inter, slope = correlation.LeastSquares(xs, lcps)
    return -slope


def GenerateRdt(p, lam1, lam2):
    """Generate an RDT from a mixture of exponential distributions.

    With prob p, generate a negative value with param lam2;
    otherwise generate a positive value with param lam1.
    """
    if random.random() < p:
        return -random.expovariate(lam2)
    else:
        return random.expovariate(lam1)


def GenerateSample(n, p, lam1, lam2):
    """Generates a sample of RDTs."""
    xs = [GenerateRdt(p, lam1, lam2) for i in xrange(n)]
    return xs


def GenerateCdf(n=10000, p=0.35, lam1=0.79, lam2=5.0):
    """Generates a sample of RDTs and returns its CDF."""
    xs = GenerateSample(n, p, lam1, lam2)
    cdf = Cdf.MakeCdfFromList(xs)
    return cdf


def Step(size, rdt, interval):
    """Computes the final volume of a tumor with given size and rdt.

    size: initial volume in mL (cm^3)
    rdt: reciprocal doubling time in doublings per year
    interval in years
    """
    doublings = rdt * interval
    final = size * 2**doublings
    return final


def RandomStep(size, interval):
    """Chooses a random RDT and returns tumor size at end of interval."""
    rdt = GenerateRdt(p=0.35, lam1=0.79, lam2=5.0)
    final = Step(size, rdt, interval)
    return final


"""Cache maps from size bin to a list of sequences that could be
observed in that bin.
"""
cache = {}

def BinToCm(y, factor=10):
    """Computes the linear dimension for a given bin."""
    return math.exp(y / factor)

def CmToBin(x, factor=10):
    """Computes the bin for a given linear dimension."""
    return factor * math.log(x)


def LinearMeasure(volume, exp=1.0/3.0):
    """Converts a colume to a linear measure."""
    return volume ** exp


def AddToCache(final, seq):
    """Adds a sequence to the bin the corresponds to final."""
    cm = LinearMeasure(final)
    bin = round(CmToBin(cm))
    cache.setdefault(bin, []).append(seq)


def ExtendSequence(t, interval):
    """Generates a new random value and adds it to the end of t.

    Side-effect: adds sub-sequences to the cache.

    t: sequence of values so far
    interval: timestep in years
    """
    initial = t[-1]
    final = RandomStep(initial, interval)
    res = t + (final,)
    AddToCache(final, res)
    
    return final, res

def MakeSequence(n=60, v0=0.1, interval=INTERVAL):
    """Simulate the growth of a tumor.

    n: number of time steps
    v0: initial volume in mL (cm^3)
    interval: timestep in years
    """
    vs = v0,
    time = 0
    times = [time]

    for i in range(n):
        time += interval
        times.append(time)

        final, vs = ExtendSequence(vs, interval)

    return times, vs


def MakeSequences(n=10):
    """Returns a sequence of times and a list of sequences of volumes."""
    sequences = []
    for i in range(n):
        ts, vs = MakeSequence()
        sequences.append(vs)

    return ts, sequences


def PrintCache():
    """Prints the size (cm) for each bin, and the number of sequences."""
    for bin, ss in sorted(cache.iteritems()):
        size = BinToCm(bin)
        print bin, size, len(ss)
        

def PlotSequence(ts, vs, color='blue'):
    """Plots a time series of linear measurements.

    ts: sequence of times in years
    vs: sequence of columes
    color: color string
    """
    line_options = dict(color=color, linewidth=1, alpha=0.2)
    xs = [v**(1.0/3) for v in vs]
    myplot.Plot(ts, xs,
                line_options=line_options,
                yscale='log',
                clf=False)

def PlotSequences(ts, ss):
    """Plots linear measurement vs time.

    ts: sequence of times
    ss: list of sequences of volumes
    """
    for vs in ss:
        PlotSequence(ts, vs)

    myplot.Save(root='kidney4',
                xlabel='years',
                ylabel='log size cm')
    pyplot.clf()


def PlotBin(bin, color='blue'):
    "Plots the set of sequences for the given bin."""
    ss = cache[bin]
    for vs in ss:
        n = len(vs)
        age = n * INTERVAL
        ts = numpy.linspace(-age, 0, n)
        PlotSequence(ts, vs, color)


def PlotCache():
    """Plots the set of sequences for each bin."""
    # 2.01, 4.95 cm, 9.97 cm, 14.879 cm
    bins = [7.0, 16.0, 23.0, 27.0]
    colors = ['blue', 'green', 'red', 'cyan']
    cdfs = []

    for bin, color in zip(bins, colors):
        PlotBin(bin, color)

    myplot.Save(root='kidney5',
                xlabel='years',
                ylabel='log size cm')


def CdfBin(bin, name=''):
    """Forms the cdf of ages for the sequences in this bin."""
    ss = cache[bin]
    ages = []
    for vs in ss:
        n = len(vs)
        age = n * INTERVAL
        ages.append(age)

    cdf = Cdf.MakeCdfFromList(ages, name=name)
    return cdf


def CdfCache():
    """Plots the cdf of ages for each bin."""
    # 2.01, 4.95 cm, 9.97 cm, 14.879 cm
    bins = [7.0, 16.0, 23.0, 27.0]
    names = ['2 cm', '5 cm', '10 cm', '15 cm']
    cdfs = []

    for bin, name in zip(bins, names):
        cdf = CdfBin(bin, name)
        cdfs.append(cdf)
        print name, cdf.Percentile(50), cdf.Percentile(5), cdf.Percentile(95)

    myplot.Cdfs(cdfs,
                root='kidney6',
                xlabel='years',
                ylabel='CDF')


def main(script):
    random.seed(17)
    ts, ss = MakeSequences(100)
    PlotSequences(ts, ss)
    PlotCache()

    ts, ss = MakeSequences(900)
    CdfCache()

    #PrintCache()

    return

    cdf = MakeCdf()

    lam1 = FitCdf(cdf)
    fit = GenerateCdf(lam1=lam1)
    PlotCdf(cdf, fit)
    QQPlot(cdf, fit)

    
if __name__ == '__main__':
    main(*sys.argv)


