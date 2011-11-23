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
import continuous
import correlation
import erf
import myplot
import matplotlib.pyplot as pyplot

import cProfile

INTERVAL = 245/365.0
FORMATS = ['png', 'pdf']
MINSIZE = 0.2

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
                root='kidney1',
                formats=FORMATS,
                xlabel='RDT',
                ylabel='CCDF (log scale)',
                yscale='log',
                )

    # CDF, model and data
    #myplot.Cdf(fit,
    #           axis=[-2, 7, 0, 1])

    mxs, mys = ModelCdf()
    myplot.Plot(mxs, mys, 'b-', 
                line_options=dict(label='model'))

    myplot.Plot(xs, ps, 'gs',
                line_options=dict(label='data'),
                clf=False,
                root='kidney2',
                formats=FORMATS,
                xlabel='RDT (volume doublings per year)',
                ylabel='CDF',
                title='Distribution of RDT',
                axis=[-2, 7, 0, 1],
                loc=4)


def QQPlot(cdf, fit):
    """Makes a QQPlot of the values from actual and fitted distributions."""
    xs = [-1.5, 5.5]
    myplot.Plot(xs, xs, 'b-')

    xs, ps = cdf.xs, cdf.ps
    fs = [fit.Value(p) for p in ps]

    myplot.Plot(xs, fs, 'gs',
                clf=False,
                root = 'kidney3',
                formats=FORMATS,
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


def CorrelatedGenerator(cdf, rho):
    """Generates a sequence of values from cdf with correlation.

    Generates a correlated standard normal series, then transforms to
    values from cdf

    cdf: distribution to choose from
    rho: target coefficient of correlation
    """
    def Transform(x):
        p = erf.NormalCdf(x)
        y = cdf.Value(p)
        return y

    x = random.gauss(0, 1)
    yield Transform(x)

    sigma = math.sqrt(1 - rho**2);    
    while True:
        x = random.gauss(x * rho, sigma)
        yield Transform(x)


def UncorrelatedGenerator(cdf, rho):
    """Generates a sequence of values from cdf with correlation.

    Generates a correlated standard normal series, then transforms to
    values from cdf

    cdf: distribution to choose from
    rho: target coefficient of correlation
    """
    while True:
        x = cdf.Random()
        yield x


def RdtGenerator(n, rho, cdf):
    """Returns an iterator with n values from cdf and the given correlation.

    n: number of elements
    rho: coefficient of correlation
    cdf: Cdf object
    """
    if rho == 0.0:
        return UncorrelatedGenerator(cdf, rho)
    else:
        return CorrelatedGenerator(cdf, rho)


def GenerateRdt(pc, lam1, lam2):
    """Generate an RDT from a mixture of exponential distributions.

    With prob p, generate a negative value with param lam2;
    otherwise generate a positive value with param lam1.
    """
    if random.random() < pc:
        return -random.expovariate(lam2)
    else:
        return random.expovariate(lam1)


def GenerateSample(n, pc, lam1, lam2):
    """Generates a sample of RDTs."""
    xs = [GenerateRdt(pc, lam1, lam2) for i in xrange(n)]
    return xs


def GenerateCdf(n=1000, pc=0.35, lam1=0.79, lam2=5.0):
    """Generates a sample of RDTs and returns its CDF."""
    xs = GenerateSample(n, pc, lam1, lam2)
    cdf = Cdf.MakeCdfFromList(xs)
    return cdf


def ModelCdf(pc=0.35, lam1=0.79, lam2=5.0):
    x1 = numpy.arange(-2, 0, 0.1)
    y1 = [pc * (1 - continuous.ExpoCdf(-x, lam2)) for x in x1]
    x2 = numpy.arange(0, 7, 0.1)
    y2 = [pc + (1-pc) * continuous.ExpoCdf(x, lam1) for x in x2]
    return list(x1) + list(x2), y1+y2


def BinToCm(y, factor=10):
    """Computes the linear dimension for a given bin."""
    return math.exp(y / factor)

def CmToBin(x, factor=10):
    """Computes the bin for a given linear dimension."""
    return factor * math.log(x)


def Diameter(volume, factor=3/math.pi/4, exp=1/3.0):
    """Converts a volume to a diameter.

    d = 2r = 2 * (3/4/pi V)^1/3
    """
    return 2 * (factor * volume) ** exp


def Volume(diameter, factor=4*math.pi/3):
    """Converts a diameter to a volume.

    V = 4/3 pi (d/2)^3
    """
    return factor * (diameter/2.0)**3


class Cache(object):
    def __init__(self):
        """cache: maps from size bin to a list of sequences that could be
           observed in that bin.

           init_rdt: sequence of (V0, rdt) pairs
        """
        self.cache = {}
        self.initial_rdt = []

    def GetKeys(self):
        return self.cache.iterkeys()

    def GetBin(self, bin):
        return self.cache[bin]

    def Add(self, rdt, initial, final, seq):
        """Adds a sequence to the bin the corresponds to final."""
        cm = Diameter(final)
        bin = round(CmToBin(cm))
        self.cache.setdefault(bin, []).append(seq)
        self.initial_rdt.append((initial, rdt))

    def Print(self):
        """Prints the size (cm) for each bin, and the number of sequences."""
        for bin in sorted(cache.GetKeys()):
            ss = cache.GetBin(bin)
            diameter = BinToCm(bin)
            print diameter, len(ss)
        
    def Correlation(self):
        vs, rdts = zip(*self.initial_rdt)
        lvs = [math.log(v) for v in vs]
        return correlation.Corr(vs, rdts)



cache = Cache()

def ExtendSequence(t, rdt, interval):
    """Generates a new random value and adds it to the end of t.

    Side-effect: adds sub-sequences to the cache.

    t: sequence of values so far
    rdt: reciprocal doubling time in doublings per year
    interval: timestep in years
    """
    initial = t[-1]
    doublings = rdt * interval
    final = initial * 2**doublings
    res = t + (final,)
    cache.Add(rdt, initial, final, res)
    
    return final, res


def MakeSequence(iterator, v0=0.01, interval=INTERVAL, vmax=Volume(20.0)):
    """Simulate the growth of a tumor.

    n: number of time steps
    v0: initial volume in mL (cm^3)
    interval: timestep in years
    """
    vs = v0,

    for rdt in iterator:
        final, vs = ExtendSequence(vs, rdt, interval)
        if final > vmax:
            break

    return vs


def MakeSequences(n, rho, cdf):
    """Returns a sequence of times and a list of sequences of volumes."""
    sequences = []
    for i in range(n):
        iterator = RdtGenerator(n, rho, cdf)
        vs = MakeSequence(iterator)
        sequences.append(vs)

        if i % 100 == 0:
            print i

    return sequences


def PlotSequence(ts, vs, color='blue'):
    """Plots a time series of linear measurements.

    ts: sequence of times in years
    vs: sequence of columes
    color: color string
    """
    line_options = dict(color=color, linewidth=1, alpha=0.2)
    xs = [Diameter(v) for v in vs]
    myplot.Plot(ts, xs,
                line_options=line_options,
                yscale='log',
                clf=False)

def PlotSequences(ss):
    """Plots linear measurement vs time.

    ts: sequence of times
    ss: list of sequences of volumes
    """
    pyplot.clf()
    line_options = dict(color='gray', linewidth=1, linestyle='dashed')
    myplot.Plot([0, 40], [10, 10], line_options=line_options)
    for vs in ss:
        n = len(vs)
        age = n * INTERVAL
        ts = numpy.linspace(0, age, n)
        PlotSequence(ts, vs)

    myplot.Save(root='kidney4',
                formats=FORMATS,
                axis=[0, 40, MINSIZE, 20],
                title='Simulations of tumor growth',
                xlabel='tumor age (years)',
                yticks=MakeTicks([0.2, 0.5, 1, 2, 5, 10, 20]),
                ylabel='diameter (cm, log scale)')


def PlotBin(bin, color='blue'):
    "Plots the set of sequences for the given bin."""
    ss = cache.GetBin(bin)
    for vs in ss:
        n = len(vs)
        age = n * INTERVAL
        ts = numpy.linspace(-age, 0, n)
        PlotSequence(ts, vs, color)


def PlotCache():
    """Plots the set of sequences for each bin."""
    # 9.97 cm
    bins = [23.0]
    colors = ['blue', 'green', 'red', 'cyan']
    cdfs = []

    pyplot.clf()
    for bin, color in zip(bins, colors):
        PlotBin(bin, color)

    myplot.Save(root='kidney5',
                formats=FORMATS,
                title='History of simulated tumors',
                axis=[-40, 1, MINSIZE, 20],
                xlabel='years',
                ylabel='diameter (cm, log scale)')


def CdfBin(bin, name=''):
    """Forms the cdf of ages for the sequences in this bin."""
    ss = cache.GetBin(bin)
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

    myplot.Cdfs(cdfs,
                root='kidney6',
                title='Distribution of age for several diameters',
                formats=FORMATS,
                xlabel='tumor age (years)',
                ylabel='CDF',
                loc=4)


def PrintCI(fp, cm, ps):
    """Writes a line in the LaTeX table.

    fp: file pointer
    cm: diameter in cm
    ts: tuples of percentiles
    """
    fp.write('%0.1f' % round(cm, 1))
    for p in reversed(ps):
        fp.write(' & %0.1f ' % round(p, 1))
    fp.write(r'\\' '\n')


def PrintTable(fp, xs, ts):
    """Writes the data in a LaTeX table.

    fp: file pointer
    xs: diameters in cm
    ts: sequence of tuples of percentiles
    """
    fp.write(r'\begin{tabular}{|r||r|r|r|r|r|}' '\n')
    fp.write(r'\hline' '\n')
    fp.write(r'Diameter   & \multicolumn{5}{c|}{Percentiles of age} \\' '\n')
    fp.write(r'(cm)   & 5th & 25th & 50th & 75th & 95th \\' '\n')
    fp.write(r'\hline' '\n')

    for i, (cm, ps) in enumerate(zip(xs, ts)):
        #print cm, ps
        if i % 3 == 0:
            PrintCI(fp, cm, ps)

    fp.write(r'\hline' '\n')
    fp.write(r'\end{tabular}' '\n')


def FitLine(xs, ys, fxs):
    """Fits a line to the xs and ys, and returns fitted values for fxs.

    Applies a log transform to the xs.

    xs: diameter in cm
    ys: age in years
    fxs: diameter in cm
    """
    lxs = [math.log(x) for x in xs]
    inter, slope = correlation.LeastSquares(lxs, ys)
    res = correlation.Residuals(lxs, ys, inter, slope)
    R2 = correlation.CoefDetermination(ys, res)

    lfxs = [math.log(x) for x in fxs]
    fys = [inter + slope * x for x in lfxs]
    return fys


def ConfidenceIntervalFigure(xscale='linear'):
    """Plots the confidence interval for each bin."""
    xs = []
    ts = []
    percentiles = [95, 75, 50, 25, 5]
    min_size = 0.3
    
    # loop through the bins, accumulate
    # xs: sequence of sizes in cm
    # ts: sequence of percentile tuples
    for i, bin in enumerate(sorted(cache.GetKeys())):
        cm = BinToCm(bin)
        if cm < min_size or cm > 20.0:
            continue
        xs.append(cm)
        cdf = CdfBin(bin)      
        ps = [cdf.Percentile(p) for p in percentiles]
        ts.append(ps)

    # dump the results into a table
    fp = open('kidney_table.tex', 'w')
    PrintTable(fp, xs, ts)
    fp.close()

    # make the figure
    linewidths = [1, 2, 3, 2, 1]
    alphas = [0.3, 0.5, 1, 0.5, 0.3]
    labels = ['95th', '75th', '50th', '25th', '5th']

    # transpose the ts so we have sequences for each percentile rank
    pyplot.clf()
    yys = zip(*ts)

    for ys, linewidth, alpha, label in zip(yys, linewidths, alphas, labels):
        line_options = dict(color='blue', linewidth=linewidth, 
                            alpha=alpha, label=label, markersize=2)

        # plot the data points
        myplot.Plot(xs, ys, 'bo',
                    line_options=line_options,
                    clf=False,
                    legend=False)

        # plot the fit lines
        fxs = [min_size, 20.0]
        fys = FitLine(xs, ys, fxs)

        myplot.Plot(fxs, fys,
                    clf=False,
                    line_options=line_options,
                    legend=False)

        # put a label at the end of each line
        x, y = fxs[-1], fys[-1]
        pyplot.text(x*1.05, y, label, 
                    horizontalalignment='left',
                    verticalalignment='center')

    # make the figure
    myplot.Save(root='kidney7',
                formats=FORMATS,
                title='Confidence interval for age vs diameter',
                xlabel='diameter (cm, log scale)',
                ylabel='tumor age (years)',
                xscale=xscale,
                xticks=MakeTicks([0.5, 1, 2, 5, 10, 20]),
                axis=[0.25, 35, 0, 45],
                legend=False)


def MakeTicks(xs):
    """Makes a pair of sequences for use as pyplot ticks.

    xs: sequence of floats

    Returns (xs, labels), where labels is a sequence of strings.
    """
    labels = [str(x) for x in xs]
    return xs, labels


def TestCorrelation(cdf):
    """Tests the correlated generator.

    Makes sure that the sequence has the right distribution and correlation.
    """
    n = 10000
    rho = 0.4

    iterator = CorrelatedGenerator(cdf, rho)
    xs = [iterator.next() for i in range(n)]
    
    rho2 = correlation.SerialCorr(xs)
    print rho, rho2
    cdf2 = Cdf.MakeCdfFromList(xs)

    myplot.Cdfs([cdf, cdf2], show=True)


def main(script):
    random.seed(17)

    cdf = MakeCdf()

    lam1 = FitCdf(cdf)
    fit = GenerateCdf(lam1=lam1)

    #TestCorrelation(fit)

    PlotCdf(cdf, fit)
    QQPlot(cdf, fit)

    rho = 0.0

    ss = MakeSequences(100, rho, fit)

    PlotSequences(ss)
    PlotCache()

    ss = MakeSequences(1900, rho, fit)
    print 'V0-RDT correlation', cache.Correlation()

    CdfCache()

    ConfidenceIntervalFigure(xscale='log')
    #cache.Print()


if __name__ == '__main__':
    profile = False
    if profile:
        import cProfile
        cProfile.run('main(*sys.argv)')
    else:
        main(*sys.argv)


