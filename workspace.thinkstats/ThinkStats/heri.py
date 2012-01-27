"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib.pyplot as pyplot
import myplot
import csv

import Cdf
import correlation
import math
import random
import thinkstats

import rpy2.robjects as robjects
r = robjects.r


def ReadData(filename):
    """Reads a CSV file of SAT scores.

    Args:
      filename: string filename

    Returns:
      list of (score, number) pairs
    """
    fp = open(filename)
    res = []

    for line in fp:
        try:
            t = [float(x) for x in line.split()]
            res.append(t)
        except ValueError:
            pass

    return res


def Pvalue(filename='heri.0', delta=0.033, n=100000):
    data = ReadData(filename)
    
    count = 0
    for i in range(n):
        xs, ys = zip(*data)
        ys = list(ys)
        random.shuffle(ys)
        inter, slope = correlation.LeastSquares(xs, ys)

        if abs(slope) > delta:
            count += 1
    return float(count) / n


def Fit(xs, ys):
    """Find the linear least squares fit between xs and ys."""
    inter, slope = correlation.LeastSquares(xs, ys)
    print '(inter, slope):', inter, slope

    res = correlation.Residuals(xs, ys, inter, slope)
    R2 = correlation.CoefDetermination(ys, res)

    print 'inter', inter
    print 'slope', slope
    print 'R^2', R2
    print

    return inter, slope, R2


def PlotDiffs(filename='heri.0', root='heri1', flag=False):
    pyplot.clf()

    data = ReadData(filename)
    xs, ys = zip(*data)

    if flag:
        RunFit(xs, ys)
    
    pyplot.plot(xs, ys, 'b.:', markersize=15)
    myplot.Save(root=root,
                title='Yearly changes',
                xlabel='',
                ylabel='percentage points',
                axis=[1972, 2013, -1.2, 2.1])


def RunFit(xs, ys):
    inter, slope, R2 = Fit(xs, ys)
    fxs = [min(xs), max(xs)]
    fys = [inter + slope*x for x in fxs]
    pyplot.plot(fxs, fys, 'r-', linewidth=2)
    print 'Mean diff', thinkstats.Mean(ys)
    print 'Current rate:', fys[-1]


def PlotProbs(filename='p.heri.31'):
    pyplot.clf()
    for x in [1975.5, 1984.5, 1998.5, 2006.5]:
        xs = [x, x]
        ys = [0, 1]
        pyplot.plot(xs, ys, color='0.8', linewidth=10)

    data = ReadData(filename)
    xs, ys = zip(*data)
    pyplot.plot(xs, ys, 'bo-', color='blue', linewidth=2, markersize=6)
    myplot.Save(root='heri2',
                title='Location of changepoints',
                xlabel='',
                ylabel='cumulative probability',
                axis=[1972, 2010, 0, 1])


def RunModel(model, print_flag=False):
    """Submits model to r.lm and returns the result."""
    model = r(model)
    res = r.lm(model)
    if print_flag:
        PrintSummary(res)
    return res


def Regress(model, ys, ts, print_flag=False, drop=1):
    t2 = [t**2 for t in ts]

    # put the data into the R environment
    robjects.globalEnv['ts'] = robjects.FloatVector(ts[:-drop])
    robjects.globalEnv['t2'] = robjects.FloatVector(t2[:-drop])
    robjects.globalEnv['ys'] = robjects.FloatVector(ys[:-drop])

    res = RunModel(model, print_flag)
    coeffs = GetCoefficients(res)
    return coeffs


def EvalFit(coeffs, ts):
    betas = [GetEst(coeff) for coeff in reversed(coeffs)]
    fys = [Horner(betas, t) for t in ts]
    return fys


def MakeFit(model, ys, ts):
    coeffs = Regress(model, ys, ts)
    fys = EvalFit(coeffs, ts)
    return fys


def Residuals(model, ys, ts):
    coeffs = Regress(model, ys, ts, print_flag=True)
    fys = EvalFit(coeffs, ts)
    residuals = [fy - y for fy, y in zip(fys, ys)]
    return residuals


def Horner(betas, t):
    total = 0
    for beta in betas:
        total = total * t + beta
    return total


def MakeSampleError(model, ys, ts, n=100):
    # make a model of the residuals
    residuals = Residuals(model, ys, ts)
    mu, var = thinkstats.MeanVar(residuals)
    sig = math.sqrt(var)

    # make the best fit
    fys = MakeFit(model, ys, ts)

    # resample residuals and generate fake fits
    fits = []
    for i in range(n):
        fake_ys = [fy + random.gauss(mu, sig) for fy in fys]
        fake_fys = MakeFit(model, fake_ys, ts)
        fits.append(fake_fys)

    # plot the 90% CI in each column
    data = zip(*fits)
    cdfs = [Cdf.MakeCdfFromList(ys) for ys in data]
    max_fys = [cdf.Percentile(95) for cdf in cdfs]
    min_fys = [cdf.Percentile(5) for cdf in cdfs]

    return min_fys, max_fys


def PlotData(ys, ts):
    pyplot.plot(ts, ys, 'bo-', linewidth=2, markersize=10)


def GetData(filename):
    data = ReadData(filename)
    ts, ys = zip(*data)
    return ts, ys


def MakePlot(filename='heri.1'):
    # do the analysis without the last data point
    ts, ys = GetData(filename)
    shift = ts[0]
    tshift = [t-shift for t in ts]

    # plot sample error cone
    model = 'ys ~ ts + t2'
    min_fys, max_fys = MakeSampleError(model, ys, tshift)
    pyplot.fill_between(ts, min_fys, max_fys, color='0.9')

    # plot fit
    fys = MakeFit(model, ys, tshift)
    pyplot.plot(ts, fys, color='red', linewidth=2)

    # plot data (including the last data point
    PlotData(ys, ts)

    myplot.Save(root='heri4',
                title='',
                xlabel='',
                ylabel='Percent None',
                axis=[1968, 2013, 0, 27])


def PlotResiduals(ts, ys):
    residuals = Residuals(ts, ys)
    pyplot.clf()
    pyplot.plot(ts, residuals, 'bo-', linewidth=2)
    myplot.Save(root='heri5',
                title='',
                xlabel='',
                ylabel='Residuals',
                axis=[1968, 2012, -6, 6])



def GetEst(coeff):
    name, est, stderr = coeff
    return est


def PrintSummary(res):
    """Prints results from r.lm (just the parts we want)."""
    flag = False
    lines = r.summary(res)
    lines = str(lines)

    for line in lines.split('\n'):
        # skip everything until we get to coefficients
        if line.startswith('Coefficients'):
            flag = True
        if flag:
            print line
    print


def GetCoefficients(res):
    """Prints results from r.lm (just the parts we want)."""
    flag = False
    lines = r.summary(res)
    lines = str(lines)

    coeffs = []
    for line in lines.split('\n'):
        line = line.strip()
        if flag:
            t = line.split()
            if len(t) < 5:
                break
            name, est, stderr = t[0], float(t[1]), float(t[2])
            coeffs.append((name, est, stderr))

        # skip everything until we get to the coefficients
        if line.startswith('Estimate'):
            flag = True

    return coeffs


def main(script):

    MakePlot()


if __name__ == '__main__':
    import sys
    main(*sys.argv)

