"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib.pyplot as pyplot
import myplot
import csv

import correlation
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

    Regress(model_number=0)
    Regress(model_number=1)


def PlotProbs(filename='p.heri.31'):
    pyplot.clf()
    for x in [1975.5, 1984.5, 1998.5, 2006.5]:
        xs = [x, x]
        ys = [0, 1]
        pyplot.plot(xs, ys, color='0.8', linewidth=10)

    data = ReadData(filename)
    xs, ys = zip(*data)
    pyplot.plot(xs, ys, 'bo-', color='blue', linewidth=2, markersize=10)
    myplot.Save(root='heri2',
                title='Location of changepoints',
                xlabel='',
                ylabel='cumulative probability',
                axis=[1972, 2010, 0, 1])


def PlotLinearFit(coeffs, ts):
    inter = GetEst(coeffs[0])
    slope = GetEst(coeffs[1])
    ys = [inter + slope * t for t in ts]
    pyplot.plot(fts, fys, color='gray', linewidth=2)


def PlotRandomLinearFit(coeffs, ts, n=10):
    for i in range(n):
        inter = RandomEst(coeffs[0])
        slope = RandomEst(coeffs[1])
        print inter, slope
        ys = [inter + slope * t for t in ts]
        pyplot.plot(ts, ys, color='gray', linewidth=2)


def RandomEst(coeff):
    name, est, stderr = coeff
    print est, stderr
    return random.gauss(est, stderr/100)


def GetEst(coeff):
    name, est, stderr = coeff
    return est


def EvalLinearModel(coeffs, ts):
    inter = GetEst(coeffs[0])
    slope = GetEst(coeffs[1])
    print inter, slope
    print type(slope)
    ys = [inter + slope * t for t in ts]
    return ts, ys

    

def RunModel(model, print_flag=True):
    """Submits model to r.lm and returns the result."""
    model = r(model)
    res = r.lm(model)
    if print_flag:
        PrintSummary(res)
    return res


def Regress(filename='heri.1', model_number=0):
    data = ReadData(filename)
    ts, ys = zip(*data)
    t2 = [t**2 for t in ts]

    # put the data into the R environment
    robjects.globalEnv['ts'] = robjects.FloatVector(ts)
    robjects.globalEnv['t2'] = robjects.FloatVector(t2)
    robjects.globalEnv['ys'] = robjects.FloatVector(ys)

    # run the models
    models = ['ys ~ ts',
              'ys ~ ts + t2']
    model = models[model_number]
    res = RunModel(model)

    coeffs = GetCoefficients(res)
    PlotRandomLinearFit(coeffs, ts)


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

    PlotDiffs(root='heri3', flag=True)


if __name__ == '__main__':
    import sys
    main(*sys.argv)

