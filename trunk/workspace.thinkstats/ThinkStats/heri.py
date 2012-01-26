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

    Regress(xs, ys)

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


def RunModel(model, print_flag=True):
    """Submits model to r.lm and returns the result."""
    model = r(model)
    res = r.lm(model)
    if print_flag:
        PrintSummary(res)
    return res


def Regress(xs, ys, model_number=0):

    # put the data into the R environment
    robjects.globalEnv['xs'] = robjects.FloatVector(xs)
    robjects.globalEnv['ys'] = robjects.FloatVector(ys)

    # run the models
    models = ['ys ~ xs']
    model = models[model_number]
    print model
    RunModel(model)


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


def main(script):

    PlotDiffs(root='heri3', flag=True)


if __name__ == '__main__':
    import sys
    main(*sys.argv)

