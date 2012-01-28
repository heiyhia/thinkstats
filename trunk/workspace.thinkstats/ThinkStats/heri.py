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


def Regress(model, ys, ts, print_flag=False, drop=1):
    """Run a linear regression using rpy2.

    Returns a list of coefficients (which are rpy2.RVectors)
    Use GetEst to extract the estimated coefficient from a coeff.
    """
    t2 = [t**2 for t in ts]

    # put the data into the R environment
    robjects.globalEnv['ts'] = robjects.FloatVector(ts[:-drop])
    robjects.globalEnv['t2'] = robjects.FloatVector(t2[:-drop])
    robjects.globalEnv['ys'] = robjects.FloatVector(ys[:-drop])

    model = r(model)
    res = r.lm(model)
    if print_flag:
        PrintSummary(res)

    coeffs = GetCoefficients(res)
    return coeffs


def GetEst(coeff):
    """Extracts the estimated coefficient from a coeff."""
    name, est, stderr = coeff
    return est


def GetCoefficients(res):
    """Extracts coefficients from r.lm.

    This is an awful function.  It actually generates a text representation
    of the results and then parses it.  Ack!

    Maybe the rpy2 interface (or it's documentation) will improve at
    some point so this nonsense is no longer necessary.
    """
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


def MakeFit(model, ys, ts):
    """Fit a model to the data and return the fitted values."""
    coeffs = Regress(model, ys, ts)
    fys = EvalFit(coeffs, ts)
    return fys


def Residuals(model, ys, ts):
    """Fit a model to the data and return the residuals."""
    coeffs = Regress(model, ys, ts, print_flag=True)
    fys = EvalFit(coeffs, ts)
    residuals = [fy - y for fy, y in zip(fys, ys)]
    return residuals


def EvalFit(coeffs, ts):
    """Evaluate a fitted model at a sequence of locations.

    coeffs: a list of coefficients as returned by rpy2
    ts: locations to evaluate the model

    Returns a list of fitted values.
    """
    betas = [GetEst(coeff) for coeff in reversed(coeffs)]
    fys = [Horner(betas, t) for t in ts]
    return fys


def Horner(betas, t):
    """Use Horner's method to evaluate a polynomial.

    betas: coefficients in decreasing order of power.
    t: where to evaluate
    """
    total = 0
    for beta in betas:
        total = total * t + beta
    return total


def MakeErrorModel(model, ys, ts, n=100):
    """Makes a model that captures sample error and residual error.

    model: string representation of the regression model
    ys:    dependent variable
    ts:    explanatory variable
    n:     number of simulations to run

    Returns a pair of models, where each model is a pair of rows.
    """
    # estimate mean and stddev of the residuals
    residuals = Residuals(model, ys, ts)
    mu, var = thinkstats.MeanVar(residuals)
    sig = math.sqrt(var)

    # make the best fit
    fys = MakeFit(model, ys, ts)

    # resample residuals and generate hypothetical fits
    fits = []
    for i in range(n):
        fake_ys = [fy + random.gauss(mu, sig) for fy in fys]
        fake_fys = MakeFit(model, fake_ys, ts)
        fits.append(fake_fys)

    # find the 90% CI in each column
    columns = zip(*fits)

    #sample_error = MakeIntervals(columns)
    #columns = AddResidualError(columns, mu, sig)
    #total_error = MakeIntervals(columns)

    sample_error = MakeStddev(columns)
    total_error = MakeStddev(columns, mu, var)

    return sample_error, total_error


def MakeStddev(columns, mu2=0, var2=0):
    """Finds a confidence interval for each column.

    Returns two rows: the low end of the intervals and the high ends.
    """
    stats = [thinkstats.MeanVar(ys) for ys in columns]
    
    min_fys = [mu1 + mu2 - 2 * math.sqrt(var1 + var2) for mu1, var1 in stats]
    max_fys = [mu1 + mu2 + 2 * math.sqrt(var1 + var2) for mu1, var1 in stats]
    return min_fys, max_fys


def MakeIntervals(columns, low=5, high=95):
    """Finds a confidence interval for each column.

    Returns two rows: the low end of the intervals and the high ends.
    """
    cdfs = [Cdf.MakeCdfFromList(ys) for ys in columns]
    min_fys = [cdf.Percentile(low) for cdf in cdfs]
    max_fys = [cdf.Percentile(high) for cdf in cdfs]
    return min_fys, max_fys


def AddResidualError(columns, mu, sig):
    """Adds Gaussian noise to the data in columns.

    columns: list of columns, where each column is a set of y-values
             for a given t-value
    mu, sig: parameters of the noise
    """
    return [[y + random.gauss(mu, sig) for y in col]
            for col in columns]


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


def MakePlot(filename, model, ylabel, axis):
    """Generates a plot with the data, a fitted model, and error bars."""
    pyplot.clf()

    # do the analysis without the last data point
    data = ReadData(filename)
    ts, ys = zip(*data)

    # shift the times to start at 0 (but use the originals for plots)
    shift = ts[0]
    tshift = [t-shift for t in ts]

    # plot the error models
    sample_error, total_error = MakeErrorModel(model, ys, tshift)
    pyplot.fill_between(ts, *total_error, color='0.9')
    pyplot.fill_between(ts, *sample_error, color='0.7')

    # plot the estimated fit
    fys = MakeFit(model, ys, tshift)
    pyplot.plot(ts, fys, color='red', linewidth=2)

    # plot the data
    pyplot.plot(ts, ys, 'bo-', linewidth=2, markersize=8)

    myplot.Save(root=filename,
                title='',
                xlabel='',
                ylabel=ylabel,
                axis=axis)


def PlotResiduals(ts, ys):
    residuals = Residuals(ts, ys)
    pyplot.clf()
    pyplot.plot(ts, residuals, 'bo-', linewidth=2)
    myplot.Save(root='heri5',
                title='',
                xlabel='',
                ylabel='Residuals',
                axis=[1968, 2012, -6, 6])


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
    MakePlot(filename='heri.0', 
             model='ys ~ ts', 
             ylabel='Change in Religion None',
             axis=[1967, 2013, -3, 3])
    MakePlot(filename='heri.1',
             model='ys ~ ts + t2',
             ylabel='Religion None',
             axis=[1967, 2013, 0, 28])
    MakePlot(filename='heri.2',
             model='ys ~ ts + t2',
             ylabel='Attendance None',
             axis=[1967, 2013, 0, 28])


if __name__ == '__main__':
    import sys
    main(*sys.argv)

