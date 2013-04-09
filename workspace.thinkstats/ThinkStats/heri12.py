"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
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

FORMATS = ['pdf', 'png']

def Regress(model, ys, ts, print_flag=False):
    """Run a linear regression using rpy2.

    Returns a list of coefficients (which are rpy2.RVectors)
    Use GetEst to extract the estimated coefficient from a coeff.
    """
    t2 = [t**2 for t in ts]

    # put the data into the R environment
    robjects.globalenv['ts'] = robjects.FloatVector(ts)
    robjects.globalenv['t2'] = robjects.FloatVector(t2)
    robjects.globalenv['ys'] = robjects.FloatVector(ys)

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


def MakeFit(model, ys, ts, extrap=1):
    """Fit a model to the data and return the fitted values."""
    coeffs = Regress(model, ys, ts)

    fts = ts + [ts[-1] + 1]
    fys = EvalFit(coeffs, fts)
    return fts, fys


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
    fts, fys = MakeFit(model, ys, ts)

    # resample residuals and generate hypothetical fits
    fits = []
    for i in range(n):
        fake_ys = [fy + random.gauss(mu, sig) for fy in fys[:-1]]
        _, fake_fys = MakeFit(model, fake_ys, ts)
        fits.append(fake_fys)

    # find the 90% CI in each column
    columns = zip(*fits)

    sample_error = MakeStddev(columns)
    total_error = MakeStddev(columns, mu, var)

    return fts, sample_error, total_error


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
    """Reads a CSV file of data from HERI scores.

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

    return zip(*res)


def MakePlot(ts, ys, model):
    """Generates a plot with the data, a fitted model, and error bars."""
    pyplot.clf()

    # shift the times to start at 0 (but use the originals for plots)
    shift = ts[0]
    tshift = [t-shift for t in ts]

    # plot the error models
    fts, sample_error, total_error = MakeErrorModel(model, ys, tshift)
    fts = [t+shift for t in fts]

    pyplot.fill_between(fts, *total_error, color='0.8', alpha=0.5, linewidth=0)
    pyplot.fill_between(fts, *sample_error, color='0.6', alpha=0.5, linewidth=0)

    # plot the estimated fit
    fts, fys = MakeFit(model, ys, tshift)
    fts = [t+shift for t in fts]

    pyplot.plot(fts, fys, color='red', linewidth=2, alpha=0.5)


def PlotResiduals(ts, ys):
    residuals = Residuals(ts, ys)
    pyplot.clf()
    pyplot.plot(ts, residuals, 'bo-', linewidth=2)
    myplot.Save(root='heri5',
                formats=FORMATS,
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


def PlotReligious(filename):
    """Make plots showing percentage None broken down by college type.

    filename: string filename
    """
    fp = open(filename)
    reader = csv.reader(fp)

    header = reader.next()

    years = []
    rows = []
    for t in reader:
        year = int(t[0])
        data = t[1:6]
        try:
            data = [float(x) for x in data]
        except ValueError:
            continue
        years.append(year)
        rows.append(data)

    cols = zip(*rows)
    labels = header[1:]

    PlotReligiousSubset(years, cols, labels, 0, 3)
    PlotReligiousSubset(years, cols, labels, 3, 5)


def PlotReligiousSubset(years, cols, labels, i, j):
    """Helper function that factors out common plotting code.

    years: sequence of years
    cols: list of columns to plot
    labels: list of labels (corresponding to cols)
    i,j: slice indices of the columns to plot
    """
    pyplot.clf()
    options = dict(linewidth=3, markersize=0, alpha=0.7)
    for col, label in zip(cols[i:j], labels[i:j]):
        pyplot.plot(years, col, label=label, **options)

    root = 'heri.religious.%d.%d' % (i, j)
    myplot.Save(root=root,
                formats=FORMATS,
                xlabel='Year',
                ylabel='% None',
                title='Religious preference')


def PlotReligious2(filename):
    """Makes a plot of the number/fraction of students at religious colleges.

    filename: string filename
    """
    fp = open(filename)
    reader = csv.reader(fp)

    header = reader.next()

    years = [int(x) for x in header[1:]]

    labels = []
    cols = []
    for t in reader:
        label = t[0]
        if label == '':
            break
        col = [float(x)/1000000.0 for x in t[1:]]
        labels.append(label)
        cols.append(col)

    PlotReligiousScale2(years, cols, labels, flag='raw')

    cols = PercentTotal(cols)
    PlotReligiousScale2(years, cols, labels, flag='percent')


def PercentTotal(cols):
    """Converts the data in cols to percentages of total.

    Modifies the columns.

    cols: sequence of columns
    """
    totals = []
    for row in zip(*cols):
        print row
        total = sum(row)
        totals.append(total)

    for col in cols:
        for i in range(len(col)):
            col[i] /= totals[i] / 100
    
    return cols


def PlotReligiousScale2(years, cols, labels, flag):
    """Helper function that factors out common plotting code.

    years: sequence of years
    cols: list of columns to plot
    labels: list of labels (corresponding to cols)
    flag: string 'raw' or 'percent'
    """
    pyplot.clf()
    options = dict(linewidth=3, markersize=0, alpha=0.7)
    for col, label in zip(cols, labels):
        pyplot.plot(years, col, label=label, **options)

    root = 'heri.religious2.%s' % flag
    ylabel = dict(raw='Enrollment (millions)',
                  percent='Enrollment (percent of total)')[flag]
    axis = dict(raw=[1977, 2010, 0, 16],
                  percent=[1977, 2010, 0, 100])[flag]

    myplot.Save(root=root,
                formats=FORMATS,
                xlabel='Year',
                ylabel=ylabel,
                title='Enrollment by college type',
                axis=axis)


def main(script):

    # make plots showing enrollment at private, religious and non-religious
    # colleges
    #PlotReligious2('heri.religious2.csv')
    #return
    upper = 2014

    ts, ys = ReadData('heri.0')
    MakePlot(ts, ys, model='ys ~ ts')

    options = dict(linewidth=3, markersize=0, alpha=0.7)
    pyplot.plot(ts, ys, color='purple', label='Change in no religion',
                **options)

    myplot.Save(root='heri.0',
                formats=FORMATS,
                ylabel='Percentage points',
                loc=2,
                axis=[1967, upper, -3, 3])

    ts, ys = ReadData('heri.1')
    MakePlot(ts, ys, model='ys ~ ts + t2')

    pyplot.plot(ts, ys, 'bs-', label='No religion', **options)
    
    # add the actual value from 2012
    myplot.Plot([2012], [23.8], 'bs')

    myplot.Save(root='heri.1',
                formats=FORMATS,
                ylabel='Percent',
                loc=2,
                axis=[1967, upper, 0, 30])

    ts, ys = ReadData('heri.2')
    MakePlot(ts, ys, model='ys ~ ts + t2')

    pyplot.plot(ts, ys, 'go-', label='No attendance', **options)

    # add the actual value from 2012
    myplot.Plot([2012], [100 - 73.2], 'gs')

    myplot.Save(root='heri.2',
                formats=FORMATS,
                ylabel='Percent',
                loc=2,
                axis=[1967, upper, 0, 30])

    print (2011 - 1973) * 0.03548 - 0.36036 

if __name__ == '__main__':
    import sys
    main(*sys.argv)

