"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import datetime
import math
import random
import sys

import matplotlib.pyplot as pyplot

import bayes
import correlation
import Cdf
import myplot
import Pmf
import rankit
import thinkstats


def ReadScale(filename='sat_scale.csv', col=2):
    """Reads a CSV file of SAT scales (maps from raw score to standard score).

    Args:
      filename: string filename
      col: which column to start with (0=Reading, 2=Math, 4=Writing)

    Returns:
      list of (raw score, standardize score) pairs
    """
    def ParseRange(s):
        t = [int(x) for x in s.split('-')]
        return 1.0 * sum(t) / len(t)

    fp = open(filename)
    reader = csv.reader(fp)
    raws = []
    scores = []

    for t in reader:
        try:
            raw = int(t[col])
            raws.append(raw)
            score = ParseRange(t[col+1])
            scores.append(score)
        except:
            pass

    raws.sort()
    scores.sort()
    return thinkstats.Interpolator(raws, scores)


def ReadRanks(filename='sat_ranks.csv'):
    """Reads a CSV file of SAT scores.

    Args:
      filename: string filename

    Returns:
      list of (score, number) pairs
    """
    fp = open(filename)
    reader = csv.reader(fp)
    res = []

    for t in reader:
        try:
            score = int(t[0])
            number = int(t[1])
            res.append((score, number))
        except ValueError:
            pass

    return res


def Summarize(pmf):
    mu, var = pmf.Mean(), pmf.Var()
    sigma = math.sqrt(var)
    print 'mu, sigma', mu, sigma
    return mu, sigma


def ApplyLogistic(pmf, inter=-2.5, slope=10):
    mu, sigma = Summarize(pmf)
    new = Pmf.Pmf()
    for val, prob in sorted(pmf.Items()):
        z = inter + slope * StandardScore(val, mu, sigma)
        
        prob_admit = bayes.Logistic(z)
        new.Incr(val, prob * prob_admit)

        print val, z, prob_admit

    new.Normalize()
    mu, sigma = Summarize(new)
    return new


def StandardScore(val, mu, sigma):
    return (val-mu) / sigma


def SummarizeR(r, sigma):
    """Prints summary statistics about a correlation."""

    r2 = r**2
    var = sigma**2

    rmse_without = sigma
    rmse_with = math.sqrt(var * (1 - r2))

    reduction = 1.0 - (rmse_with / rmse_without)
    print 'r, R^2', r, r2
    print 'RMSE (reduction)', rmse_with, reduction


def ApplyLogit(pmf, denom):
    new = Pmf.Pmf()
    for val, prob in pmf.Items():
        if val > 0 and val < denom:
            x = bayes.Logit(val/denom)
            new.Incr(x, prob)
    return new


def ReverseScale(pmf, scale):
    """Applies the reverse scale to the values of a PMF.

    Args:
        pmf: Pmf object
        scale: Interpolator object

    Returns:
        new Pmf
    """
    new = Pmf.Pmf()
    for val, prob in pmf.Items():
        raw = scale.Reverse(val)
        new.Incr(raw, prob)
    return new


def SamplePmf(pmf, total, fraction=0.001):
    """Generates a sample from a PMF by making copies of each value."""
    t = []
    for val, prob in pmf.Items():
        n = int(prob * total * fraction)
        t.extend([val] * n)
    return t


def MakeNormalPlot(ys, root=None, lineoptions={}, **options):
    """Makes a normal probability plot.
    
    Args:
        ys: sequence of values
        lineoptions: dictionary of options for pyplot.plot        
        options: dictionary of options for myplot.Plot
    """
    n = len(ys)
    ys.sort()
    xs = [random.normalvariate(0.0, 1.0) for i in range(n)]
    xs.sort()

    inter, slope = correlation.LeastSquares(xs, ys)
    print 'inter, slope', inter, slope
    x_fit = [-4, 4]
    y_fit = [inter + slope*x for x in x_fit]

    pyplot.clf()
    pyplot.plot(x_fit, y_fit, 'r-')
    pyplot.plot(sorted(xs), sorted(ys), 'b.', markersize=2, **lineoptions)
 
    myplot.Plot(root,
                xlabel = 'Standard normal values',
                legend=False,
                **options)


def ShiftValues(pmf, shift):
    """Shifts the values in a PMF by shift.  Returns a new PMF."""
    new = Pmf.Pmf()
    for val, prob in pmf.Items():
        if val >= shift:
            x = val - shift
            new.Incr(x, prob)
    return new


def DivideValues(pmf, denom):
    """Divides the values in a PMF by denom.  Returns a new PMF."""
    new = Pmf.Pmf()
    for val, prob in pmf.Items():
        if val >= 0:
            x = 1.0 * val / denom
            new.Incr(x, prob)
    return new


def ProbBigger(pmf1, pmf2):
    """Returns the probability that a value from one pmf exceeds another."""
    total = 0.0
    for v1, p1 in pmf1.Items():
        for v2, p2 in pmf2.Items():
            if v1 > v2:
                total += p1 * p2
    return total


class Exam:
    """Encapsulates information about an exam.

    Contains the distribution of scaled scores and an
    Interpolator that maps between scaled and raw scores.
    """
    def __init__(self):
        self.scale = ReadScale()
        self.scores = ReadRanks()
        self.hist = Pmf.MakeHistFromDict(dict(self.scores))

        self.scaled = Pmf.MakePmfFromHist(self.hist)

        self.raw = ReverseScale(self.scaled, self.scale)
        self.shift = 0

        self.max_score = max(self.raw.Values())
        self.log = ApplyLogit(self.raw, denom=self.max_score)
        
        self.prior = DivideValues(self.raw, denom=self.max_score)

    def NormalPlot(self):
        """Makes a normal probability plot for the raw scores."""
        total = self.hist.Total()
        raw_sample = SamplePmf(self.raw, total, fraction=0.01)
        MakeNormalPlot(raw_sample, 
                       root='sat_normal',
                       ylabel='Raw scores (math)',)

    def Shift(self, shift):
        """Shifts the values in the raw distribution."""
        self.shift = 0
        self.raw = ShiftValues(self.raw, shift=shift)

    def Update(self, score):
        """Computes the posterior distribution based on score."""
        raw = self.scale.Reverse(score) - self.shift
        evidence = raw, self.max_score-raw

        updater = bayes.Binomial()
        posterior = updater.Posterior(self.prior, evidence)

        return posterior

    def ScaledCredibleInterval(self, score):
        """Computes the credible interval for someone with the given score."""
        posterior = self.Update(score)
        low, high = bayes.CredibleInterval(posterior, 90)
        scale_low = self.scale.Lookup(low * self.max_score)
        scale_high = self.scale.Lookup(high * self.max_score)
        return scale_low, scale_high
        
    def PlotPosteriors(self, low, high):
        """Computes posteriors for the given scores and plots them."""
        posterior_low = self.Update(low)
        posterior_high = self.Update(high)

        prob_bigger = ProbBigger(posterior_high, posterior_low)
        print "prob_bigger:", prob_bigger

        cdf1 = Cdf.MakeCdfFromPmf(posterior_low, 'posterior %d' % low)
        cdf2 = Cdf.MakeCdfFromPmf(posterior_high, 'posterior %d' % high)

        myplot.Cdfs([cdf1, cdf2],
                    xlabel='P', 
                    ylabel='CDF', 
                    axis=[0.6, 1.0, 0.0, 1.0],
                    root='sat_posteriors')


def main(script):

    SummarizeR(r=0.53, sigma=0.71)

    exam = Exam()

    exam.NormalPlot()

    low = 700
    high = low+70

    low_range = exam.ScaledCredibleInterval(low)
    print low, low_range, (low_range[1] - low_range[0]) / 2.0

    high_range = exam.ScaledCredibleInterval(high)
    print high, high_range, (high_range[1] - high_range[0]) / 2.0

    exam.PlotPosteriors(low, high)


def PlotScaledDist():
    cdf1 = Cdf.MakeCdfFromPmf(pmf, 'scaled')
    myplot.Cdfs([cdf1],
               xlabel='score', 
               ylabel='CDF', 
               show=False)


def PlotRawDist():
    cdf2 = Cdf.MakeCdfFromPmf(raw, 'raw')
    myplot.Cdfs([cdf2],
               xlabel='score', 
               ylabel='CDF', 
               show=False)


def InferLogit(pmf):
    admitted = ApplyLogistic(pmf)

    cdf1 = Cdf.MakeCdfFromPmf(pmf)
    cdf2 = Cdf.MakeCdfFromPmf(admitted)

    quartiles = cdf2.Percentile(25), cdf2.Percentile(50), cdf2.Percentile(75)
    print 'quartiles', quartiles

    myplot.Cdfs([cdf1, cdf2],
               xlabel='score', 
               ylabel='CDF', 
               show=True)
    
if __name__ == '__main__':
    main(*sys.argv)
