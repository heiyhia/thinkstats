"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import numpy
import cPickle
import random

import brfss
import erf
import myplot
import Pmf
import thinkstats

import matplotlib.pyplot as pyplot

spread = 4
digits = 2
normal_pmf = erf.FixedPointNormalPmf(spread, digits, log=True)


def BatchUpdate(suite, t, batch_size=20):
    """Updates a suite of hypotheses based on new evidence.

    Modifies the suite directly; if you want to keep the original, make
    a copy.

    Args:
        suite: Pmf object
        evidence: whatever kind of object Likelihood expects
        batch_size: number of element to process before normalizing
    """
    for i in range(0, len(t), batch_size):
        print i, i+batch_size
        batch = t[i:i+batch_size]
        Update(suite, batch)


def LogUpdate(suite, evidence):
    """Updates a suite of hypotheses based on new evidence.

    Modifies the suite directly; if you want to keep the original, make
    a copy.

    Args:
        suite: Pmf object
        evidence: whatever kind of object Likelihood expects
    """
    for hypo in suite.Values():
        likelihood = LogLikelihood(evidence, hypo)
        suite.Incr(hypo, likelihood)
    print suite.Total()


def LogLikelihood(evidence, hypo):
    """Computes the likelihood of the evidence assuming the hypothesis is true.

    Args:
        evidence: a tuple of (number of heads, number of tails)
        hypo: float probability of heads

    Returns:
        probability of tossing the given number of heads and tails with a
        coin that has p probability of heads
    """
    t = evidence
    mu, sigma = hypo

    total = Summation(t, mu)
    return -len(t) * math.log(sigma) - total / 2 / sigma**2


def Summation(t, mu, cache={}):
    try:
        return cache[t, mu]
    except KeyError:
        total = sum((x-mu)**2 for x in t)
        cache[t, mu] = total
        return total

def MakeUniformPrior(t, num_points, spread=3.0):
    """Makes a prior distribution for mu and sigma based on a sample.

    t: sample
    num_points: number of values in each dimension
    spread: number of standard errors to include

    Returns: Pmf that maps from (mu, sigma) to prob.
    """
    # estimate mean and stddev of t
    n = len(t)
    xbar, S2 = thinkstats.MeanVar(t)
    sighat = math.sqrt(S2)

    print xbar, sighat

    # compute standard error for mu and the range of ms
    stderr_xbar = sighat / math.sqrt(n)
    mspread = spread * stderr_xbar
    ms = numpy.linspace(xbar-mspread, xbar+mspread, num_points)

    # compute standard error for sigma and the range of ss
    stderr_sighat = sighat / math.sqrt(2 * (n-1))
    sspread = spread * stderr_sighat
    ss = numpy.linspace(sighat-sspread, sighat+sspread, num_points)

    # populate the PMF
    pmf = Pmf.Pmf()
    for m in ms:
        for s in ss:
            pmf.Set((m, s), 1)
    return ms, ss, pmf


def PlotPosterior(xs, ys, suite, pcolor=True, contour=False):
    """Makes a contour plot.
    
    xs: sequence of values
    ys: sequence of values
    suite: Pmf that maps (x, y) to z
    """
    X, Y = numpy.meshgrid(xs, ys)
    func = lambda x, y: suite.Prob((x, y))
    prob = numpy.vectorize(func)
    Z = prob(X, Y)

    if pcolor:
        pyplot.pcolor(X, Y, Z)
    if contour:
        pyplot.contour(X, Y, Z)
    pyplot.show()


def DumpHeights(data_dir='.', n=10000):
    resp = brfss.Respondents()
    resp.ReadRecords(data_dir, n)

    d = {1:[], 2:[]}
    [d[r.sex].append(r.htm3) for r in resp.records if r.htm3 != 'NA']

    fp = open('bayes_height_data.pkl', 'wb')
    cPickle.dump(d, fp)
    fp.close()


def LoadHeights():
    fp = open('bayes_height_data.pkl', 'r')
    d = cPickle.load(fp)
    fp.close()
    return d


def EstimateParameters(t, num_points=21):
    xs, ys, suite = MakeUniformPrior(t, num_points)
    suite.Log()

    LogUpdate(suite, tuple(t))

    suite.Exp()
    suite.Normalize()

    PlotPosterior(xs, ys, suite)


def main():
    if False:
        random.seed(16)
        t = [random.gauss(3, 5) for i in range(100000)]
        EstimateParameters(t)
        return

    DumpHeights(n=10000)
    d = LoadHeights()

    for key, t in d.iteritems():
        print key, len(t)
        EstimateParameters(t)
        break


if __name__ == '__main__':
    main()
