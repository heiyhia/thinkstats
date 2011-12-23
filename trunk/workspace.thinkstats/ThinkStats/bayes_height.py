"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import numpy
import random

import erf
import myplot
import Pmf
import thinkstats

import matplotlib.pyplot as pyplot

spread = 4
digits = 2
normal_pmf = erf.FixedPointNormalPmf(spread, digits)


def MakeUniformSuite(low, high, steps):
    """Makes a PMF that represents a suite of hypotheses with equal p.
    
    Args:
        low: low end of range
        high: high end of range
        steps: number of values

    Returns:
        Pmf object
    """
    hypos = [low + (high-low) * i / (steps-1.0) for i in range(steps)]
    pmf = Pmf.MakePmfFromList(hypos)
    return pmf


def Update(suite, evidence):
    """Updates a suite of hypotheses based on new evidence.

    Modifies the suite directly; if you want to keep the original, make
    a copy.

    Args:
        suite: Pmf object
        evidence: whatever kind of object Likelihood expects
    """
    for hypo in suite.Values():
        likelihood = Likelihood(evidence, hypo)
        suite.Mult(hypo, likelihood)
    suite.Normalize()


def Likelihood(evidence, hypo):
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
    return NormalLikelihood(t, mu, sigma)


def NormalLikelihood(t, mu, sigma):
    prod = 1.0
    for x in t:
        prod *= normal_pmf.NormalProb((x-mu) / sigma)
    return prod


def MakeUniformPrior(t, num_points):
    n = len(t)
    xbar, S2 = thinkstats.MeanVar(t)
    sighat = math.sqrt(S2)

    stderr_xbar = sighat / math.sqrt(n)
    spread = 2
    mspread = spread * stderr_xbar
    ms = numpy.linspace(xbar-mspread, xbar+mspread, num_points)

    stderr_sighat = sighat / math.sqrt(2 * (n-1))
    sspread = spread * stderr_sighat
    
    ss = numpy.linspace(sighat-sspread, sighat+sspread, num_points)

    pmf = Pmf.Pmf()
    for m in ms:
        for s in ss:
            pmf.Set((m, s), 1)
    return ms, ss, pmf


def ContourPlot(xs, ys, suite):
    X, Y = numpy.meshgrid(xs, ys)
    func = lambda x, y: suite.Prob((x, y))
    prob = numpy.vectorize(func)
    Z = prob(X, Y)

    pyplot.contour(X, Y, Z)
    pyplot.show()


def main():
    random.seed(17)
    t = [random.gauss(3, 4) for i in range(20)]

    xs, ys, suite = MakeUniformPrior(t, 11)

    Update(suite, t)
    suite.name = 'posterior'

    for hypo, p in suite.Items():
        print hypo, p

    print suite.Total()

    ContourPlot(xs, ys, suite)

if __name__ == '__main__':
    main()
