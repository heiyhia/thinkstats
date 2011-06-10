"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import myplot
import matplotlib.pyplot as pyplot
import Pmf
import math
import random

def log2(x, denom=math.log(2)):
    return math.log(x) / denom

def Entropy(suite):
    total = 0
    for x, p in suite.Items():
        if p:
            term = -p * log2(p)
            total += term
    return total

def EntropyOfP(suite):
    total = 0
    mu = suite.Mean()
    for p in [mu, 1-mu]:
        term = -p * log2(p)
        total += term
    return total

def RelativeEntropy(prior, posterior):
    total = 0
    for x, p in prior.Items():
        if p == 0:
            continue
        q = posterior.Prob(x)
        if q == 0:
            total = float('inf')
        else:
            term = -p * (log2(p) - log2(q))
            total += term
    return total

def CrossEntropy(prior, posterior):
    total = 0
    for x, p in prior.Items():
        q = posterior.Prob(x)
        if q == 0:
            total = float('inf')
        else:
            term = -p * log2(q)
            total += term
    return total

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

    nc = suite.Total()
    ic = -log2(nc)
    suite.Normalize()
    return nc, ic

def Likelihood(evidence, hypo):
    """Computes the likelihood of the evidence assuming the hypothesis is true.

    Args:
        evidence: a tuple of (number of blue, number of green)
        hypo: float probability of blue

    Returns:
        probability of tossing the given number of blue and green with a
        coin that has p probability of blue
    """
    blue, green = evidence
    p = hypo
    return math.pow(p, blue) * math.pow(1-p, green)


def MakeCurves(suite, sequence):
    suite = suite.Copy()

    d = dict(b=(1,0), g=(0,1))

    total = 0
    
    e0 = []
    e1 = []
    e2 = []

    for i, x in enumerate(sequence):
        entropy = Entropy(suite)
        entropyOfP = EntropyOfP(suite)

        evidence = d[x]
        nc, ic = Update(suite, evidence)
        new_entropy = Entropy(suite)
        total += ic

        print x, ic, 2 * log2(entropy/new_entropy)

        e0.append((i, ic))
        e1.append((i, entropy))
        e2.append((i, entropyOfP))

    return e0, e1, e2


def PlotCurves(curves, root):
    pyplot.clf()
    
    for i, curve in enumerate(curves):
        alpha = 1 if i==0 else 0.2
        color = 'red' if i==0 else 'blue'
        xs, ys = zip(*curve)
        pyplot.plot(xs, ys, color=color, alpha=alpha)

    myplot.Save(root=root,
                xlabel='# observations',
                ylabel='entropy (bits)')


def quick(suite, sequence):
    d = dict(b=(1,0), g=(0,1))

    for i, x in enumerate(sequence):
        nc, ic = Update(suite, d[x])

    suite.Print()

    suite.Mult(2.0/3, 0)
    suite.Normalize()

    suite.Print()
    
def main():
    suite = Pmf.MakePmfFromList([0, 1.0/3, 2.0/3, 1])
    quick(suite, 'b')
    return

    sequence = list('bbbbbbbbbbbbbbbbbggg')
    random.shuffle(sequence)
    e0, e1, e2 = MakeCurves(suite, sequence)

    #MakeEntropyPlots(suite)

def MakeEntropyPlots(suite):
    e1s = []
    e2s = []

    for i in range(20):
        e0, e1, e2 = MakeCurves(suite, sequence)
        e1s.append(e1)
        e2s.append(e2)
        random.shuffle(sequence)

    PlotCurves(e1s, root='entropy1')
    PlotCurves(e2s, root='entropy2')


if __name__ == '__main__':
    main()
