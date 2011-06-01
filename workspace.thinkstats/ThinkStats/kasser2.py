"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import myplot
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

def main():
    suite = Pmf.MakePmfFromList([0, 1.0/3, 2.0/3, 1])
    print 'entropy', Entropy(suite)

    d = dict(b=(1,0), g=(0,1))

    total = 0
    sequence = list('bbbbbbbbbbbbbbbbbggg')
    random.shuffle(sequence)

    for x in sequence:
        prior = suite.Copy()
        evidence = d[x]
        nc, ic = Update(suite, evidence)
        total += ic
        entropy = Entropy(suite)
        relative = RelativeEntropy(prior, suite)
        print x, 'info, total, entropy', ic, total, entropy

    print 'total', total

    return

    suite.Print()

    # plot the posterior distributions
    myplot.Pmf(suite, 
               xlabel='P(urn i)',
               ylabel='Posterior probability',
               show=True)

if __name__ == '__main__':
    main()
