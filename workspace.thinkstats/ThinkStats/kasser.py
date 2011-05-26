"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import myplot
import Pmf
from math import pow


    
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
    pmf = Pmf.MakePmfFromList(hypos, name='prior')
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
        print hypo, likelihood
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
    heads, tails = evidence
    p = hypo
    return pow(p, heads) * pow(1-p, tails)

def TwoUp(suite):
    total = 0.0
    for p, prob in suite.Items():
        total += p**2 * prob
    return total

def main():
    prior = MakeUniformSuite(0.0, 1.0, 101)

    post1 = prior.Copy(name='post1')
    evidence = 1,1
    Update(post1, evidence)
    print 'P(two heads)', TwoUp(post1)

    post2 = prior.Copy(name='post2')
    evidence = 49, 49
    Update(post2, evidence)
    print 'P(two heads)', TwoUp(post2)

    # plot the posterior distributions
    myplot.Pmfs([prior, post1, post2], 
                root='kasser',
                xlabel='P(heads)',
                ylabel='Probability',
                show=True)

if __name__ == '__main__':
    main()
