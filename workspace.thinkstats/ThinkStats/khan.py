"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import Pmf
import Cdf
import random
import myplot
import thinkstats

import matplotlib.pyplot as pyplot
import numpy
import scipy
import scipy.stats

def ChiSquared(expected, observed):
    """Compute the Chi-squared statistic for two tables.
    
    Args:
      expected: Hist of expected values
      observed: Hist of observed values
      
    Returns:
      float chi-squared statistic
    """
    total = 0.0
    for x, exp in expected.Items():
        obs = observed.Freq(x)
        total += (obs - exp)**2 / exp
    return total


def Simulate(pa, q, n):
    """Run a simulation.

    pa: probability of showing version A
    q:  probability of success for both A and B
    n:  number of trials
    
    Returns: Hist that maps (version, outcome) to frequency 
    """
    hist = Pmf.Hist()
    pairs = []
    for i in range(1, n+1):
        version = Flip(pa, 'A', 'B')
        outcome = Flip(q, 'Y', 'N')
        hist.Incr((version, outcome))

        
        expected = Expected(pa, q, i)
        chi2 = ChiSquared(expected, hist)
        pvalue = Pvalue(chi2, df=3)

        pairs.append((i, pvalue))

    return pairs
        
def Flip(p, y='Y', n='N'):
    return y if random.random() <= p else n

def Expected(pa, q, n):
    versions = Pmf.MakePmfFromDict(dict(A=pa, B=1-pa))
    outcomes = Pmf.MakePmfFromDict(dict(Y=q, N=1-q))

    hist = Pmf.Hist()
    for version, pp in versions.Items():
        for outcome, qq in outcomes.Items():
            hist.Incr((version, outcome), pp * qq * n)
    return hist

def Pvalue(chi2, df):
    return scipy.stats.chi2.cdf(chi2, df)

def main():
    MakeSpaghetti()

def Crosses(ps, thresh):
    for p in ps:
        if p <= thresh:
            return True
    return False

def MakeSpaghetti(iters=1000, lines=100, n=300, thresh=0.05):
    pyplot.plot([1,n], [thresh, thresh], color='red', alpha=1, linewidth=2)
    
    count = 0.0
    for i in range(iters):
        pairs = Simulate(0.5, 0.5, n)
        xs, ps = zip(*pairs)
        if Crosses(ps, thresh):
            count += 1

        if i < lines:
            pyplot.plot(xs, ps, alpha=0.2)
    
    print iters, count / iters

    myplot.Save(root='khan',
                xlabel='Number of trials',
                ylabel='p-value',
                title='A-B test random walk',
                )

def CheckCdf2():
    df = 3
    t = [SimulateChi2() for i in range(1000)]
    t2 = [scipy.stats.chi2.cdf(x, df) for x in t]
    cdf = Cdf.MakeCdfFromList(t2)

    myplot.Cdf(cdf, show=True)

def CheckCdf():
    xs, ys = Chi2Cdf(df=3, high=15)
    pyplot.plot(xs, ys)

    t = [SimulateChi2() for i in range(1000)]
    cdf = Cdf.MakeCdfFromList(t)

    myplot.Cdf(cdf, clf=False, show=True)


def Chi2Cdf(df=2, high=5, n=100):
    xs = numpy.linspace(0, high, n)
    ys = scipy.stats.chi2.cdf(xs, df)
    return xs, ys

def SimulateChi2(pa=0.5, q=0.5, n=100):
    expected = Expected(pa, q, n)

    simulated = Simulate(pa, q, n)

    chi2 = ChiSquared(expected, simulated)
    return chi2

if __name__ == "__main__":
    main()
