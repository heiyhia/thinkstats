"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import numpy
import math
import scipy

import thinkbayes
import thinkplot

FORMATS = ['png']


def FakeData(n, p):
    """Generate a binomial variate with parameters n, p.

    n: int number of values
    p: float prob
    
    returns: tuple of successes, failures
    """
    k = numpy.random.binomial(n, p)
    return k, n-k


def MakePosterior(data, preload=20, steps=201, name=""):
    """Makes the posterior distribution of CTR given the data.

    Where CTR is "click through rate"

    data: sequence of values in [0, 1]
    preload: number of preloaded "no"
    steps: number of values to put in the posterior dist
    name: name to give the pmf

    returns: Pmf of CTR
    """
    yes, no = data
    pmf = thinkbayes.Beta(yes+1, no+1).MakePmf(steps, name)
    return pmf


def SampleCtr(n=1, preload=20):
    """Draw a sample from a beta distribution of CTRs.

    n: sample size
    preload: number of preloaded "no"

    returns: sequence of CTR
    """
    beta = thinkbayes.Beta(1, preload+1)
    ps = beta.Sample(n)
    return ps


def PosteriorProb(data1, data2):
    """Computes the posterior probability of A>B, given data.

    data1: sequence of [0, 1]
    data2: sequence of [0, 1]

    returns: float prob
    """
    pmf1 = MakePosterior(data1)
    pmf2 = MakePosterior(data2)

    gt = (pmf1 > pmf2)
    eq = (pmf1 == pmf2)

    return gt + eq/2


def PlotPosteriors():
    thinkbayes.RandomSeed(18)

    data1 = FakeData(100, 0.03)
    data2 = FakeData(100, 0.05)

    pmf1 = MakePosterior(data1, name="headline a")
    pmf2 = MakePosterior(data2, name="headline b")

    lt = pmf1 < pmf2
    eq = pmf1 == pmf2
    gt = pmf1 > pmf2

    print lt + eq/2
    print gt + eq/2

    cdf1 = pmf1.MakeCdf()
    cdf2 = pmf2.MakeCdf()
    thinkplot.PrePlot(num=2)
    thinkplot.Cdfs([cdf1, cdf2])
    thinkplot.Show(axis=[0, 0.2, 0, 1])


def RunSimulation(ctr1, ctr2):
    """Generate data and compute posterior prob of A>B.

    ctr1: float CTR for A
    ctr2: float CTR for B

    returns: float posterior prob A>B
    """
    data1 = FakeData(100, ctr1)
    data2 = FakeData(100, ctr2)

    p = PosteriorProb(data1, data2)

    #print data1[0], data2[0], p
    return p


def RunPvalue1(ctr1, ctr2):
    """Generate data and return the p-value of the observed difference.

    ctr1: float CTR for A
    ctr2: float CTR for B

    returns: float p-value
    """
    def proportion(data):
        yes, no = data
        n = yes + no
        p = float(yes) / n
        return p, n

    data1 = FakeData(100, ctr1)
    data2 = FakeData(100, ctr2)

    diff = abs(data1[0] - data2[0])

    p1, n1 = proportion(data1)
    p2, n2 = proportion(data2)
    p = float(p1 * n1 + p2 * n2) / (n1 + n2)
    term = 1.0/n1 + 1.0/n2
    se = math.sqrt(p * (1-p) * term)
    z = -abs(p1 - p2) / se
    pval = 2 * scipy.stats.norm.cdf(z)

    print p, se, z, pval

    return pval


def RunPvalue2(ctr1, ctr2):
    """Generate data and return the p-value of the observed difference.

    ctr1: float CTR for A
    ctr2: float CTR for B

    returns: float p-value
    """
    data1 = FakeData(100, ctr1)
    data2 = FakeData(100, ctr2)

    diff = abs(data1[0] - data2[0])

    pval = 1 - sample_diff.Prob(diff)
    print data1[0], data2[0], pval
    print p, se, z, pval_z

    return pval


def SampleDistOfDiff(ctr1, ctr2, n=2000):
    diffs = []
    for i in range(n):
        data1 = FakeData(100, ctr1)
        data2 = FakeData(100, ctr2)
        diff = abs(data1[0] - data2[0])
        diffs.append(diff)

    cdf = thinkbayes.MakeCdfFromList(diffs, name='diff')
    return cdf
    

def PredDist(ctr1, ctr2, n=100):
    """Predictive posterior distribution of prob A>B.

    ctr1: float CTR for A
    ctr2: float CTR for B
    n: number of simulations to run

    returns: Cdf of posterior probs
    """
    data1 = FakeData(100, ctr1)
    data2 = FakeData(100, ctr2)

    cdf1 = MakePosterior(data1).MakeCdf()
    cdf2 = MakePosterior(data2).MakeCdf()

    ctr1s = cdf1.Sample(n)
    ctr2s = cdf2.Sample(n)

    ps = [RunSimulation(q1, q2) for q1, q2 in zip(ctr1s, ctr2s)]

    pred = thinkbayes.MakeCdfFromList(ps)
    return pred


def SamplePredDist(ctr1, ctr2, n=30):    
    """Computes the sample distribution of p.

    Where p is the predictive posterior probability of A>B.

    ctr1: CTR of A
    ctr2: CTR of B
    n: number of iterations

    returns: Cdf of p
    """
    ps = []
    for i in range(n):
        pred = PredDist(ctr1, ctr2)
        p = pred.Mean()
        ps.append(p)

    sample_pred = thinkbayes.MakeCdfFromList(ps, name='pred means')
    return sample_pred
    

def SampleDist(ctr1, ctr2, n=365):
    """Computes the sample distribution of p.

    Where p is the posterior probability of A>B.

    ctr1: CTR of A
    ctr2: CTR of B
    n: number of iterations

    returns: Cdf of p
    """
    ps = [RunSimulation(ctr1, ctr2) for i in range(n)]
    cdf = thinkbayes.MakeCdfFromList(ps, name='posterior probs')
    return cdf


def SampleDistPval(ctr1, ctr2, n=365):
    """Computes the sample distribution of p-values.

    ctr1: CTR of A
    ctr2: CTR of B
    n: number of iterations

    returns: Cdf of p-values
    """
    ps = [RunPvalue1(ctr1, ctr2) for i in range(n)]
    cdf = thinkbayes.MakeCdfFromList(ps, name='p-values')
    return cdf


def RunScenario(ctr1, ctr2):
    """
    """
    pred = PredDist(ctr1, ctr2)
    p = pred.Mean()

    odds = thinkbayes.Odds(p)
    
    actual_p = RunSimulation(ctr1, ctr2)
    if random.random() < p:
        payout = odds
    else:
        payout = 0
            

def main():
    ctr1 = 0.05
    ctr2 = 0.05

    global sample_diff
    sample_diff = SampleDistOfDiff(ctr1, ctr2, n=2000)
    #thinkplot.Cdf(sample_diff)
    #thinkplot.Show()

    sample_pval = SampleDistPval(ctr1, ctr2)

    thinkplot.Cdf(sample_pval)
    thinkplot.Save(root='abtest3',
                   xlabel='p-value',
                   ylabel='CDF',
                   formats=FORMATS)

    return

    sample_dist = SampleDist(ctr1, ctr2)

    thinkplot.Cdf(sample_dist)
    thinkplot.Save(root='abtest1',
                   xlabel='prob A > B',
                   ylabel='CDF',
                   formats=FORMATS)

    return

    sample_pred = SamplePredDist(ctr1, ctr2)

    thinkplot.Cdf(sample_dist)
    thinkplot.Cdf(sample_pred)
    thinkplot.Save(root='abtest2',
                   xlabel='prob A > B',
                   ylabel='CDF',
                   formats=FORMATS)

    return

    # plot the prior distribution of CTR
    ps = SampleCtr(100)
    cdf = thinkbayes.MakeCdfFromList(ps)
    thinkplot.Cdf(cdf)
    thinkplot.Show()







if __name__ == "__main__":
    main()
