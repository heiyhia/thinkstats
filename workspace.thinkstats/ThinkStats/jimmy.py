"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import descriptive
import itertools
import Pmf

import Cdf
import random
import thinkstats


def SimulateSample(cdf, num_nuts):
    """Generates a Hist of simulated nuts.
    
    Args:
      cdf: 
      num_nuts: number of times to nuts

    Returns:
      Hist object
    """
    t = cdf.Sample(num_nuts)
    hist = Pmf.MakeHistFromList(t)
    return hist


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


def Test(expected, observed, num_trials=1000):
    """Run a simulation to estimate the p-value of the observed values.

    Args:
      expected: Hist of expected values
      observed: Hist of observed values
      num_trials: how many simulations to run

    Returns:
      float p-value
    """

    # compute the chi-squared stat
    threshold = ChiSquared(expected, observed)
    print 'chi-squared', threshold

    print 'simulated %d trials' % num_trials
    chi2s = []
    count = 0
    num_nuts = observed.Total()
    cdf = Cdf.MakeCdfFromHist(expected)

    for _ in range(num_trials):
        simulated = SimulateSample(cdf, num_nuts)
        chi2 = ChiSquared(expected, simulated)
        chi2s.append(chi2)
        if chi2 >= threshold:
            count += 1
            
    print 'max chi2', max(chi2s)
    
    pvalue = 1.0 * count / num_trials
    print 'p-value', pvalue

    return pvalue


def ConvertToCount(sample, count_per):
    """Convert from weight to count.

    sample: Hist or Pmf that maps from category to weight in pounds
    count_per: dictionary that maps from category to count per ounce
    """
    for value, count in sample.Items():
        sample.Mult(value, 16 * count_per[value])

def MakeVat(expected, num_nuts, factor=10, stir=0.0):
    t = []
    for value, freq in expected.Items():
        t.extend([value] * freq * factor)

    if stir == -1:
        random.shuffle(t)
    else:
        [RandomSwap(t) for i in xrange(int(num_nuts*stir))]

    
    return t

def RandomSwap(t):
    i, j = [random.randint(len(t)) for i in range(2)]
    t[i], t[j] = t[j], t[i]

def PercentAdjacent(t):
    for i in range(len(t) - 1):
        if t[i] == t[i+1]:
            

def main():
    # make a Hist of observed values
    count_per = dict(cashew=17, brazil=7, almond=22, peanut=28)
    count_per = dict(cashew=12, brazil=5, almond=18, peanut=24)
    sample = dict(cashew=6, brazil=3, almond=5, peanut=6)


    observed = Pmf.MakeHistFromDict(sample)
    ConvertToCount(observed, count_per)

    advertised = dict(cashew=40, brazil=15, almond=20, peanut=25)
    expected = Pmf.MakePmfFromDict(advertised)
    ConvertToCount(expected, count_per)
    expected.Normalize(observed.Total())

    for value, e in expected.Items():
        o = observed.Freq(value)
        print value, e, o, 100 * (o - e)/ e, '%'

    # compute the chi-squared statistic
    print ChiSquared(expected, observed)

    Test(expected, observed, num_trials=1000)

if __name__ == "__main__":
    main()
