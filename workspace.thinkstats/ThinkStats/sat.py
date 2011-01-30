"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import datetime
import math
import sys

import Cdf
import myplot
import Pmf
import thinkstats


def ReadScale(filename='sat_scale.csv'):
    """Reads a CSV file of SAT scales (maps from raw score to standard score.

    Args:
      filename: string filename

    Returns:
      list of (raw score, standardize score) pairs
    """
    fp = open(filename)
    reader = csv.reader(fp)
    raws = []
    scores = []

    for t in reader:
        if t[0] == 'Raw':
            continue

        raws.append(int(t[0]))
        scores.append(int(t[-1]))

    raws.sort()
    scores.sort()
    return thinkstats.Interpolator(raws, scores)


def ReadScores(filename='SATPercentileRanks2009.csv'):
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
            score, number = [int(x) for x in t]
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
        
        prob_admit = Logistic(z)
        new.Incr(val, prob * prob_admit)

        print val, z, prob_admit

    new.Normalize()
    mu, sigma = Summarize(new)
    return new

def StandardScore(val, mu, sigma):
    return (val-mu) / sigma

def Logistic(z):
    return 1 / (1 + math.exp(-z))

def main(script):

    scale = ReadScale()
    print scale.xs
    print scale.ys
    print scale.Lookup(53)
    print scale.Reverse(760)
    return

    # read 'em and sort 'em
    scores = ReadScores()
    pmf = Pmf.MakePmfFromDict(dict(scores))
    pmf.Normalize()

    admitted = ApplyLogistic(pmf)

    cdf1 = Cdf.MakeCdfFromPmf(pmf)
    cdf2 = Cdf.MakeCdfFromPmf(admitted)

    quartiles = cdf2.Percentile(25), cdf2.Percentile(50), cdf2.Percentile(75)
    print 'quartiles', quartiles
    return

    myplot.Cdfs([cdf1, cdf2],
               xlabel='score', 
               ylabel='CDF', 
               show=True)
    
if __name__ == '__main__':
    main(*sys.argv)
