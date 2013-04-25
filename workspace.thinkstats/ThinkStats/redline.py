"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import thinkbayes

import matplotlib.pyplot as pyplot
import myplot
import numpy

import csv
import math
import random
import sys
import time

from math import log

formats = ['pdf']

"""
Notation guide:

z: time between trains
x: time since the last train
y: time until the next train

z': distribution of z as seen by a random arrival

"""


# interarrival_times in seconds, collected by Kai Austin and Brendan Ritter
# using http://developer.mbta.com/Data/Red.txt

interarrival_times = [
    365, 47, 146, 545, 132, 270, 443, 190, 372, 167, 
    375, 128, 455, 262, 233, 386, 561, 386, 241, 562, 489, 
    455, 410, 489, 342, 347, 170, 260, 148, 375, 312, 265, 
    201, 245, 567, 633, 627, 643, 159, 827, 773, 159, 1100, 160, 
    148, 187, 290, 353, 133, 180, 355, 151, 558, 220, 232, 
    353, 199, 160, 172
]

def BiasPmf(pmf, name='', invert=False):
    """Returns the Pmf with oversampling proportional to value.

    If pmf is the distribution of true values, the result is the
    distribution that would be seen if values are oversampled in
    proportion to their values; for example, if you ask students
    how big their classes are, large classes are oversampled in
    proportion to their size.

    If invert=True, computes in inverse operation; for example,
    unbiasing a sample collected from students.

    Args:
      pmf: Pmf object.
      name: string name for the new Pmf.
      invert: boolean

     Returns:
       Pmf object
    """
    new_pmf = pmf.Copy(name=name)

    for x, p in pmf.Items():
        if invert:
            new_pmf.Mult(x, 1.0/x)
        else:
            new_pmf.Mult(x, x)
        
    new_pmf.Normalize()
    return new_pmf


def UnbiasPmf(pmf, name):
    """Returns the Pmf with oversampling proportional to 1/value.

    Args:
      pmf: Pmf object.
      name: string name for the new Pmf.

     Returns:
       Pmf object
    """
    return BiasPmf(pmf, name, invert=True)


def WeibullSample(n, lam, k):
    return [random.weibullvariate(lam, k) for i in range(n)]


def ActualInterarrivals(n, lam=12, k=1.5):
    sample = WeibullSample(n, lam, k)
    pmf = thinkbayes.MakePmfFromList(sample)
    pmf.name = 'actual'
    
def BiasedInterarrivals(n, lam=12, k=1.5):
    actual_pmf = ActualInterarrivals(1000)
    observed_pmf = BiasPmf(actual_pmf, 'observed')


def MakeUniformPmf(low, high):
    """Make a uniform Pmf.

    low: lowest value (inclusive)
    high: highest value (inclusize)
    """
    pmf = thinkbayes.Pmf()
    for x in range(low, high+1, 10):
        pmf.Set(x, 1)
    pmf.Normalize()
    return pmf    
    

def PmfOfWaitTime(pmf_z):
    pmf_zp = BiasPmf(pmf_z, name="z'")

    meta_pmf = thinkbayes.Pmf()
    for interarrival, prob in pmf_zp.Items():
        uniform = MakeUniformPmf(0, interarrival)
        meta_pmf.Set(uniform, prob)

    pmf_y = thinkbayes.MakeMixture(meta_pmf, name='y')
    return pmf_zp, pmf_y


def main(script):
    pdf_z = thinkbayes.EstimatedPdf(interarrival_times)
    low, high = 0, 1300
    xs = range(low, high+1, 10)
    pmf_z = pdf_z.MakePmf(xs, name='z')

    pmf_zp, pmf_y = PmfOfWaitTime(pmf_z)

    myplot.Pmf(pmf_z)
    myplot.Pmf(pmf_zp)
    myplot.Pmf(pmf_y)
    myplot.Show()

    print 'Mean z', pmf_z.Mean()
    print 'Mean zp', pmf_zp.Mean()
    print 'Mean y', pmf_y.Mean()

    # arrival rate of passengers in passengers / second
    lam_pass = 0.05
    
    
    
    return



    n = 300
    lam = 12
    k = 0.9
    sample = WeibullSample(n, lam, k)
    cdf = thinkbayes.MakeCdfFromList(sample)
    myplot.Cdf(cdf)
    myplot.Show()

    print 'mean', cdf.Mean()
    print 'median', cdf.Percentile(50)
    print '25th', cdf.Percentile(25)
    print '75th', cdf.Percentile(75)

if __name__ == '__main__':
    main(*sys.argv)
