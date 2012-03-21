"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import bayes
import copy
import string
import matplotlib.pyplot as pyplot
import myplot
import Pmf
import Cdf
import math
import random
import sys

lowercase = string.ascii_lowercase


class Beta(object):
    """Represents a Beta distribution.

    See http://en.wikipedia.org/wiki/Beta_distribution

    alpha and beta are the parameters 

    """
    def __init__(self, yes, no, name=''):
        """Initializes a Beta distribution.

        x, -1 yields a distribution always near 1

        yes and no are the number of successful/unsuccessful trials.
        """
        if yes == -1:
            yes, no = 0, 999999

        self.alpha = yes+1
        self.beta = no+1

    def Update(self, yes, no):
        """Updates a Beta distribution."""
        self.alpha += yes
        self.beta += no
        
    def Mean(self):
        """Computes the mean of this distribution."""
        return float(self.alpha) / (self.alpha + self.beta)

    def Random(self):
        """Generates a random variate from this distribution."""
        return random.betavariate(self.alpha, self.beta)

    def Pdf(self, p):
        """Computes the PDF at p."""
        return math.pow(p, self.alpha-1) * math.pow(1-p, self.beta-1)
        
    def Pmf(self, steps=1001):
        """Returns the PDF of this distribution."""
        ps = [i / (steps-1.0) for i in xrange(steps)]
        probs = [self.Pdf(p) for p in ps]
        pmf = Pmf.MakePmfFromDict(dict(zip(ps, probs)))
        return pmf

    def Cdf(self, steps=1001):
        """Returns the CDF of this distribution."""
        pmf = self.Pmf(steps=steps)
        cdf = Cdf.MakeCdfFromPmf(pmf)
        return cdf

    def ConditionalCdf(self, fraction, steps=1001):
        """Generates a CDF conditioned on p <= fraction."""
        ps = [fraction * i / (steps-1.0) for i in xrange(steps)]
        probs = [self.Pdf(p) for p in ps]
        cdf = Cdf.MakeCdfFromItems(zip(ps, probs))
        return cdf                                                    


def main(script, *args):
    beta = Beta(1, 0)
    pmf = beta.Pmf()
    myplot.Pmf(pmf, show=True,
               xlabel='Probability of sunrise: p',
               ylabel='Probability density',
               title='Beta distribution')

    cdf = beta.Cdf()
    print cdf.Percentile(5)
    print cdf.Percentile(95)

    print cdf.Prob(0.5)

if __name__ == '__main__':
    main(*sys.argv)
