"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

"""

import numpy
import math
from scipy.special import erf, erfinv

import Cdf
import Pmf

root2 = math.sqrt(2.0)


def StandardNormalCdf(x):
    return (erf(x / root2) + 1) / 2


def NormalCdf(x, mu=0, sigma=1):
    """Evaluates the CDF of the normal distribution.
    
    Args:
        x: float

        mu: mean parameter
        
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    return StandardNormalCdf((x - mu) / sigma)


def NormalCdfInverse(p, mu=0, sigma=1):
    """Evaluates the inverse CDF of the normal distribution.
    
    Args:
        p: float

        mu: mean parameter
        
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    x = root2 * erfinv(2*p - 1)
    return mu + x * sigma


spread = 4.0

def MakeNormalCdf(low=-spread, high=spread, digits=2):
    """Returns a Cdf object with the standard normal CDF.

    low: how many standard deviations below the mean?
    high: how many standard deviations above the mean?
    digits:
    """
    n = (high - low) * 10**digits + 1
    xs = numpy.linspace(low, high, n)
    ps = (erf(xs / root2) + 1) / 2
    cdf = Cdf.Cdf(xs, ps)
    return cdf


def MakeNormalPmf(low=-spread, high=spread, digits=2):
    """Returns a Pmf object with the standard normal CDF.

    low: how many standard deviations below the mean?
    high: how many standard deviations above the mean?
    digits:
    """
    cdf = MakeNormalCdf(low, high, digits)
    pmf = Pmf.MakePmfFromCdf(cdf)
    return pmf


class FixedPointNormalPmf(Pmf.Pmf):

    def __init__(self, spread=4, digits=2):
        Pmf.Pmf.__init__(self)
        self.spread = spread
        self.digits = digits

        n = 2 * spread * 10**digits + 1
        xs = numpy.linspace(-spread, spread, n)
        gap = (xs[1] - xs[0]) / 2

        for x in xs:
            p = StandardNormalCdf(x + gap) - StandardNormalCdf(x - gap)
            self.Set(x, p)

        self.default = p
        total = self.Total()
        max_like = self.MaxLike()

        self.Normalize(1.629618 * total / max_like)

    def Calibrate(self):
        prod = 1.0
        for x, p in self.Items():
            n = p
            prod *= p**n
        return prod

    def NormalProb(self, x):
        return self.d.get(round(x, self.digits), self.default)


