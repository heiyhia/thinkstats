"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import first_baby
import Pmf
import Cdf

class Pregnancies(first_baby.Pregnancies):

    def Process(self):
        self.lengths = [p.prglength for p in self.records]
        self.mu = Mean(self.lengths)
        self.var = Var(self.lengths, self.mu)

        self.lengths.sort()
        self.trim = TrimmedMean(self.lengths)

    def MakeDist(self, name=''):
        """Makes the histogram, PMF and CDF of self.lengths."""
        self.hist = Pmf.MakeHist(self.lengths, name=name)
        self.pmf = self.hist.MakePmf()
        self.cdf = Cdf.MakeCdfFromDict(self.hist.GetDict(), name=name)


def Mean(t):
    """Computes the mean of a sequence of numbers.

    Args:
        t: sequence of numbers

    Returns:
        float
    """
    return float(sum(t)) / len(t)


def TrimmedMean(t, p=0.01):
    """Computes the trimmed mean of a sequence of numbers.

    Args:
        t: sequence of numbers

        p: fraction of values to trim off each end

    Returns:
        float
    """
    n = int(p * len(t))
    t = t[n:-n]
    return Mean(t)


def Var(t, mu=None):
    """Computes the mean of a sequence of numbers.

    Args:
        t: sequence of numbers

    Returns:
        float
    """
    if mu is None:
        mu = Mean(t)

    # compute the squared deviations and return their mean.
    dev2 = [(x - mu)**2 for x in t]
    var = Mean(dev2)
    return var

