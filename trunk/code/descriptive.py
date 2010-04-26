"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import first_baby
import Pmf

class Pregnancies(first_baby.Pregnancies):

    def Process(self):
        self.lengths = [p.prglength for p in self.records]
        self.mu = Mean(self.lengths)
        self.var = Var(self.lengths, self.mu)

    def MakePmf(self, name=''):
        self.hist = Pmf.MakeHist(self.lengths, name=name)
        self.pmf = self.hist.MakePmf()


def Mean(t):
    """Computes the mean of a sequence of numbers.

    Args:
        t: sequence of numbers

    Returns:
        float
    """
    return float(sum(t)) / len(t)


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

