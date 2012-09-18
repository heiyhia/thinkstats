"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from thinkbayes import Pmf
from dice import Dice
import myplot


class Train(Dice):
    """The likelihood function for the train problem is the same as
    for the Dice problem."""

    def __init__(self, hypos, alpha=1.0):
        """Initializes the hypotheses with a power law distribution.

        hypos: sequence of hypotheses
        """
        Pmf.__init__(self)
        for hypo in hypos:
            self.Set(hypo, hypo**(-alpha))
        self.Normalize()



def Mean(suite):
    total = 0
    for hypo, prob in suite.Items():
        total += hypo * prob
    return total


def MakePosterior(high, dataset):
    hypos = xrange(1, high+1)
    suite = Train(hypos)
    suite.name = str(high)

    for data in dataset:
        suite.Update(data)

    myplot.Pmf(suite)
    return suite


def Percentile(pmf, p):
    """Computes a percentile of a given Pmf.

    p: float probability
    """
    total = 0
    for val, prob in pmf.Items():
        total += prob
        if total >= p:
            return val    


def main():
    dataset = [30, 60, 90]

    for high in [500, 1000, 2000]:
        suite = MakePosterior(high, dataset)
        print high, suite.Mean()

    myplot.Save(root='train3',
                xlabel='Number of trains',
                ylabel='Probability')

    interval = Percentile(suite, 0.05), Percentile(suite, 0.95)
    print interval


if __name__ == '__main__':
    main()
