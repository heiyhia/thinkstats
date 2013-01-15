"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from dice import Dice
import myplot


class Train(Dice):
    """The likelihood function for the train problem is the same as
    for the Dice problem."""


def main():
    hypos = xrange(1, 1001)
    suite = Train(hypos)

    suite.Update(60)
    print suite.Mean()

    myplot.PrePlot(1)
    myplot.Pmf(suite)
    myplot.Save(root='train1',
                xlabel='Number of trains',
                ylabel='Probability',
                formats=['pdf', 'eps'])


if __name__ == '__main__':
    main()
