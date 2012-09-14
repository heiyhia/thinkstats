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


def Mean(suite):
    total = 0
    for hypo, prob in suite.Items():
        total += hypo * prob
    return total


def main():
    hypos = (range(1, 10) + 
             range(10, 100, 10) + 
             range(100, 1000, 100) +
             range(1000, 10000, 1000))

    suite = Train(hypos)

    suite.Update(30)

    myplot.Pmf(suite)
    myplot.Save(root='train2',
                xlabel='Number of trains',
                ylabel='Probability',
                xscale='log')

    print Mean(suite)
    print suite.Mean()

if __name__ == '__main__':
    main()
