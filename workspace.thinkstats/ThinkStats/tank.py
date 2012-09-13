"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from dice import Dice
import myplot

class Tank(Dice):
    """The likelihood function for the tank problem is the same as
    for the Dice problem."""


def main():
    hypos = xrange(1, 200)
    suite = Tank(hypos)

    suite.Update(30)
    print 'After one tank'

    myplot.Pmf(suite)
    myplot.Show(xlabel='Number of tanks',
                ylabel='Probability')


if __name__ == '__main__':
    main()
