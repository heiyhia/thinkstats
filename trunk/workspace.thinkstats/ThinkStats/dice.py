"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from thinkbayes import Suite


class Dice(Suite):
    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer number of sides on the die
        data: integer die roll
        """
        if hypo < data:
            return 0
        else:
            return 1.0/hypo


def main():
    suite = Dice([4, 6, 8, 20])

    for roll in [7]:
        suite.Update(roll)

    print 'Posterior'
    suite.Print()

if __name__ == '__main__':
    main()
