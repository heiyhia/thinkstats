"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from thinkbayes import Suite


class Dice(Suite):
    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer number of sides on the die
        data: integer die roll
        """
        num_sides = hypo
        outcome = data

        # TODO: write this method!
        return 1


def main():
    suite = Dice([4, 6, 8, 12, 20])

    suite.Update(6)
    print 'After one 6'
    suite.Print()


if __name__ == '__main__':
    main()
