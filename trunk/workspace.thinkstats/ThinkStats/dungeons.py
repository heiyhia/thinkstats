"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import random

import thinkbayes
import myplot

class Die(thinkbayes.Pmf):
    def __init__(self, sides, name=''):
        d = dict((i, 1) for i in xrange(1, sides+1))
        thinkbayes.Pmf.__init__(self, d, name)
        self.Normalize()


def main():
    random.seed(17)

    d6 = Die(6, 'd6')

    dice = [d6] * 3
    three = thinkbayes.SampleSum(dice, 1000)
    three.name = 'sample'
    three.Print()

    three_exact = d6 + d6 + d6
    three_exact.name = 'exact'
    three_exact.Print()

    myplot.Pmf(three)
    myplot.Pmf(three_exact)
    myplot.Save(root='dungeons1',
                xlabel='Sum of three d6',
                ylabel='Probability',
                axis=[2, 19, 0, 0.15],
                formats=['pdf', 'eps'])


if __name__ == '__main__':
    main()
