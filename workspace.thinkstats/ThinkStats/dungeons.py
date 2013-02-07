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


def PmfMax(pmf1, pmf2):
    """Computes the distribution of the max of values drawn from two Pmfs.

    pmf1, pmf2: Pmf objects

    returns: new Pmf
    """
    res = thinkbayes.Pmf()
    for v1, p1 in pmf1.Items():
        for v2, p2 in pmf2.Items():
            res.Incr(max(v1, v2), p1*p2)
    return res
    

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

    myplot.PrePlot(num=2)
    myplot.Pmf(three)
    myplot.Pmf(three_exact, linestyle='dashed')
    myplot.Save(root='dungeons1',
                xlabel='Sum of three d6',
                ylabel='Probability',
                axis=[2, 19, 0, 0.15],
                formats=['pdf', 'eps'])

    myplot.Clf()
    myplot.PrePlot(num=1)
    
    # compute the distribution of the best attribute the hard way
    best_attr2 = PmfMax(three_exact, three_exact)
    best_attr4 = PmfMax(best_attr2, best_attr2)
    best_attr6 = PmfMax(best_attr4, best_attr2)
    # myplot.Pmf(best_attr6)

    # and the easy way
    best_attr_cdf = three_exact.Max(6)
    best_attr_cdf.name = ''
    best_attr_pmf = thinkbayes.MakePmfFromCdf(best_attr_cdf)
    best_attr_pmf.Print()

    myplot.Pmf(best_attr_pmf)
    myplot.Save(root='dungeons2',
                xlabel='Sum of three d6',
                ylabel='Probability',
                axis=[2, 19, 0, 0.23],
                formats=['pdf', 'eps'])
    

if __name__ == '__main__':
    main()
