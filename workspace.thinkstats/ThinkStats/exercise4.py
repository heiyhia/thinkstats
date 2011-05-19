"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import myplot
import Pmf
import Cdf

def BiasPmf(pmf, name):
    """Returns the Pmf with oversampling proportional to value.

    If pmf is the distribution of true values, the result is the
    distribution that would be seen if values are oversampled in
    proportion to their values; for example, if you ask students
    how big their classes are, large classes are oversampled in
    proportion to their size.

    Args:
      pmf: Pmf object.

     Returns:
       Pmf object
    """
    new_pmf = pmf.Copy()
    new_pmf.name = name

    for x, p in pmf.Items():
        new_pmf.Mult(x, x)

        # NOTE: this function is incomplete; you need to finish it!
        # See http://greenteapress.com/thinkstats/Pmf.html for information
        # about Pmf operations.
    new_pmf.Normalize()
    return new_pmf


def ClassSizes():

    # start with the actual distribution of class sizes from the book
    d = {
         7: 8, 
         12: 8, 
         17: 14, 
         22: 4, 
         27: 6, 
         32: 12, 
         37: 8, 
         42: 3, 
         47: 2, 
    }

    # form the pmf
    pmf = Pmf.MakePmfFromDict(d, 'actual')
    pmf.Normalize()
    print 'mean', pmf.Mean()
    print 'var', pmf.Var()
    
    biased = BiasPmf(pmf, name='biased')
    print 'mean', biased.Mean()
    print 'var', biased.Var()

    # plot the actual pmf
    myplot.Pmfs([pmf, biased], 
               show=True, 
               xlabel='Class size',
               ylabel='PMF')

 
def main():
    ClassSizes()


if __name__ == '__main__':
    main()
