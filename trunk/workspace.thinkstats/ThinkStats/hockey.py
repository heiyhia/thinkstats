"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import numpy

import myplot
from thinkbayes import Suite


def norm_pdf(x, mu, sigma):
    """Computes the unnormalized PDF of the normal distribution.

    x: value
    mu: mean
    sigma: standard deviation
    
    returns: float probability density (unnormalized)
    """
    z = (x - mu) / sigma
    p = math.exp(-z**2/2)
    return p


class Hockey(Suite):
    def __init__(self):
        Suite.__init__(self)
        for x in numpy.arange(1.5, 3.9, 0.05):
            p = norm_pdf(x, 2.7, 0.3)
            self.Set(x, p)
            
    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        Evaluates the Poisson PMF for lambda, k and t.

        I dropped the k! term in the denominator because it does not
        depend on lambda.

        hypo: goal scoring rate in goals per game
        data: goals scored in one period
        """
        t = 1.0 / 3.0
        lam = hypo
        k = data
        like = (lam*t)**k * math.exp(lam*t)
        return like


def main():
    suite1 = Hockey()
    suite1.name = 'bruins'
    suite1.Update(1)
    myplot.Pmf(suite1)
    
    suite2 = Hockey()
    suite2.name = 'sabres'
    suite2.Update(0)
    myplot.Pmf(suite2)
    
    myplot.Show()


if __name__ == '__main__':
    main()
