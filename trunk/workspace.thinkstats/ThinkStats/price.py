"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import myplot
import thinkbayes

class Price(thinkbayes.Suite):

    def __init__(self, error_sigma):
        """Constructs the suite.

        hypos: sequence of hypotheses
        error_sigma: standard deviation of the distribution of error
        """
        thinkbayes.Suite.__init__(self)
        pmf = thinkbayes.MakeGaussianPmf(35000, 7500, 4)
        for val, prob in pmf.Items():
            self.Set(val, prob)

        self.error_sigma = error_sigma

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: actual price
        data: my guess
        """
        actual_price = hypo
        my_guess = data

        x = my_guess - actual_price
        like = thinkbayes.EvalGaussianPdf(mu=0, sigma=self.error_sigma, x=x)

        return like


def main():
    sig1 = 500
    sig2 = 3000
    error_sigma = math.sqrt(sig1**2 + sig2**2)
    print error_sigma

    myplot.Clf()

    suite = Price(error_sigma)
    suite.name = 'prior'
    myplot.Pmf(suite)

    suite.Update(15000)
    suite.name = 'posterior'
    myplot.Pmf(suite)

    myplot.Show()


if __name__ == '__main__':
    main()
