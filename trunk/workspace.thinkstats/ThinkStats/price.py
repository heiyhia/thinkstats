"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import myplot
import thinkbayes

FORMATS = ['png']


class Price(thinkbayes.Suite):

    def __init__(self, error_sigma):
        """Constructs the suite.

        error_sigma: standard deviation of the distribution of error
        """
        thinkbayes.Suite.__init__(self)
        pmf = thinkbayes.MakeGaussianPmf(35000, 7500, num_sigmas=4)

        # copy items from pmf to self
        for val, prob in pmf.Items():
            self.Set(val, prob)

        # store error_sigma for use in Likelihood
        self.error_sigma = error_sigma

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: actual price
        data: my guess
        """
        actual_price = hypo
        my_guess = data

        error = my_guess - actual_price
        like = thinkbayes.EvalGaussianPdf(
            mu=0, 
            sigma=self.error_sigma,
            x=error)

        return like


def main():
    
    # compute the stddev of the total error
    error_snowmobile = 500
    error_trip = 3000
    error_total = math.sqrt(error_snowmobile**2 + error_trip**2)
    print error_total

    myplot.Clf()

    suite = Price(error_total)
    suite.name = 'prior'
    myplot.Pmf(suite)

    my_guess = 15000
    suite.Update(my_guess)
    suite.name = 'posterior'
    myplot.Pmf(suite)

    print 'Posterior mean', suite.Mean()

    myplot.Save(root='price1',
                xlabel='price',
                formats=FORMATS)


if __name__ == '__main__':
    main()
