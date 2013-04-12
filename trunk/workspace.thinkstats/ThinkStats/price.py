"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import math
import myplot
import numpy
import scipy.stats
import thinkbayes

FORMATS = ['png']

def ReadData(filename='showcases.2011.csv'):
    """Reads a CSV file of data.

    Args:
      filename: string filename

    Returns: sequence of (val1 val2 bid1 bid2 diff1 diff2) tuples
    """
    fp = open(filename)
    reader = csv.reader(fp)
    res = []

    for t in reader:
        heading = t[0]
        data = t[1:]
        try:
            data = [int(x) for x in data]
            print heading, data[0], len(data)
            res.append(data)
        except ValueError:
            pass

    fp.close()
    return zip(*res)
    

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


class Pdf(object):
    def Density(x):
        """Returns the Pdf evaluated at x."""

    def MakePmf(self, xs):
        pmf = thinkbayes.Pmf()
        for x in xs:
            pmf.Set(x, self.Density(x))
        pmf.Normalize()
        return pmf


class GaussianPdf(Pdf):
    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma
        
    def Density(self):
        return thinkbayes.EvalGaussianPdf(x, self.mu, self.sigma)


class EstimatedPdf(Pdf):
    def __init__(self, seq):
        """Estimates the density function based on a sample.

        seq: sequence of data
        """
        xs = numpy.array(seq, dtype=numpy.double)
        self.kde = scipy.stats.gaussian_kde(xs)

    def Density(self, x):
        return self.kde.evaluate(x)


class ReturnCalculator(object):

    def __init__(self, error_sigma=3041):
        """Constructs the calculator.

        error_sigma: the standard deviation of your opponent's error dist
        """
        self.error_sigma = error_sigma
        pmf = thinkbayes.MakeGaussianPmf(0, error_sigma, 4)
        RemoveNegatives(pmf)
        self.opponent_cdf = thinkbayes.MakeCdfFromPmf(pmf) 

    def PlotExpectedReturns(self, suite):
        low, high = 0, 40000
        bids = numpy.linspace(low, high, 51)

        returns = [self.ExpectedReturn(bid, suite) for bid in bids]

        for bid, ret in zip(bids, returns):
            print bid, ret

        myplot.Clf()
        myplot.Plot(bids, returns)
        myplot.Show()

    def ExpectedReturn(self, bid, suite):
        """Computes the expected return of a given bid.

        bid: your bid
        suite: posterior distribution of prices
        """
        total = 0
        for price, prob in sorted(suite.Items()):
            roi = self.Roi(bid, price)
            total += prob * roi
        return total

    def Roi(self, bid, price):
        """Computes the return of a bid, given the actual price.
        """
        # if you overbid, you get nothing
        if bid > price:
            return 0

        # otherwise compute the probability of winning
        diff = price - bid
        prob = self.ProbWin(price - bid)

        # if you are within 250 dollars, you win both showcases
        if diff > 250:
            return price * prob
        else:
            return 2 * price * prob

    def ProbWin(self, diff):
        """Computes the probability of winning for a given error.

        diff: how much your bid was off by
        """
        prob = 1 - self.opponent_cdf.Prob(diff)
        return prob


def RemoveNegatives(pmf):
    """Remove negative values from the Pmf."""
    for val in pmf.Values():
        if val < 0:
            pmf.Remove(val)
    pmf.Normalize()


def main():
    data = ReadData(filename='showcases.2011.csv')
    data += ReadData(filename='showcases.2012.csv')

    cols = zip(*data)
    val1, val2, bid1, bid2, diff1, diff2 = cols

    kde_val1 = EstimatedPdf(val1)

    low, high = 0, 75000
    xs = numpy.linspace(low, high, 51)
    ys = kde_val1.kde.evaluate(xs)

    pmf = kde_val1.MakePmf(xs)
    myplot.Pmf(pmf)
    myplot.Show()

    return


    cdf_val1 = thinkbayes.MakeCdfFromList(val1)
    cdf_val2 = thinkbayes.MakeCdfFromList(val2)

    print cdf_val1.Mean()
    print cdf_val2.Mean()

    cdf_diff1 = thinkbayes.MakeCdfFromList(diff1)
    cdf_diff2 = thinkbayes.MakeCdfFromList(diff2)

    print 'Prob diff1 <= -1', cdf_diff1.Prob(-1)
    print 'Prob diff2 <= -1', cdf_diff2.Prob(-1)


    myplot.Cdf(cdf_diff1)
    myplot.Cdf(cdf_diff2)
    myplot.Show()

    return

    for line in data:
        print line
    return


    calc = ReturnCalculator()

    # test ProbWin
    for diff in [0, 100, 1000, 10000]:
        print diff, calc.ProbWin(diff)
    print

    # test Roi
    price = 18000
    for bid in [15000, 16000, 17000, 17500, 17800, 18001]:
        print bid, calc.Roi(bid, price)
    print

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

    calc.PlotExpectedReturns(suite)


if __name__ == '__main__':
    main()
