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

    def __init__(self, pmf, player):
        """Constructs the suite.

        pmf: prior distribution of value
        player: Player object containing the Pdf of underness
        """
        thinkbayes.Suite.__init__(self)

        # copy items from pmf to self
        for val, prob in pmf.Items():
            self.Set(val, prob)

        self.player = player

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: actual price
        data: my guess
        """
        actual_price = hypo
        my_guess = data

        underness = my_guess - actual_price
        like = self.player.UndernessDensity(underness)

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

    def MakePmf(self, xs):
        ps = self.kde.evaluate(xs)
        pmf = thinkbayes.MakePmfFromItems(zip(xs, ps))
        return pmf


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


class Player(object):

    low, high = 0, 75000
    xs = numpy.linspace(low, high, 51)
    diffs = numpy.linspace(-30000, 30000, 51)

    def __init__(self, val, bid, diff):
        self.val = val
        self.big = bid
        self.diff = diff

        self.cdf_val = thinkbayes.MakeCdfFromList(val)
        self.cdf_diff = thinkbayes.MakeCdfFromList(diff)

        self.kde_val = EstimatedPdf(val)
        self.kde_diff = EstimatedPdf(diff)
        # self.pmf_diff = kde_diff.MakePmf(self.xs)

    def UndernessDensity(self, underness):
        return self.kde_diff.Density(underness)

    def PmfVal(self):
        return self.kde_val.MakePmf(self.xs)

    def PmfDiff(self):
        return self.kde_diff.MakePmf(self.diffs)

    def CdfDiff(self):
        return self.cdf_diff

    def MakePosterior(estimated_value):
        pmf = thinkbayes.MakePmfFromList(self.val)
        self.prior = thinkbayes.Price(pmf)
        self.posterior = self.prior.Copy()
        self.posterior.Update(estimated_value)

    def OptimalBid(self, opponent):
        """
        
        Precondition: self.posterior has been computed
        """

def MakePlots(player1, player2):
    myplot.Clf()
    myplot.PrePlot(num=2)
    pmf1 = player1.PmfVal()
    pmf1.name = 'Player 1 prior'
    pmf2 = player2.PmfVal()
    pmf2.name = 'Player 2 prior'
    myplot.Pmfs([pmf1, pmf2])
    myplot.Save(root='price1',
                xlabel='price ($)',
                formats=FORMATS)

    myplot.Clf()
    myplot.PrePlot(num=2)
    cdf1 = player1.CdfDiff()
    cdf1.name = 'Player 1 undernesss'
    cdf2 = player2.CdfDiff()
    cdf2.name = 'Player 2 undernesss'
    myplot.Cdfs([cdf1, cdf2])
    myplot.Save(root='price2',
                xlabel='underness ($)',
                formats=FORMATS)


def main():
    data = ReadData(filename='showcases.2011.csv')
    data += ReadData(filename='showcases.2012.csv')

    cols = zip(*data)
    val1, val2, bid1, bid2, diff1, diff2 = cols

    player1 = Player(val1, bid1, diff1)
    player2 = Player(val2, bid2, diff2)

    MakePlots(player1, player2)
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
