"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import math
import myplot
import numpy
import thinkbayes

import matplotlib.pyplot as pyplot

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
            # print heading, data[0], len(data)
            res.append(data)
        except ValueError:
            pass

    fp.close()
    return zip(*res)
    

class Price(thinkbayes.Suite):

    def __init__(self, pmf, player, name=''):
        """Constructs the suite.

        pmf: prior distribution of value
        player: Player object
        """
        thinkbayes.Suite.__init__(self, name=name)

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

        error = actual_price - my_guess
        like = self.player.ErrorDensity(error)

        return like


class GainCalculator(object):

    def __init__(self, player, opponent):
        """Constructs the calculator.

        player: Player
        opponent: Player
        """
        self.player = player
        self.opponent = opponent

    def ExpectedGains(self, low=0, high=75000, n=101):
        bids = numpy.linspace(low, high, n)

        gains = [self.ExpectedGain(bid) for bid in bids]

        #for bid, ret in zip(bids, gains):
        #    print bid, ret

        return bids, gains

    def ExpectedGain(self, bid):
        """Computes the expected return of a given bid.

        bid: your bid
        suite: posterior distribution of prices
        """
        suite = self.player.posterior
        total = 0
        for price, prob in sorted(suite.Items()):
            roi = self.Gain(bid, price)
            total += prob * roi
        return total

    def Gain(self, bid, price):
        """Computes the return of a bid, given the actual price.
        """
        # if you overbid, you get nothing
        if bid > price:
            return 0

        # otherwise compute the probability of winning
        diff = price - bid
        prob = self.ProbWin(diff)

        # if you are within 250 dollars, you win both showcases
        if diff <= 250:
            return 2 * price * prob
        else:
            return price * prob

    def ProbWin(self, diff):
        """Computes the probability of winning for a given diff.

        diff: how much your bid was off by
        """
        prob = (self.opponent.ProbOverbid() + 
                self.opponent.ProbWorseThan(diff))
        return prob


def RemoveNegatives(pmf):
    """Remove negative values from the Pmf."""
    for val in pmf.Values():
        if val < 0:
            pmf.Remove(val)
    pmf.Normalize()


class Player(object):

    low, high = 0, 75000
    n = 101
    xs = numpy.linspace(low, high, n)
    diffs = numpy.linspace(-30000, 30000, n)

    def __init__(self, val, bid, diff):
        """Construct the Player.

        val: sequence of values
        bid: sequence of bids
        diff: sequence of underness (negative means over)
        """
        self.val = val
        self.big = bid
        self.diff = diff

        self.cdf_val = thinkbayes.MakeCdfFromList(val)
        self.cdf_diff = thinkbayes.MakeCdfFromList(diff)

        self.kde_val = thinkbayes.EstimatedPdf(val)

        mu = 0
        sigma = numpy.std(self.diff)
        self.pdf_error = thinkbayes.GaussianPdf(mu, sigma)

    def ErrorDensity(self, error):
        """Density of the given error in the distribution of error.

        error: how much the bid is under the actual price
        """
        return self.pdf_error.Density(error)

    def PmfVal(self):
        """Returns a new Pmf of values.

        A discrete version of the estimated Pdf.
        """
        return self.kde_val.MakePmf(self.xs)

    def CdfDiff(self):
        """Returns a reference to the Cdf of differences (underness).
        """
        return self.cdf_diff

    def ProbOverbid(self):
        """Returns the probability this player overbids.
        """
        return self.cdf_diff.Prob(-1)

    def ProbWorseThan(self, diff):
        """Probability this player's diff is greater than the given diff.

        diff: how much the oppenent is off by (always positive)
        """
        return 1 - self.cdf_diff.Prob(diff)

    def MakeBeliefs(self, estimate):
        """Makes a posterior distribution based on estimated value.

        Sets attributes prior and posterior.

        estimate: what the player thinks the showcase is worth        
        """
        pmf = self.PmfVal()
        self.prior = Price(pmf, self, name='prior')
        self.posterior = self.prior.Copy(name='posterior')
        self.posterior.Update(estimate)

    def OptimalBid(self, estimate, opponent):
        """Computes the bid that maximizes expected return.
        
        estimate: what the player thinks the showcase is worth 
        opponent: Player

        Returns: (optimal bid, expected gain)
        """
        self.MakeBeliefs(estimate)
        calc = GainCalculator(self, opponent)
        bids, gains = calc.ExpectedGains()
        gain, bid = max(zip(gains, bids))
        return bid, gain

    def PlotBeliefs(self, root):
        """Plots prior and posterior beliefs.

        root: string filename root for saved figure
        """
        myplot.Clf()
        myplot.PrePlot(num=2)
        myplot.Pmfs([self.prior, self.posterior])
        myplot.Save(root=root,
                    xlabel='price ($)',
                    formats=FORMATS)


def MakePlots(player1, player2):
    """Generates two plots.

    price1 shows the priors for the two players
    price2 shows the distribution of diff for the two players
    """

    # plot the prior distribution of value for both players
    myplot.Clf()
    myplot.PrePlot(num=2)
    pmf1 = player1.PmfVal()
    pmf1.name = 'showcase 1'
    pmf2 = player2.PmfVal()
    pmf2.name = 'showcase 2'
    myplot.Pmfs([pmf1, pmf2])
    myplot.Save(root='price1',
                xlabel='price ($)',
                formats=FORMATS)

    # plot the historical distribution of underness for both players
    myplot.Clf()
    myplot.PrePlot(num=2)
    cdf1 = player1.CdfDiff()
    cdf1.name = 'player 1'
    cdf2 = player2.CdfDiff()
    cdf2.name = 'player 2'

    print 'Player 1 overbids', player1.ProbOverbid()
    print 'Player 2 overbids', player2.ProbOverbid()

    myplot.Cdfs([cdf1, cdf2])
    myplot.Save(root='price2',
                xlabel='underness ($)',
                formats=FORMATS)




def MakePlayers():
    """Reads data and makes player objects."""
    data = ReadData(filename='showcases.2011.csv')
    data += ReadData(filename='showcases.2012.csv')

    cols = zip(*data)
    val1, val2, bid1, bid2, diff1, diff2 = cols

    player1 = Player(val1, bid1, diff1)
    player2 = Player(val2, bid2, diff2)

    return player1, player2


def PlotExpectedGains(estimate1=20000, estimate2=40000):
    """Plots expected gains as a function of bid.

    estimate1: player1's estimate of the value of showcase 1
    estimate2: player2's estimate of the value of showcase 2
    """
    player1, player2 = MakePlayers()
    MakePlots(player1, player2)

    player1.MakeBeliefs(estimate1)
    player2.MakeBeliefs(estimate2)

    print 'Player 1 mean', player1.posterior.Mean()
    print 'Player 2 mean', player2.posterior.Mean()
    print 'Player 1 mle', player1.posterior.MaximumLikelihood()
    print 'Player 2 mle', player2.posterior.MaximumLikelihood()

    player1.PlotBeliefs('price3')
    player2.PlotBeliefs('price4')

    calc1 = GainCalculator(player1, player2)
    calc2 = GainCalculator(player2, player1)

    myplot.Clf()
    myplot.PrePlot(num=2)

    bids, gains = calc1.ExpectedGains()
    myplot.Plot(bids, gains, label='Player 1')
    print 'Player 1 optimal bid', max(zip(gains, bids))

    bids, gains = calc2.ExpectedGains()
    myplot.Plot(bids, gains, label='Player 2')
    print 'Player 2 optimal bid', max(zip(gains, bids))

    myplot.Save(root='price5',
                xlabel='price ($)',
                formats=FORMATS)


def PlotOptimalBid():
    """Plots optimal bid vs estimated value.
    """
    player1, player2 = MakePlayers()
    estimates = numpy.linspace(15000, 60000, 21)

    res = []
    for estimate in estimates:
        player1.MakeBeliefs(estimate)

        mean = player1.posterior.Mean()
        mle = player1.posterior.MaximumLikelihood()

        calc = GainCalculator(player1, player2)
        bids, gains = calc.ExpectedGains()
        gain, bid = max(zip(gains, bids))

        res.append((estimate, mean, mle, gain, bid))

    estimates, means, mles, gains, bids = zip(*res)
    
    myplot.PrePlot(num=3)
    pyplot.plot([15000,60000], [15000,60000], color='gray')
    myplot.Plot(estimates, means, label='mean')
    #myplot.Plot(estimates, mles, label='MLE')
    myplot.Plot(estimates, bids, label='bid')
    myplot.Plot(estimates, gains, label='gain')
    myplot.Save(root='price6',
                xlabel='estimated price ($)',
                formats=FORMATS)


def TestCode():
    """Check some intermediate results."""
    # test ProbWin
    for diff in [0, 100, 1000, 10000, 20000]:
        print diff, calc1.ProbWin(diff)
    print

    # test Return
    price = 20000
    for bid in [17000, 18000, 19000, 19500, 19800, 20001]:
        print bid, calc1.Gain(bid, price)
    print


def main():
    PlotExpectedGains()
    PlotOptimalBid()



if __name__ == '__main__':
    main()
