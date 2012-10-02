"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import numpy

import thinkbayes
import myplot


def EvalGaussianPdf(mu, sigma, x):
    """Computes the unnormalized PDF of the normal distribution.

    mu: mean
    sigma: standard deviation
    x: value
    
    returns: float probability density (unnormalized)
    """
    z = (x - mu) / sigma
    p = math.exp(-z**2/2)
    return p


def MakeGaussianPmf(mu, sigma, num_sigmas, n=201):
    """Makes a PMF discrete approx to a Gaussian distribution.
    
    mu: float mean
    sigma: float standard deviation
    num_sigmas: how many sigmas to extend in each direction
    n: number of values in the Pmf

    returns: normalized Pmf
    """
    pmf = thinkbayes.Pmf()
    low = mu - num_sigmas*sigma
    high = mu + num_sigmas*sigma

    for x in numpy.linspace(low, high, n):
        p = EvalGaussianPdf(mu, sigma, x)
        pmf.Set(x, p)
    pmf.Normalize()
    return pmf


def EvalPoissonPmf(lam, k):
    """Computes the Poisson PMF.

    lam: parameter lambda in events per unit time
    k: number of events

    returns: float probability
    """
    return (lam)**k * math.exp(lam) / math.factorial(k)


def MakePoissonPmf(lam, high):
    """Makes a PMF discrete approx to a Poisson distribution.

    lam: parameter lambda in events per unit time
    high: upper bound of the Pmf

    returns: normalized Pmf
    """
    pmf = thinkbayes.Pmf()
    for k in xrange(0, high+1):
        p = EvalPoissonPmf(lam, k)
        pmf.Set(k, p)
    pmf.Normalize()
    return pmf


def EvalExponentialPdf(lam, x):
    """Computes the exponential PDF.

    lam: parameter lambda in events per unit time
    x: value

    returns: float probability density
    """
    return lam * math.exp(-lam * x)


def MakeExponentialPmf(lam, high, n=200):
    """Makes a PMF discrete approx to an exponential distribution.

    lam: parameter lambda in events per unit time
    high: upper bound
    n: number of values in the Pmf

    returns: normalized Pmf
    """
    pmf = thinkbayes.Pmf()
    for x in numpy.linspace(0, high, n):
        p = EvalExponentialPdf(lam, x)
        pmf.Set(x, p)
    pmf.Normalize()
    return pmf


class Hockey(thinkbayes.Suite):
    def __init__(self):
        thinkbayes.Suite.__init__(self)

        pmf = MakeGaussianPmf(2.7, 0.3, 6)
        for x, p in pmf.Items():
            self.Set(x, p)
            
    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        Evaluates the Poisson PMF for lambda, k and t.

        I dropped the k! term in the denominator because it does not
        depend on lambda.

        hypo: goal scoring rate in goals per game
        data: goals scored in one period
        """
        t = 1.0
        lam = hypo
        k = data
        like = (lam*t)**k * math.exp(lam*t)
        return like


def MakeGoalPmf(suite):
    """Makes the distribution of goals scored, given distribution of lam.

    suite: distribution of goal-scoring rate

    returns: Pmf of goals per game
    """
    metapmf = thinkbayes.Pmf()
    high = 10

    for lam, prob in suite.Items():
        pmf = MakePoissonPmf(lam, high)
        metapmf.Set(pmf, prob)

    mix = thinkbayes.MakeMixture(metapmf, name=suite.name)
    return mix


def MakeGoalTimePmf(suite):
    """Makes the distribution of time til first goal.

    suite: distribution of goal-scoring rate

    returns: Pmf of goals per game
    """
    metapmf = thinkbayes.Pmf()

    for lam, prob in suite.Items():
        pmf = MakeExponentialPmf(lam, high=2, n=2001)
        metapmf.Set(pmf, prob)

    mix = thinkbayes.MakeMixture(metapmf, name=suite.name)
    return mix


def main():
    suite1 = Hockey()
    suite1.name = 'bruins'
    suite1.UpdateSet([0, 2, 8, 4])
    goal_dist1 = MakeGoalPmf(suite1)
    time_dist1 = MakeGoalTimePmf(suite1)
    
    suite2 = Hockey()
    suite2.name = 'sabres'
    suite2.UpdateSet([1, 3, 1, 0])
    goal_dist2 = MakeGoalPmf(suite2)
    time_dist2 = MakeGoalTimePmf(suite2)
    
    myplot.Clf()
    myplot.Pmf(suite1)
    myplot.Pmf(suite2)
    myplot.Save(root='hockey1',
                xlabel='Goals per game',
                ylabel='Probability',
                formats=['pdf', 'eps'])

    myplot.Clf()
    myplot.Pmf(goal_dist1)
    myplot.Pmf(goal_dist2)
    myplot.Save(root='hockey2',
                xlabel='Goals',
                ylabel='Probability',
                formats=['pdf', 'eps'])

    myplot.Clf()
    myplot.Pmf(time_dist1)
    myplot.Pmf(time_dist2)    
    myplot.Save(root='hockey3',
                xlabel='Games until goal',
                ylabel='Probability',
                formats=['pdf', 'eps'])

    diff = goal_dist1 - goal_dist2
    p_win = diff.ProbGreater(0)
    p_loss = diff.ProbLess(0)
    p_tie = diff.Prob(0)

    print p_win, p_loss, p_tie

    p_overtime = thinkbayes.PmfProbLess(time_dist1, time_dist2)
    p_adjust = thinkbayes.PmfProbEqual(time_dist1, time_dist2)
    p_overtime += p_adjust / 2
    print p_overtime 

    print p_overtime * p_tie
    p_win += p_overtime * p_tie
    print 'p_win', p_win

    # win the next two
    p_series = p_win**2

    # split the next two, win the third
    p_series += 2 * p_win * (1-p_win) * p_win

    print 'p_series', p_series


if __name__ == '__main__':
    main()
