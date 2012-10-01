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


def EvalPoissonPmf(lam, t, k):
    """Computes the Poisson PMF.

    lam: parameter lambda in events per unit time
    t: duration in units of time
    k: number of events

    returns: float probability
    """
    return (lam*t)**k * math.exp(lam*t) / math.factorial(k)


def MakePoissonPmf(lam, t, high):
    """Makes a PMF discrete approx to a Poisson distribution.

    lam: parameter lambda in events per unit time
    t: duration in units of time
    high: upper bound of the Pmf

    returns: normalized Pmf
    """
    pmf = thinkbayes.Pmf()
    for k in xrange(0, high+1):
        p = EvalPoissonPmf(lam, t, k)
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

        pmf = MakeGaussianPmf(2.7, 0.3, 5)
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
    pmfs = thinkbayes.Pmf()
    t = 1.0
    high = 10

    for lam, prob in suite.Items():
        pmf = MakePoissonPmf(lam, t, high)
        pmfs.Set(pmf, prob)

    mix = thinkbayes.MakeMixture(pmfs, name=suite.name)
    return mix


def MakeGoalTimePmf(suite):
    """Makes the distribution of time til first goal.

    suite: distribution of goal-scoring rate

    returns: Pmf of goals per game
    """
    pmfs = thinkbayes.Pmf()
    t = 1.0
    high = 5   # max number of games without a goal

    for lam, prob in suite.Items():
        pmf = MakeExponentialPmf(lam, high, 501)
        pmfs.Set(pmf, prob)

    mix = thinkbayes.MakeMixture(pmfs, name=suite.name)
    return mix


def main():
    suite1 = Hockey()
    suite1.name = 'bruins'
    suite1.UpdateSet([5, 3, 1])
    goal_dist1 = MakeGoalPmf(suite1)
    time_dist1 = MakeGoalTimePmf(suite1)
    
    suite2 = Hockey()
    suite2.name = 'sabres'
    suite2.UpdateSet([1, 2, 3])
    goal_dist2 = MakeGoalPmf(suite2)
    time_dist2 = MakeGoalTimePmf(suite2)

    #myplot.Pmf(suite1)
    #myplot.Pmf(suite2)    
    myplot.Pmf(goal_dist1)
    myplot.Pmf(goal_dist2)    
    #myplot.Pmf(time_dist1)
    #myplot.Pmf(time_dist2)    
    myplot.Show()

    diff = goal_dist1 - goal_dist2
    p_win, p_loss, p_tie = diff.ProbGreater(0), diff.Prob(0), diff.ProbLess(0)

    print p_win, p_loss, p_tie

    p_overtime = thinkbayes.PmfProbLess(time_dist1, time_dist2)
    p_tie = thinkbayes.PmfProbEqual(time_dist1, time_dist2)
    print p_overtime 
    print p_overtime + p_tie


if __name__ == '__main__':
    main()
