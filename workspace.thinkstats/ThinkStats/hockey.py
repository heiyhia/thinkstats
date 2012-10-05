"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math

import columns
import thinkbayes
import myplot


class Hockey(thinkbayes.Suite):
    def __init__(self, name=''):
        thinkbayes.Suite.__init__(self, name=name)

        mu = 2.7
        sigma = 0.3
        sigma = 0.9
        pmf = thinkbayes.MakeGaussianPmf(mu, sigma, 6)
        for x, p in pmf.Items():
            self.Set(x, p)
            
    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        Evaluates the Poisson PMF for lambda and k.

        hypo: goal scoring rate in goals per game
        data: goals scored in one period
        """
        lam = hypo
        k = data
        like = thinkbayes.EvalPoissonPmf(lam, k)
        return like


def MakeGoalPmf(suite):
    """Makes the distribution of goals scored, given distribution of lam.

    suite: distribution of goal-scoring rate

    returns: Pmf of goals per game
    """
    metapmf = thinkbayes.Pmf()
    high = 10

    for lam, prob in suite.Items():
        pmf = thinkbayes.MakePoissonPmf(lam, high)
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
        pmf = thinkbayes.MakeExponentialPmf(lam, high=2, n=2001)
        metapmf.Set(pmf, prob)

    mix = thinkbayes.MakeMixture(metapmf, name=suite.name)
    return mix


class Game(object):
    """Represents a game."""
    convert = dict()

    def clean(self):
        self.goals = self.pd1 + self.pd2 + self.pd3


def ReadHockeyData(filename='hockey_data.csv'):
    game_list = columns.read_csv(filename, Game)
    games = {}
    
    for game in game_list:
        if game.season != 2012:
            continue
        key = game.game
        games.setdefault(key, []).append(game)

    return games


def ProcessHockeyData(games):
    pairs = {}

    for key, pair in games.iteritems():
        t1, t2 = pair
        key = t1.team, t2.team
        entry = t1.total, t2.total
        print key, entry
        pairs.setdefault(key, []).append(entry)


def main():
    games = ReadHockeyData()
    ProcessHockeyData(games)
    return

    suite1 = Hockey('bruins')
    suite1.UpdateSet([0, 2, 8, 4])
    goal_dist1 = MakeGoalPmf(suite1)
    time_dist1 = MakeGoalTimePmf(suite1)
    
    suite2 = Hockey('canucks')
    suite2.UpdateSet([1, 3, 1, 0])
    goal_dist2 = MakeGoalPmf(suite2)
    time_dist2 = MakeGoalTimePmf(suite2)
 
    print 'MLE bruins', thinkbayes.MaximumLikelihood(suite1)
    print 'MLE canucks', thinkbayes.MaximumLikelihood(suite2)
   
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
    print 'p_overtime', p_overtime 

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
