"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import math
import numpy
import sys

import thinkbayes
import myplot


def ReadScale(filename='sat_scale.csv', col=2):
    """Reads a CSV file of SAT scales (maps from raw score to standard score).

    Args:
      filename: string filename
      col: which column to start with (0=Reading, 2=Math, 4=Writing)

    Returns: thinkbayes.Interpolator object
    """
    def ParseRange(s):
        t = [int(x) for x in s.split('-')]
        return 1.0 * sum(t) / len(t)

    fp = open(filename)
    reader = csv.reader(fp)
    raws = []
    scores = []

    for t in reader:
        try:
            raw = int(t[col])
            raws.append(raw)
            score = ParseRange(t[col+1])
            scores.append(score)
        except:
            pass

    raws.sort()
    scores.sort()
    return thinkbayes.Interpolator(raws, scores)


def ReadRanks(filename='sat_ranks.csv'):
    """Reads a CSV file of SAT scores.

    Args:
      filename: string filename

    Returns:
      list of (score, freq) pairs
    """
    fp = open(filename)
    reader = csv.reader(fp)
    res = []

    for t in reader:
        try:
            score = int(t[0])
            freq = int(t[1])
            res.append((score, freq))
        except ValueError:
            pass

    return res


def DivideValues(pmf, denom):
    """Divides the values in a Pmf by denom.

    Returns a new Pmf.
    """
    new = thinkbayes.Pmf()
    denom = float(denom)
    for val, prob in pmf.Items():
        x = val / denom
        new.Set(x, prob)
    return new


class Exam(object):
    """Encapsulates information about an exam.

    Contains the distribution of scaled scores and an
    Interpolator that maps between scaled and raw scores.
    """
    def __init__(self):
        self.scale = ReadScale()

        scores = ReadRanks()
        hist = thinkbayes.MakeHistFromDict(dict(scores))
        score_dist = thinkbayes.MakePmfFromHist(hist)

        self.raw = self.ReverseScale(score_dist)
        self.max_score = max(self.raw.Values())
        self.prior = DivideValues(self.raw, denom=self.max_score)
        
        center = -0.05
        width = 1.8
        self.difficulties = MakeDifficulties(center, width, self.max_score)

    def CompareScores(self, a_score, b_score):
        """Computes posteriors for two test scores and the likelihood ratio.

        a_score, b_score: scales SAT scores
        """
        a_sat = constructor(self, a_score)
        b_sat = constructor(self, b_score)

        a_sat.PlotPosteriors(b_sat)

        top = TopLevel(['Alice', 'Bob'])
        top.Update((a_score, b_score))
        top.Print()

        ratio = top.Prob('Alice') / top.Prob('Bob')
        
        print 'Likelihood ratio', ratio

        posterior = ratio / (ratio + 1)
        print 'Posterior', posterior

    def MakeRawScoreDist(self):
        """Makes the distribution of raw scores for given difficulty.

        Assumes that efficacy is Gaussian(0, 1)
        """
        efficacies = thinkbayes.MakeGaussianPmf(0, 1.5, 3)

        #diff_pmf = thinkbayes.MakePmfFromList(self.difficulties)
        #myplot.Pmf(diff_pmf)
        #myplot.Pmf(efficacies)
        #myplot.Show()

        pmfs = thinkbayes.Pmf()
        for efficacy, prob in efficacies.Items():
            scores = self.PmfCorrect(efficacy)
            pmfs.Set(scores, prob)

        mix = thinkbayes.MakeMixture(pmfs)
        return mix

    def CalibrateDifficulty(self):
        """Make a plot showing the model distribution of raw scores."""
        myplot.Clf()

        cdf = thinkbayes.MakeCdfFromPmf(self.raw, name='data')
        myplot.Cdf(cdf)

        pmf = self.MakeRawScoreDist()
        cdf = thinkbayes.MakeCdfFromPmf(pmf, name='model')
        myplot.Cdf(cdf)
        
        myplot.Save(root='sat_calibrate',
                    xlabel='raw score',
                    ylabel='CDF',
                    formats=['pdf', 'eps'])

    def PmfCorrect(self, efficacy):
        pmf = PmfCorrect(efficacy, self.difficulties)
        return pmf

    def Lookup(self, raw):
        """Looks up a raw score and returns a scaled score."""
        return self.scale.Lookup(raw)
        
    def Reverse(self, score):
        """Looks up a scaled score and returns a raw score.

        Since we ignore the penalty, negative scores round up to zero.
        """
        raw = self.scale.Reverse(score)
        return raw if raw > 0 else 0
        
    def ReverseScale(self, pmf):
        """Applies the reverse scale to the values of a PMF.

        Args:
            pmf: Pmf object
            scale: Interpolator object

        Returns:
            new Pmf
        """
        new = thinkbayes.Pmf()
        for val, prob in pmf.Items():
            raw = self.Reverse(val)
            new.Incr(raw, prob)
        return new



class Sat(thinkbayes.Suite):
    """Represents the distribution of efficacy for a test-taker."""

    def __init__(self, exam, score):
        thinkbayes.Suite.__init__(self)

        self.exam = exam
        self.score = score

        # start with the prior distribution
        for p_correct, prob in exam.prior.Items():
            self.Set(p_correct, prob)

        # update based on an exam score
        self.Update(score)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of a test score, given efficacy."""
        p_correct = hypo
        score = data
        raw = self.exam.Reverse(score)

        yes, no = raw, self.exam.max_score - raw
        like = thinkbayes.EvalBinomialPmf(p_correct, yes, no)
        return like

    def PlotPosteriors(self, other):
        """Plots posterior distributions of efficacy.

        self, other: Sat objects.
        """
        cdf1 = thinkbayes.MakeCdfFromPmf(self, 'posterior %d' % self.score)
        cdf2 = thinkbayes.MakeCdfFromPmf(other, 'posterior %d' % other.score)

        myplot.Clf()
        myplot.Cdfs([cdf1, cdf2])
        myplot.Save(xlabel='p_correct', 
                    ylabel='CDF', 
                    axis=[0.7, 1.0, 0.0, 1.0],
                    root='sat_posteriors_p_corr',
                    formats=['pdf', 'eps'])


class Sat2(Sat):

    def __init__(self, exam, score):
        thinkbayes.Suite.__init__(self)

        self.exam = exam
        self.score = score

        # start with the prior distribution
        efficacies = thinkbayes.MakeGaussianPmf(0, 1.5, 3)
        for efficacy, prob in efficacies.Items():
            self.Set(efficacy, prob)

        # update based on an exam score
        self.Update(score)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of a test score, given efficacy."""
        efficacy = hypo
        score = data
        raw = self.exam.Reverse(score)

        pmf = self.exam.PmfCorrect(efficacy)
        like = pmf.Prob(raw)
        return like
    
    def PlotPosteriors(self, other):
        """Plots posterior distributions of efficacy.

        self, other: Sat objects.
        """
        cdf1 = thinkbayes.MakeCdfFromPmf(self, 'posterior %d' % self.score)
        cdf2 = thinkbayes.MakeCdfFromPmf(other, 'posterior %d' % other.score)

        myplot.Clf()
        myplot.Cdfs([cdf1, cdf2])
        myplot.Save(xlabel='efficacy', 
                    ylabel='CDF', 
                    axis=[0, 4.6, 0.0, 1.0],
                    root='sat_posteriors_eff',
                    formats=['pdf', 'eps'])


def PlotPriorDist(pmf):
    """Plot the prior distribution of p_correct."""
    myplot.Clf()
    cdf1 = thinkbayes.MakeCdfFromPmf(pmf, 'prior')
    myplot.Cdf(cdf1)
    myplot.Save(root='sat_prior',
                xlabel='p_correct', 
                ylabel='CDF',
                formats=['pdf', 'eps'])


def PlotRawDist(raw):
    cdf2 = thinkbayes.MakeCdfFromPmf(raw, 'raw')
    myplot.Cdf(cdf2)
    myplot.Show(xlabel='score', 
               ylabel='CDF')


class TopLevel(thinkbayes.Suite):

    def Update(self, data):
        a_score, b_score = data

        exam = Exam()
        a_sat = constructor(exam, a_score)
        b_sat = constructor(exam, b_score)

        a_like = thinkbayes.PmfProbGreater(a_sat, b_sat)
        b_like = thinkbayes.PmfProbGreater(b_sat, a_sat)
        
        self.Mult('Alice', a_like)
        self.Mult('Bob', b_like)

        # self.Normalize()


def ProbCorrect(efficacy, difficulty, a=1.0):
    """Returns the probability that a person gets a question right.

    efficacy: personal ability to answer questions
    difficulty: how hard the question is

    Returns: float prob
    """
    return 1 / (1 + math.exp(-a * (efficacy - difficulty)))


def BinaryPmf(p):
    """Makes a Pmf with values 1 and 0.
    
    p: probability given to 1
    
    Returns: Pmf object
    """
    pmf = thinkbayes.Pmf()
    pmf.Set(1, p)
    pmf.Set(0, 1-p)
    return pmf


def PmfCorrect(efficacy, difficulties):
    """Computes the distribution of correct responses.

    efficacy: personal ability to answer questions
    difficulties: list of difficulties, one for each question

    Returns: new Pmf object
    """
    pmf0 = thinkbayes.Suite([0])

    ps = [ProbCorrect(efficacy, difficulty) for difficulty in difficulties]
    pmfs = [BinaryPmf(p) for p in ps]
    dist = sum(pmfs, pmf0)
    return dist


def MakeDifficulties(center, width, n):
    """Makes a list of n difficulties with a given center and width.

    Returns: list of n floats between center-width and center+width
    """
    low, high = center-width, center+width
    return numpy.linspace(low, high, n)


def TestEfficacy():
    efficacy = 0
    difficulties = [0] * 54
    dist = PmfCorrect(efficacy, difficulties)
    dist.name = 'constant'
    myplot.Pmf(dist)

    difficulties = MakeDifficulties(-1, 1, 54)
    dist = PmfCorrect(efficacy, difficulties)
    dist.name = 'uniform'
    myplot.Pmf(dist)

    myplot.Show()


# which version of the Sat class to use, Sat or Sat2
constructor = Sat2

def main(script):
    global constructor

    exam = Exam()

    PlotPriorDist(exam.prior)
    exam.CalibrateDifficulty()

    constructor = Sat
    exam.CompareScores(780, 740)

    constructor = Sat2
    exam.CompareScores(780, 740)


if __name__ == '__main__':
    main(*sys.argv)
