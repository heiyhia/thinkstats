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
    """Divides the values in a PMF by denom.  Returns a new PMF."""
    new = thinkbayes.Pmf()
    denom = float(denom)
    for val, prob in pmf.Items():
        if val >= 0:
            x = val / denom
            new.Set(x, prob)
    return new


class Exam(object):
    """Encapsulates information about an exam.

    Contains the distribution of scaled scores and an
    Interpolator that maps between scaled and raw scores.
    """
    def __init__(self, diff_std=1.0):
        self.scale = ReadScale()

        scores = ReadRanks()
        hist = thinkbayes.MakeHistFromDict(dict(scores))
        score_dist = thinkbayes.MakePmfFromHist(hist)

        raw = self.ReverseScale(score_dist)
        self.max_score = max(raw.Values())
        self.prior = DivideValues(raw, denom=self.max_score)
        
        low, high = -diff_std * 2, diff_std * 2
        self.difficulties = MakeDifficulties(low, high, self.max_score)

    def MakeRawScoreDist(self, eff_std):
        efficacies = thinkbayes.MakeGaussianPmf(0, eff_std, 3)

        pmfs = thinkbayes.Pmf()
        for efficacy, prob in efficacies.Items():
            scores = self.PmfCorrect(efficacy)
            pmfs.Set(scores, prob)

        mix = thinkbayes.MakeMixture(pmfs)
        return mix

    def PmfCorrect(self, efficacy):
        pmf = PmfCorrect(efficacy, self.difficulties)
        return pmf

    def Lookup(self, raw):
        """Looks up a raw score and returns a scaled score."""
        return self.scale.Lookup(raw)
        
    def Reverse(self, score):
        """Looks up a scaled score and returns a raw score."""
        return self.scale.Reverse(score)
        
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
        for efficacy, prob in exam.prior.Items():
            self.Set(efficacy, prob)

        # update based on an exam score
        self.Update(score)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of a test score, given efficacy."""
        efficacy = hypo
        score = data

        raw = self.exam.Reverse(score)
        yes, no = raw, self.exam.max_score - raw

        like = thinkbayes.EvalBinomialPmf(efficacy, yes, no)
        return like

    def ScaledConfidenceInterval(self):
        """Computes the credible interval for someone with the given score."""
        low, high = thinkbayes.ConfidenceInterval(self, 90)
        scale_low = self.exam.Lookup(low * self.exam.max_score)
        scale_high = self.exam.Lookup(high * self.exam.max_score)
        return scale_low, scale_high


class Sat2(Sat):
    def Likelihood(self, hypo, data):
        """Computes the likelihood of a test score, given efficacy."""
        efficacy = hypo
        score = data

        pmf = self.exam.PmfCorrect(efficacy)
        raw = self.exam.Reverse(score)
        like = pmf.Prob(raw)
        return like
    

def PlotPosteriors(sat1, sat2):
    """Plots posterior distributions of efficacy.

    sat1, sat2: Sat objects.
    """
    cdf1 = thinkbayes.MakeCdfFromPmf(sat1, 'posterior %d' % sat1.score)
    cdf2 = thinkbayes.MakeCdfFromPmf(sat2, 'posterior %d' % sat2.score)

    myplot.Cdfs([cdf1, cdf2])
    myplot.Save(xlabel='efficacy', 
                ylabel='CDF', 
                axis=[0.6, 1.0, 0.0, 1.0],
                root='sat_posteriors')


def PlotPriorDist(pmf):
    cdf1 = thinkbayes.MakeCdfFromPmf(pmf, 'prior')
    myplot.Cdf(cdf1)
    myplot.Save(root='sat_prior',
                xlabel='efficacy', 
                ylabel='CDF')


def PlotRawDist(raw):
    cdf2 = thinkbayes.MakeCdfFromPmf(raw, 'raw')
    myplot.Cdf(cdf2)
    myplot.Show(xlabel='score', 
               ylabel='CDF')


class TopLevel(thinkbayes.Suite):

    def Update(self, data):
        a_score, b_score = data

        exam = Exam()
        a_sat = Sat(exam, a_score)
        b_sat = Sat(exam, b_score)

        a_like = thinkbayes.PmfProbGreater(a_sat, b_sat)
        b_like = thinkbayes.PmfProbGreater(b_sat, a_sat)
        
        self.Mult('Alice', a_like)
        self.Mult('Bob', b_like)

        self.Normalize()


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


def MakeDifficulties(low, high, n):
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


def main(script):
    exam = Exam()
    #PlotPriorDist(exam.prior)

    cdf = thinkbayes.MakeCdfFromPmf(exam.prior)
    myplot.Cdf(cdf)

    pmf = exam.MakeRawScoreDist(1.6)
    pmf = DivideValues(pmf, 54)

    cdf = thinkbayes.MakeCdfFromPmf(pmf, name='mix')
    myplot.Cdf(cdf)

    myplot.Show()
    return

    a_score, b_score = 780, 740

    a_sat = Sat(exam, a_score)
    low_range = a_sat.ScaledConfidenceInterval()
    print a_score, low_range, (low_range[1] - low_range[0]) / 2.0

    b_sat = Sat(exam, b_score)
    high_range = b_sat.ScaledConfidenceInterval()
    print b_score, high_range, (high_range[1] - high_range[0]) / 2.0

    PlotPosteriors(a_sat, b_sat)

    top = TopLevel(['Alice', 'Bob'])
    top.Update((a_score, b_score))
    top.Print()


if __name__ == '__main__':
    main(*sys.argv)
