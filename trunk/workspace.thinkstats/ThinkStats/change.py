"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import thinkbayes

import matplotlib.pyplot as pyplot
import myplot
import numpy

import csv
import math
import random
import sys
import time

from math import log

formats = ['pdf', 'eps', 'png']


class Series(object):
    """Represents a time series of values."""

    def __init__(self, t):
        self.n = len(t)
        self.t = t
        self.cumsum = numpy.cumsum(t)

    def Split(self, n):
        """Returns the number and sum of elements on each half of a split.

        n: how many elements in the first half

        Returns: (n, sum of first n elements, m, sum of remaining m elements)
        """
        assert n > 0
        m = self.n - n
        s1 = self.cumsum[n-1]
        s2 = self.cumsum[-1] - s1
        return n, s1, m, s2


class Expo(thinkbayes.Suite):
    """Suite of hypotheses about parameter of an exponential distribution."""

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under a hypothesis.

        hypo: hypothetical value of the parameter
        data: tuple of (number of elements, their sum)

        Returns: float likelihood
        """
        like = math.exp(self.LogLikelihood(hypo, data))
        return like

    def LogLikelihood(self, hypo, data):
        """Computes the log-likelihood of the data under a hypothesis.

        hypo: hypothetical value of the parameter
        data: tuple of (number of elements, their sum)

        Returns: float log-likelihood
        """
        lam = hypo
        n, s = data
        return n * log(lam) - lam * s
            

class Erlang(thinkbayes.Suite):
    """Suite of hypotheses about the parameter of an Erlang distribution."""

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under a hypothesis.

        hypo: hypothetical value of the parameter
        data: tuple of (number of elements, their sum)

        Returns: float likelihood
        """
        like = math.exp(self.LogLikelihood(hypo, data))
        return like

    def LogLikelihood(self, hypo, data):
        """Computes the log-likelihood of the data under a hypothesis.

        hypo: hypothetical value of the parameter
        data: tuple of (number of elements, their sum)

        Returns: float log-likelihood
        """
        lam = hypo
        n, s = data

        if n == 1:
            return n * log(lam) - lam * s
        else:
            logn = log(n-1)
            return n * log(lam) - lam * s + (n-1) * (log(s) - logn + 1)
            


class Split(thinkbayes.Suite):
    """Suite of hypotheses about the location of a changepoint."""

    def __init__(self, n, low, high, num_lams=101):
        """Initializes the suite.

        n: number of elements in the series
        low, high: range of values for lam
        num_lams: number of values for lam
        """
        hypos = range(1, n)
        thinkbayes.Suite.__init__(self, hypos)
        self.lams = numpy.linspace(low, high, num_lams)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under a hypothesis.

        hypo: hypothetical location of the changepoint
        data: Series object

        Returns: float likelihood
        """
        n = hypo
        series = data
        n, s1, m, s2 = series.Split(n)

        like = self.HalfLike(n, s1) * self.HalfLike(m, s2)
        return like

    def HalfLike(self, n, s):
        """Computes the likelihood of half of a split.

        n: number of elements
        s: sum of the elements

        Returns: float likelihood
        """
        expo = Expo(self.lams)
        like = expo.Update((n, s))
        return like


class Split2(Split):
    """Suite of hypotheses about the location of a changepoint."""

    def Update(self, data):
        """Updates the suite based on the data.

        data: Series object
        """
        series = data
        splits = [series.Split(n) for n in self.Values()]
        ns, s1s, ms, s2s = zip(*splits)

        likes1 = self.HalfLike(ns, s1s)
        likes2 = self.HalfLike(ms, s2s)
        
        scale, likes = ShiftExp(numpy.log(likes1) + numpy.log(likes2))
        #likes = likes1 * likes2

        for hypo, like in zip(self.Values(), likes):
            self.Mult(hypo, like)
        return self.Normalize()

    def HalfLike(self, ns, ss):
        """Computes a vector of likelihoods for half of a split.

        ns: vector of number of elements
        ss: vector of sums of elements

        Returns: numpy vector of likelihoods, one for each value of n
        """
        loglams = numpy.log(self.lams)
        logarray = numpy.outer(ns, loglams) - numpy.outer(ss, self.lams)
        scale, array = ShiftExp(logarray)
        likes = array.sum(axis=1)
        return likes


def OuterSum(a, b):
    """Computes the outer sum of two vectors.

    a, b: numpy vectors

    Returns: matrix with the height of a, width of b.
    """
    return a[:, numpy.newaxis] + b


class Split3(Split):
    """Suite of hypotheses about the location of a changepoint."""

    def Update(self, data):
        """Updates the suite based on the data.

        data: Series object
        """
        series = data
        splits = [series.Split(n) for n in self.Values()]
        ns, s1s, ms, s2s = zip(*splits)

        logs1 = self.HalfLike(ns, s1s)
        logs2 = self.HalfLike(ms, s2s)

        num = len(ns)
        scales = numpy.zeros(num)
        likes = numpy.zeros(num)
        
        for i in range(num):
            terms = OuterSum(logs1[i], logs2[i])
            scales[i], array = ShiftExp(terms)
            likes[i] = array.sum()

        _, scales = ShiftExp(scales)
        likes /= scales

        for hypo, like in zip(self.Values(), likes):
            self.Mult(hypo, like)

        return self.Normalize()

    def HalfLike(self, ns, ss):
        """Computes a vector of likelihoods for half of a split.

        ns: vector of number of elements
        ss: vector of sums of elements

        Returns: numpy vector of likelihoods, one for each value of n
        """
        loglams = numpy.log(self.lams)
        logarray = numpy.outer(ns, loglams) - numpy.outer(ss, self.lams)
        return logarray


def ShiftExp(array):
    """Scales an array of logs and then exponentiates.

    Returns the log of the scale factor.

    array: numpy matrix

    Returns: scale, result matrix
    """
    scale = -array.max()
    array += scale
    return scale, numpy.exp(array)


def MakeSeries(n, lam1, m, lam2):
    t1 = [random.expovariate(lam1) for i in range(n)]
    t2 = [random.expovariate(lam2) for i in range(m)]

    t1.extend(t2)
    return t1


def main(script, *args):
    #random.seed(22)

    lam1 = 1
    lam2 = 1.5
    
    n = 300
    t = MakeSeries(n, lam1, n, lam2)
    series = Series(t)

    low, high = 0.01, 5.01

    for cons in [Split2, Split3]:
        print cons.__name__
        split = cons(series.n, low, high)
        split.Update(series)
        cdf = thinkbayes.MakeCdfFromPmf(split)
        cdf.name = cons.__name__
        myplot.Cdf(cdf)

    myplot.Show()


def ExpoErlangDemo():
    num = 10

    lam1 = 1
    lam2 = 2
    t = MakeSeries(num, lam1, num, lam2)
    series = Series(t)
    n, s1, m, s2 = series.Split(num)

    print n, s1, m, s2

    low, high = 0.01, 5.01
    lams = numpy.linspace(low, high, 101)

    expo = Expo(lams)
    expo.name = 'expo'
    expo.Update((n, s1))
    
    erlang = Erlang(lams)
    erlang.name = 'erlang'
    erlang.Update((n, s1))
    
    myplot.Pmf(expo)
    myplot.Pmf(erlang)
    myplot.Show()


def ExpoDemo():
    num = 10

    lam1 = 1
    lam2 = 2
    t = MakeSeries(num, lam1, num, lam2)
    series = Series(t)
    n, s1, m, s2 = series.Split(num)

    print n, s1, m, s2

    low, high = 0.01, 5.01
    lams = numpy.linspace(low, high, 101)

    expo = Expo(lams)
    expo.Update((n, s1))
    
    expo2 = Expo(lams)
    expo2.Update((m, s2))
    
    myplot.Pmf(expo)
    myplot.Pmf(expo2)
    myplot.Show()


def ErlangDemo():
    num = 10

    lam1 = 1
    lam2 = 2
    t = MakeSeries(num, lam1, num, lam2)
    series = Series(t)
    n, s1, m, s2 = series.Split(num)

    print n, s1, m, s2

    low, high = 0.01, 5.01
    lams = numpy.linspace(low, high, 101)

    erlang = Erlang(lams)
    erlang.Update((n, s1))
    
    erlang2 = Erlang(lams)
    erlang2.Update((m, s2))
    
    myplot.Pmf(erlang)
    myplot.Pmf(erlang2)
    myplot.Show()


if __name__ == '__main__':
    main(*sys.argv)
