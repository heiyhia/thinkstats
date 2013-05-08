"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import thinkbayes

import matplotlib.pyplot as pyplot
import thinkplot
import numpy

import csv
import math
import random
import sys
import time

from math import log

FORMATS = ['pdf']

"""
Notation guide:

z: time between trains
x: time since the last train
y: time until the next train

zb: distribution of z as seen by a random arrival

"""


# gap_times in seconds, collected by Kai Austin and Brendan Ritter
# using http://developer.mbta.com/Data/Red.txt

# actual time between trains in seconds

observed_gap_times = [
    428, 705, 407, 465, 433, 425, 204, 
    506, 143, 351, 450, 598, 464, 749
]

# pairs of observation time in seconds and number of passenger arrivals

observed_passenger_data = [
    (150, 3),
    (250, 10),
    (660, 15),
    (510, 13),
    (30, 1),
]

def BiasPmf(pmf, name='', invert=False):
    """Returns the Pmf with oversampling proportional to value.

    If pmf is the distribution of true values, the result is the
    distribution that would be seen if values are oversampled in
    proportion to their values; for example, if you ask students
    how big their classes are, large classes are oversampled in
    proportion to their size.

    If invert=True, computes in inverse operation; for example,
    unbiasing a sample collected from students.

    Args:
      pmf: Pmf object.
      name: string name for the new Pmf.
      invert: boolean

     Returns:
       Pmf object
    """
    new_pmf = pmf.Copy(name=name)

    for x, p in pmf.Items():
        if invert:
            new_pmf.Mult(x, 1.0/x)
        else:
            new_pmf.Mult(x, x)
        
    new_pmf.Normalize()
    return new_pmf


def UnbiasPmf(pmf, name=''):
    """Returns the Pmf with oversampling proportional to 1/value.

    Args:
      pmf: Pmf object.
      name: string name for the new Pmf.

     Returns:
       Pmf object
    """
    return BiasPmf(pmf, name, invert=True)


def MakeUniformPmf(low, high):
    """Make a uniform Pmf.

    low: lowest value (inclusive)
    high: highest value (inclusive)
    """
    pmf = thinkbayes.Pmf()
    for x in MakeRange(low=low, high=high):
        pmf.Set(x, 1)
    pmf.Normalize()
    return pmf    
    

def MakeRange(low=10, high=1300, skip=10):
    """Makes a range representing possible gap times in seconds.

    low: where to start
    high: where to end
    skip: how many to skip
    """
    return range(low, high+skip, skip)


class WaitTimeCalculator(object):
    """Encapsulates the forward inference process.

    Given the actual distribution of gap times (z),
    computes the distribution of gaps as seen by
    a random passenger (zb), which yields the distribution
    of wait times (y) and the distribution of elapsed times (x).
    """

    def __init__(self, pmf, inverse=False):
        """Constructor.

        pmf: Pmf of either z or zb
        inverse: boolean, true if pmf is zb, false if pmf is z
        """
        if inverse:
            self.pmf_zb = pmf
            self.pmf_z = UnbiasPmf(pmf, name="z")
        else:
            self.pmf_z = pmf
            self.pmf_zb = BiasPmf(pmf, name="zb")

        # distribution of wait time
        self.pmf_y = PmfOfWaitTime(self.pmf_zb)

        # the distribution of elapsed time is the same as the
        # distribution of wait time
        self.pmf_x = self.pmf_y

    def GenerateSampleWaitTimes(self, n):
        """Generates a random sample of wait times.

        n: sample size

        Returns: sequence of values
        """
        cdf_y = thinkbayes.MakeCdfFromPmf(self.pmf_y)
        sample = cdf_y.Sample(n)
        return sample

    def GenerateSampleGaps(self, n):
        """Generates a random sample of gaps seen by passengers.

        n: sample size

        Returns: sequence of values
        """
        cdf_zb = thinkbayes.MakeCdfFromPmf(self.pmf_zb)
        sample = cdf_zb.Sample(n)
        return sample

    def GenerateSamplePassengers(self, lam, n):
        """Generates a sample wait time and number of arrivals.

        lam: arrival rate in passengers per second
        n: number of samples

        Returns: list of (k1, y, k2) tuples
        k1: passengers there on arrival
        y: wait time
        k2: passengers arrived while waiting
        """
        zs = self.GenerateSampleGaps(n)
        xs, ys = self.SplitGaps(zs)

        res = []
        for x, y in zip(xs, ys):
            k1 = numpy.random.poisson(lam * x)
            k2 = numpy.random.poisson(lam * y)
            res.append((k1, y, k2))

        return res

    def SplitGaps(self, zs):
        """Splits zs into xs and ys.

        zs: sequence of gaps
        
        Returns: tuple of sequences (xs, ys)
        """
        xs = [random.uniform(0, z) for z in zs]
        ys = [z-x for z, x in zip(zs, xs)]
        return xs, ys

    def PlotPmfs(self):
        """Plots the computed Pmfs.
        """
        print 'Mean z', self.pmf_z.Mean()
        print 'Mean zb', self.pmf_zb.Mean()
        print 'Mean y', self.pmf_y.Mean()

        thinkplot.Pmf(self.pmf_z)
        thinkplot.Pmf(self.pmf_zb)
        thinkplot.Pmf(self.pmf_y)

    def MakePlot(self, root='redline2'):
        """Plots the computed CDFs.
        """
        cdf_z = self.pmf_z.MakeCdf()
        cdf_zb = self.pmf_zb.MakeCdf()
        cdf_y = self.pmf_y.MakeCdf()

        cdfs = ScaleCdfs([cdf_y, cdf_z, cdf_zb], 1.0/60)

        thinkplot.Clf()
        thinkplot.PrePlot(3)
        thinkplot.Cdfs(cdfs)
        thinkplot.Save(root=root,
                       xlabel='Time (min)',
                       ylabel='CDF',
                       formats=FORMATS)


def PmfOfWaitTime(pmf_zb):
    """Distribution of wait time.

    pmf_zb: dist of gap time as seen by a random observer

    Returns: dist of wait time (also dist of elapsed time)
    """
    metapmf = thinkbayes.Pmf()
    for gap, prob in pmf_zb.Items():
        uniform = MakeUniformPmf(0, gap)
        metapmf.Set(uniform, prob)

    pmf_y = thinkbayes.MakeMixture(metapmf, name="y")
    return pmf_y


def ScaleCdfs(cdfs, factor):
    [cdf.Scale(factor) for cdf in cdfs]
    return cdfs


class ElapsedTimeEstimator(object):
    """Uses the number of passengers to estimate time since last train."""

    def __init__(self, wtc, lam, num_passengers):
        """Constructor.

        pmf_x: expected distribution of elapsed time
        lam: arrival rate in passengers per second
        num_passengers: # passengers seen on the platform
        """
        # prior for elapsed time
        self.prior_x = Elapsed(wtc.pmf_x, 'prior x')

        # posterior of elapsed time (based on number of passengers)
        self.post_x = self.prior_x.Copy(name='posterior x')
        self.post_x.Update((lam, num_passengers))

        # predictive distribution of wait time
        self.pmf_y = PredictWaitTime(wtc.pmf_zb, self.post_x)

    def MakePlot(self, root='redline3'):

        # observed gaps
        cdf_prior_x = self.prior_x.MakeCdf()
        cdf_post_x = self.post_x.MakeCdf()
        cdf_y = self.pmf_y.MakeCdf()

        cdfs = ScaleCdfs([cdf_prior_x, cdf_post_x, cdf_y], 1.0/60)

        thinkplot.Clf()
        thinkplot.PrePlot(3)
        thinkplot.Cdfs(cdfs)
        thinkplot.Save(root=root,
                       xlabel='Time (min)',
                       ylabel='CDF',
                       formats=FORMATS)


class ArrivalRate(thinkbayes.Suite):
    """Represents the distribution of arrival rates (lambda)."""

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        Evaluates the Poisson PMF for lambda and k.

        hypo: arrival rate in passengers per second
        data: tuple of elapsed_time and number of passengers
        """
        lam = hypo
        x, k = data
        like = thinkbayes.EvalPoissonPmf(lam * x, k)
        return like


class ArrivalRateEstimator(object):
    """Estimates arrival rate based on passengers that arrive while waiting.
    """

    def __init__(self, passenger_data):
        """Constructor

        passenger_data: sequence of (k1, y, k2) pairs
        """
        self.passenger_data = passenger_data

        # range for lambda
        low, high = 0, 5
        n = 101
        hypos = numpy.linspace(low, high, n) / 60

        self.prior_lam = ArrivalRate(hypos, name='prior')

        self.post_lam = self.prior_lam.Copy(name='posterior')

        for k1, y, k2 in passenger_data:
            self.post_lam.Update((y, k2))

        print 'Mean posterior lambda', self.post_lam.Mean()

    def MakePlot(self, root='redline1'):
        thinkplot.Clf()

        # convert units to passengers per minute
        cdf = self.post_lam.MakeCdf()
        cdf.Scale(60)

        thinkplot.Cdfs([cdf])

        thinkplot.Save(root=root,
                       xlabel='Arrival rate (passengers / min)',
                       ylabel='CDF',
                       formats=FORMATS)
                       

class Elapsed(thinkbayes.Suite):
    """Represents the distribution of elapsed time (x)."""

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        Evaluates the Poisson PMF for lambda and k.

        hypo: elapsed time since the last train
        data: tuple of arrival rate and number of passengers
        """
        x = hypo
        lam, k = data
        like = thinkbayes.EvalPoissonPmf(lam * x, k)
        return like


def PredictWaitTime(pmf_zb, pmf_x):
    """Computes the distribution of wait times.

    Enumerate all pairs of zb from pmf_zb and x from pmf_x,
    and accumulate the distribution of y = z - x.

    pmf_zb: distribution of gaps seen by random observer
    pmf_x: distribution of elapsed time
    """
    pmf_y = pmf_zb - pmf_x
    pmf_y.name = "y"
    RemoveNegatives(pmf_y)
    return pmf_y


def RemoveNegatives(pmf):
    """Removes negative values from a PMF.

    pmf: Pmf
    """
    for val in pmf.Values():
        if val < 0:
            pmf.Remove(val)
    pmf.Normalize()


class Gaps(thinkbayes.Suite):
    """Represents the distribution of gap times,
    as updated by an observed waiting time."""

    def Likelihood(self, hypo, data):
        """The likelihood of the data under the hypothesis.

        If the actual gap time is z, what is the likelihood
        of waiting y seconds?

        hypo: actual time between trains
        data: observed wait time
        """
        z = hypo
        y = data
        if y > z:
            return 0
        return 1.0 / z


class GapDirichlet(thinkbayes.Dirichlet):
    """Represents the distribution of prevalences for each
    gap time."""

    def __init__(self, xs):
        """Constructor.

        xs: sequence of possible gap times
        """
        n = len(xs)
        thinkbayes.Dirichlet.__init__(self, n)
        self.xs = xs
        self.mean_zbs = []

    def PmfMeanZb(self):
        return thinkbayes.MakePmfFromList(self.mean_zbs)

    def Preload(self, data):
        """Adds pseudocounts to the parameters.

        data: sequence of pseudocounts
        """
        thinkbayes.Dirichlet.Update(self, data)

    def Update(self, data):
        """Computes the likelihood of the data.

        data: wait time observed by random arrival (y)

        Returns: float probability
        """
        k, y = data

        print k, y
        prior = self.PredictivePmf(self.xs)
        gaps = Gaps(prior)
        gaps.Update(y)
        probs = gaps.Probs(self.xs)

        self.params += numpy.array(probs)

    def GapProbs(self, wait_time):
        """Not currently working
        """
        xs = range(len(self.params))
        ps = self.Random()
        prior = thinkbayes.MakePmfFromItems(zip(xs, ps))
        gaps = Gaps(prior)
        gaps.Update(wait_time)
        return gaps.Probs(needs_xs)


class GapDirichlet2(GapDirichlet):
    """Represents the distribution of prevalences for each
    gap time."""

    def Update(self, data):
        """Computes the likelihood of the data.

        data: wait time observed by random arrival (y)

        Returns: float probability
        """
        k, y = data

        # get the current best guess for pmf_z
        pmf_zb = self.PredictivePmf(self.xs)

        # use it to compute prior pmf_x, pmf_y, pmf_z
        wtc = WaitTimeCalculator(pmf_zb, inverse=True)

        # use the observed passengers to estimate posterior pmf_x
        elapsed = ElapsedTimeEstimator(wtc,
                                       lam=0.0333,
                                       num_passengers=k)

        # use posterior_x and observed y to estimate observed z
        obs_zb = elapsed.post_x + Floor(y)
        probs = obs_zb.Probs(self.xs)

        mean_zb = obs_zb.Mean()
        self.mean_zbs.append(mean_zb)
        print k, y, mean_zb

        # use observed z to update beliefs about pmf_z
        self.params += numpy.array(probs)


class GapTimeEstimator(object):
    """Infers gap times using passenger data."""

    def __init__(self, xs, pcounts, passenger_data):
        self.xs = xs
        self.pcounts = pcounts
        self.passenger_data = passenger_data

        self.wait_times = [y for k1, y, k2 in passenger_data]
        self.pmf_y = thinkbayes.MakePmfFromList(self.wait_times, name="y")

        n = len(self.xs)
        dirichlet = GapDirichlet2(self.xs)
        dirichlet.params /= 1.0

        dirichlet.Preload(self.pcounts)
        dirichlet.params /= 20.0

        self.prior_zb = dirichlet.PredictivePmf(self.xs, name="prior zb")
        
        for k1, y, k2 in passenger_data:
            dirichlet.Update((k1, y))

        self.pmf_mean_zb = dirichlet.PmfMeanZb()

        self.post_zb = dirichlet.PredictivePmf(self.xs, name="post zb")
        self.post_z = UnbiasPmf(self.post_zb, name="post z")

    def PlotPmfs(self):
        print 'Mean y', self.pmf_y.Mean()
        print 'Mean z', self.pmf_z.Mean()
        print 'Mean zb', self.pmf_zb.Mean()

        thinkplot.Pmf(self.pmf_y)
        thinkplot.Pmf(self.pmf_z)
        thinkplot.Pmf(self.pmf_zb)

    def MakePlot(self):

        thinkplot.Cdf(self.pmf_y.MakeCdf())
        thinkplot.Cdf(self.prior_zb.MakeCdf())
        thinkplot.Cdf(self.post_zb.MakeCdf())
        thinkplot.Cdf(self.pmf_mean_zb.MakeCdf())
        thinkplot.Show()


def Floor(x, factor=10):
    """Rounds down to the nearest multiple of factor.

    When factor=10, all numbers from 10 to 19 get floored to 10.
    """
    return int(x/factor) * factor


def MakePcounts(xs, gap_times):
    """
    """
    # use the actual gap times to make pcounts
    vals = [Floor(t) for t in gap_times]
    hist = thinkbayes.MakeHistFromList(vals)
    pcounts = hist.Freqs(xs)
    return pcounts


def TestITE():
    random.seed(17)

    xs = [60, 120, 240]
    
    gap_times = [60,60,60,60,60,120,120,120,240,240]

    # distribution of gap time (z)
    pdf_z = thinkbayes.EstimatedPdf(gap_times)
    pmf_z = pdf_z.MakePmf(xs, name="z")

    wtc = WaitTimeCalculator(pmf_z, inverse=False)

    lam = 0.0333
    n = 100
    passenger_data = wtc.GenerateSamplePassengers(lam, n)

    pcounts = MakePcounts(xs, gap_times)
    pcounts = [0, 0, 0]
    print pcounts

    ite = GapTimeEstimator(xs, pcounts, passenger_data)

    thinkplot.Clf()

    # thinkplot.Cdf(wtc.pmf_z.MakeCdf(name="actual z"))    
    thinkplot.Cdf(wtc.pmf_zb.MakeCdf(name="actual zb"))
    ite.MakePlot()


def GenerateFakeData(lam=0.0333, n=10):
    xs = MakeRange(low=10)
    pdf_z = thinkbayes.EstimatedPdf(observed_gap_times)
    pmf_z = pdf_z.MakePmf(xs, name="z")

    wtc = WaitTimeCalculator(pmf_z, inverse=False)
    passenger_data = wtc.GenerateSamplePassengers(lam, n)
    return wtc, passenger_data


def RunSimpleProcess(gap_times):
    xs = MakeRange(low=10)
    pdf_z = thinkbayes.EstimatedPdf(gap_times)
    pmf_z = pdf_z.MakePmf(xs, name="z")

    wtc = WaitTimeCalculator(pmf_z, inverse=False)    
    wtc.MakePlot()

    elapsed = ElapsedTimeEstimator(wtc,
                                   lam=0.0333,
                                   num_passengers=15)
    elapsed.MakePlot()


def main(script):
    RunSimpleProcess(observed_gap_times)
    return


    random.seed(17)
    wtc, passenger_data = GenerateFakeData(lam=0.0333, n=20)

    #TestITE()
    #return

    # wtc.MakePlot()

    xs = MakeRange(low=10)

    cdf_zb = wtc.pmf_zb.MakeCdf(name="actual zb")
    thinkplot.Cdf(cdf_zb)
    
    pcounts = MakePcounts(xs, observed_gap_times)
    pcounts = [0] * len(xs)
    print pcounts

    ite = GapTimeEstimator(xs, pcounts, passenger_data)
    ite.MakePlot()

    return

    are = ArrivalRateEstimator(passenger_data)
    are.MakePlot()


    

if __name__ == '__main__':
    main(*sys.argv)
