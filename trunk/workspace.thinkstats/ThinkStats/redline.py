"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
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

formats = ['pdf']

"""
Notation guide:

z: time between trains
x: time since the last train
y: time until the next train

z': distribution of z as seen by a random arrival

"""


# interarrival_times in seconds, collected by Kai Austin and Brendan Ritter
# using http://developer.mbta.com/Data/Red.txt

observed_interarrival_times = [
    365, 47, 146, 545, 132, 270, 443, 190, 372, 167, 
    375, 128, 455, 262, 233, 386, 561, 386, 241, 562, 489, 
    455, 410, 489, 342, 347, 170, 260, 148, 375, 312, 265, 
    201, 245, 567, 633, 627, 643, 159, 827, 773, 159, 1100, 160, 
    148, 187, 290, 353, 133, 180, 355, 151, 558, 220, 232, 
    353, 199, 160, 172
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


def WeibullSample(n, lam, k):
    return [random.weibullvariate(lam, k) for i in range(n)]


def ActualInterarrivals(n, lam=12, k=1.5):
    sample = WeibullSample(n, lam, k)
    pmf = thinkbayes.MakePmfFromList(sample)
    pmf.name = 'actual'
    
def BiasedInterarrivals(n, lam=12, k=1.5):
    actual_pmf = ActualInterarrivals(1000)
    observed_pmf = BiasPmf(actual_pmf, 'observed')


def MakeUniformPmf(low, high):
    """Make a uniform Pmf.

    low: lowest value (inclusive)
    high: highest value (inclusize)
    """
    pmf = thinkbayes.Pmf()
    for x in MakeRange(low=low, high=high):
        pmf.Set(x, 1)
    pmf.Normalize()
    return pmf    
    

def MakeRange(low=0, high=1300, skip=10):
    return range(low, high+skip, skip)


class Forward(object):
    """Encapsulates the forward inference process.

    Given the actual distribution of interarrival times (z),
    computes the distribution of interarrivals as seen by
    a random passenger (z'), which yields the distribution
    of wait times (y) and the distribution of elapsed times (x).
    """

    def __init__(self, interarrival_times):
        self.interarrival_times = interarrival_times

        # distribution of interarrival time (z)
        self.pdf_z = thinkbayes.EstimatedPdf(interarrival_times)
        self.xs = MakeRange(low=10)
        self.pmf_z = self.pdf_z.MakePmf(self.xs, name='z')

        # distribution of interarrival time as seen by a random
        # observer (z')
        self.pmf_zp = BiasPmf(self.pmf_z, name="z'")

        # distribution of elapsed time (x)
        self.pmf_x = self.PmfOfWaitTime(self.pmf_zp)

        # the distribution of wait time (y) is the same as the
        # distribution of elapsed time (x)
        self.pmf_y = self.pmf_x

    def PmfOfWaitTime(self, pmf_zp):
        """Distribution of wait time.

        pmf_zp: dist of interarrival time as seen by a random observer

        Returns: dist of wait time (also dist of elapsed time)
        """
        meta_pmf = thinkbayes.Pmf()
        for interarrival, prob in pmf_zp.Items():
            uniform = MakeUniformPmf(10, interarrival)
            meta_pmf.Set(uniform, prob)

        pmf_y = thinkbayes.MakeMixture(meta_pmf, name='y')
        return pmf_y

    def GenerateSampleWaitTimes(self, n):
        """Generates a random sample of wait times.

        n: sample size

        Returns: sequence of values
        """
        cdf_y = thinkbayes.MakeCdfFromPmf(self.pmf_y)
        sample = cdf_y.Sample(n)
        return sample

    def PlotPmfs(self):
        """Plots the computed Pmfs.
        """
        print 'Mean z', self.pmf_z.Mean()
        print 'Mean zp', self.pmf_zp.Mean()
        print 'Mean y', self.pmf_y.Mean()

        myplot.Pmf(self.pmf_z)
        myplot.Pmf(self.pmf_zp)
        myplot.Pmf(self.pmf_y)


class Interarrivals(thinkbayes.Suite):
    """Represents the distribution of interarrival times,
    as updated by an observed waiting time."""

    def Likelihood(self, hypo, data):
        """The likelihood of the data under the hypothesis.

        If the actual interarrival time is z, what is the likelihood
        of waiting y seconds?

        hypo: actual time between trains
        data: observed wait time
        """
        z = hypo
        y = data
        if y > z:
            return 0
        return 1.0 / z


class InterarrivalDirichlet(thinkbayes.Dirichlet):
    """Represents the distribution of prevalences for each
    interarrival time."""

    def __init__(self, xs):
        """Constructor.

        xs: sequence of possible interarrival times
        """
        n = len(xs)
        thinkbayes.Dirichlet.__init__(self, n)
        self.xs = xs

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
        y = data

        print y
        prior = self.PredictivePmf(self.xs)
        interarrivals = Interarrivals(prior)
        interarrivals.Update(y)
        probs = interarrivals.Probs(self.xs)

        print len(probs)
        print len(self.params)
        self.params += numpy.array(probs)

    def InterarrivalProbs(self, wait_time):
        """Not currently working
        """
        xs = range(len(self.params))
        ps = self.Random()
        prior = thinkbayes.MakePmfFromItems(zip(xs, ps))
        interarrivals = Interarrivals(prior)
        interarrivals.Update(wait_time)
        return interarrivals.Probs(needs_xs)


class Reverse(object):
    """Represents the reverse inference process.

    """

    def __init__(self, pcounts, wait_times):
        self.pcounts = pcounts
        self.wait_times = wait_times
        self.pdf_y = thinkbayes.EstimatedPdf(wait_times)
        self.xs = MakeRange(low=10)
        self.pmf_y = self.pdf_y.MakePmf(self.xs, name='y')
        self.pmf_y = thinkbayes.MakePmfFromList(wait_times, name='y')

        self.pmf_z = self.PmfOfInterarrivalTime(wait_times)
        self.pmf_zp = UnbiasPmf(self.pmf_z, name="z'")

    def PmfOfInterarrivalTime(self, wait_times):
        n = len(self.xs)
        dirichlet = InterarrivalDirichlet(self.xs)
        dirichlet.Preload(self.pcounts)
        self.prior_z = dirichlet.PredictivePmf(self.xs, name='prior')
        
        dirichlet.params /= 100000.0

        for y in wait_times:
            dirichlet.Update(y)

        pmf_z = dirichlet.PredictivePmf(self.xs, name='z')
        return pmf_z

    def PlotPmfs(self):
        print 'Mean y', self.pmf_y.Mean()
        print 'Mean z', self.pmf_z.Mean()
        print "Mean z'", self.pmf_zp.Mean()

        myplot.Pmf(self.pmf_y)
        myplot.Pmf(self.pmf_z)
        myplot.Pmf(self.pmf_zp)

    def PlotCdfs(self):
        myplot.Cdf(self.pmf_y.MakeCdf())
        myplot.Cdf(self.prior_z.MakeCdf())
        myplot.Cdf(self.pmf_z.MakeCdf())
        #myplot.Cdf(self.pmf_zp.MakeCdf())


def Floor(x, factor=10):
    """Rounds down to the nearest multiple of factor.

    When factor=10, all numbers from 10 to 19 get floored to 10.
    """
    return int(x/factor) * factor


def RunForward(interarrival_times):
    """Runs the forward inference process.
    """
    forward = Forward(interarrival_times)

    cdf_z = thinkbayes.MakeCdfFromList(interarrival_times, name='z obs')
    cdf_zp = forward.pmf_zp.MakeCdf(name="z'")
    cdf_y = forward.pmf_y.MakeCdf(name='y')
    myplot.Cdfs([cdf_z, cdf_zp, cdf_y])

    n = len(interarrival_times)
    sample_y = forward.GenerateSampleWaitTimes(n)

    cdf_y_sample = thinkbayes.MakeCdfFromList(sample_y, 'y sample')
    myplot.Cdf(cdf_y_sample)
    myplot.Show()

    return forward


def RunElapsed(interarrival_times, lam, num_passengers):
    """Estimates the elapsed time based on number of passengers
    and known lam.
    """
    forward = Forward(interarrival_times)
    cdf_z = thinkbayes.MakeCdfFromList(interarrival_times, name='z obs')

    prior = forward.pmf_x
    elapsed = Elapsed(lam, prior, 'elapsed prior')
    cdf_prior = elapsed.MakeCdf()

    elapsed.Update(num_passengers)
    elapsed.name = 'elapsed posterior'
    cdf_posterior = elapsed.MakeCdf()

    myplot.Cdfs([cdf_z, cdf_prior, cdf_posterior])
    myplot.Show()

    return forward, elapsed


class Elapsed(thinkbayes.Suite):
    """Represents the distribution of elapsed time (x)."""

    def __init__(self, lam, hypos, name=''):
        thinkbayes.Suite.__init__(self, hypos, name)
        self.lam = lam

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        Evaluates the Poisson PMF for lambda and k.

        hypo: elapsed time since the last train
        data: number of passengers on the platform
        """
        elapsed_time = hypo
        k = data
        like = thinkbayes.EvalPoissonPmf(self.lam * elapsed_time, k)
        return like

    def WaitTime(self, pmf_zp):
        """Computes the distribution of wait times.

        Enumerate all pairs of zp from pmf_zp and x from self,
        and accumulate the distribution of y = z - x.

        pmf_zp: distribution of interarrivals seen by random observer
        """
        return pmf_zp - self


def main(script):

    # forward = RunForward(observed_interarrival_times)
    elapsed = RunElapsed(observed_interarrival_times, 
                         lam=0.0333,
                         num_passengers=20)

    return

    # interarrival_times = [200, 300]*10


    # use the actual interarrival times to make pcounts
    vals = [Floor(t) for t in interarrival_times]
    hist = thinkbayes.MakeHistFromList(vals)
    for val, freq in hist.Items():
        print val, freq

    xs = MakeRange(low=10)
    pcounts = hist.Freqs(xs)
    print len(pcounts)
    print pcounts

    # do reverse inference
    reverse = Reverse(pcounts, sample_y)
    reverse.PlotCdfs()

    myplot.Show()

    return
    
    for y in [100, 200, 300, 1100]:
        inter = Interarrivals(y)
        inter.Update(y)
        myplot.Pmf(inter)
        #myplot.Cdf(inter.MakeCdf())
    myplot.Show()

    return

    

    # arrival rate of passengers in passengers / second
    lam_pass = 0.05
    
    
    
    return



    n = 300
    lam = 12
    k = 0.9
    sample = WeibullSample(n, lam, k)
    cdf = thinkbayes.MakeCdfFromList(sample)
    myplot.Cdf(cdf)
    myplot.Show()

    print 'mean', cdf.Mean()
    print 'median', cdf.Percentile(50)
    print '25th', cdf.Percentile(25)
    print '75th', cdf.Percentile(75)

if __name__ == '__main__':
    main(*sys.argv)
