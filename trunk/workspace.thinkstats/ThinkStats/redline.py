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

z': distribution of z as seen by a random arrival

"""


# interarrival_times in seconds, collected by Kai Austin and Brendan Ritter
# using http://developer.mbta.com/Data/Red.txt

# actual time between trains in seconds

observed_interarrival_times = [
    365, 47, 146, 545, 132, 270, 443, 190, 372, 167, 
    375, 128, 455, 262, 233, 386, 561, 386, 241, 562, 489, 
    455, 410, 489, 342, 347, 170, 260, 148, 375, 312, 265, 
    201, 245, 567, 633, 627, 643, 159, 827, 773, 159, 1100, 160, 
    148, 187, 290, 353, 133, 180, 355, 151, 558, 220, 232, 
    353, 199, 160, 172
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


def MakeRange(low=10, high=1300, skip=10):
    """Makes a range representing possible interarrival times in seconds.

    low: where to start
    high: where to end
    skip: how many to skip
    """
    return range(low, high+skip, skip)


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
    

class WaitTimeCalculator(object):
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
        self.xs = MakeRange()
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

    def GenerateSamplePassengers(self, lam, n, high=100):
        """Generates a sample wait time and number of arrivals.

        lam: arrival rate in passengers per second
        n: number of samples

        Returns: list of (wait time, number of passenger) pairs
        """
        wait_times = self.GenerateSampleWaitTimes(n)

        res = []
        for y in wait_times:
            pmf_k = thinkbayes.MakePoissonPmf(lam * y, high)
            k = pmf_k.Random()
            res.append((y, k))

        return res

    def PlotPmfs(self):
        """Plots the computed Pmfs.
        """
        print 'Mean z', self.pmf_z.Mean()
        print 'Mean zp', self.pmf_zp.Mean()
        print 'Mean y', self.pmf_y.Mean()

        thinkplot.Pmf(self.pmf_z)
        thinkplot.Pmf(self.pmf_zp)
        thinkplot.Pmf(self.pmf_y)

    def MakePlot(self, root='redline2'):
        """Plots the computed CDFs.
        """
        cdf_z = self.pmf_z.MakeCdf()
        cdf_zp = self.pmf_zp.MakeCdf()
        cdf_y = self.pmf_y.MakeCdf()

        cdfs = ScaleCdfs([cdf_y, cdf_z, cdf_zp], 1.0/60)

        thinkplot.Clf()
        thinkplot.PrePlot(3)
        thinkplot.Cdfs(cdfs)
        thinkplot.Save(root=root,
                       xlabel='Time (min)',
                       ylabel='CDF',
                       formats=FORMATS)


def ScaleCdfs(cdfs, factor):
    [cdf.Scale(factor) for cdf in cdfs]
    return cdfs


class ElapsedTimeEstimator(object):

    def __init__(self, wtc, lam, num_passengers):
        """Constructor

        interarrival_times: sample of actual interarrival times
        lam: arrival rate in passengers per second
        num_passengers: # passengers seen on the platform
        """
        self.wtc = wtc

        # prior for elapsed time
        self.prior = Elapsed(wtc.pmf_x, 'prior x')

        # posterior of elapsed time (based on number of passengers)
        self.posterior = self.prior.Copy(name='posterior x')
        self.posterior.Update((lam, num_passengers))

        # predictive distribution of wait time
        self.pmf_y = self.posterior.PredictWaitTime(wtc.pmf_zp)

    def MakePlot(self, root='redline3'):

        # observed interarrivals
        cdf_z = self.wtc.pmf_z.MakeCdf()
        cdf_prior = self.prior.MakeCdf()
        cdf_posterior = self.posterior.MakeCdf()
        cdf_y = self.pmf_y.MakeCdf()

        cdfs = ScaleCdfs([cdf_prior, cdf_posterior, cdf_y], 1.0/60)

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
        elapsed_time, k = data
        like = thinkbayes.EvalPoissonPmf(lam * elapsed_time, k)
        return like


class ArrivalRateEstimator(object):

    def __init__(self, passenger_data):
        """Constructor

        passenger_data: sequence of (y, k) pairs
        """
        self.passenger_data = passenger_data

        # range for lambda in passengers per second
        low, high = 0, 5
        n = 101
        hypos = numpy.linspace(low, high, n) / 60

        self.prior = ArrivalRate(hypos, name='prior')

        self.posterior = self.prior.Copy(name='posterior')

        for y, k in passenger_data:
            self.posterior.Update((y, k))

        print 'Mean posterior lambda', self.posterior.Mean()

    def MakePlot(self, root='redline1'):
        thinkplot.Clf()

        # convert units to passengers per minute
        cdf = self.posterior.MakeCdf()
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
        elapsed_time = hypo
        lam, k = data
        like = thinkbayes.EvalPoissonPmf(lam * elapsed_time, k)
        return like

    def PredictWaitTime(self, pmf_zp):
        """Computes the distribution of wait times.

        Enumerate all pairs of zp from pmf_zp and x from self,
        and accumulate the distribution of y = z - x.

        pmf_zp: distribution of interarrivals seen by random observer
        """
        pmf_y = pmf_zp - self
        pmf_y.name = 'y'
        RemoveNegatives(pmf_y)
        return pmf_y


def RemoveNegatives(pmf):
    """Removes negative values from a PMF.

    pmf: Pmf
    """
    for val in pmf.Values():
        if val < 0:
            pmf.Remove(val)


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


class InterarrivalTimeEstimator(object):
    """Represents the reverse inference process.

    """

    def __init__(self, pcounts, wait_times):
        self.pcounts = pcounts
        self.wait_times = wait_times
        self.xs = MakeRange()

        # smooth the observed wait times?
        #self.pdf_y = thinkbayes.EstimatedPdf(wait_times)
        #self.xs = MakeRange()
        #self.pmf_y = self.pdf_y.MakePmf(self.xs, name='y')

        self.pmf_y = thinkbayes.MakePmfFromList(wait_times, name='y')
        thinkplot.Cdf(self.pmf_y.MakeCdf())

        self.pmf_z = self.PmfOfInterarrivalTime(wait_times)
        self.pmf_zp = UnbiasPmf(self.pmf_z, name="z'")

    def PmfOfInterarrivalTime(self, wait_times):
        n = len(self.xs)
        dirichlet = InterarrivalDirichlet(self.xs)
        dirichlet.params /= 1.0

        # self.prior_z = dirichlet.PredictivePmf(self.xs, name='prior')
        # thinkplot.Cdf(self.prior_z.MakeCdf())

        dirichlet.Preload(self.pcounts)
        dirichlet.params /= 10000.0

        self.prior_z = dirichlet.PredictivePmf(self.xs, name='preloaded')
        thinkplot.Cdf(self.prior_z.MakeCdf())
        
        for y in wait_times:
            dirichlet.Update(y)

        pmf_z = dirichlet.PredictivePmf(self.xs, name='z')

        self.posterior_z = dirichlet.PredictivePmf(self.xs, name='posterior')
        thinkplot.Cdf(self.posterior_z.MakeCdf())

        return pmf_z

    def PlotPmfs(self):
        print 'Mean y', self.pmf_y.Mean()
        print 'Mean z', self.pmf_z.Mean()
        print "Mean z'", self.pmf_zp.Mean()

        thinkplot.Pmf(self.pmf_y)
        thinkplot.Pmf(self.pmf_z)
        thinkplot.Pmf(self.pmf_zp)

    def MakePlot(self):
        thinkplot.Cdf(self.pmf_y.MakeCdf())
        thinkplot.Cdf(self.prior_z.MakeCdf())
        thinkplot.Cdf(self.pmf_z.MakeCdf())
        #thinkplot.Cdf(self.pmf_zp.MakeCdf())


def Floor(x, factor=10):
    """Rounds down to the nearest multiple of factor.

    When factor=10, all numbers from 10 to 19 get floored to 10.
    """
    return int(x/factor) * factor


def MakePcounts(interarrival_times):

    # use the actual interarrival times to make pcounts
    vals = [Floor(t) for t in interarrival_times]
    hist = thinkbayes.MakeHistFromList(vals)

    xs = MakeRange()
    pcounts = hist.Freqs(xs)

    return pcounts


def main(script):
    cdf_z = thinkbayes.MakeCdfFromList(observed_interarrival_times, name='z')
    thinkplot.Cdf(cdf_z)
    
    pcounts = MakePcounts(observed_interarrival_times)
    print pcounts

    wait_times = [t for t, k in observed_passenger_data ]
    print wait_times

    ite = InterarrivalTimeEstimator(pcounts, wait_times)
    thinkplot.Show()

    #ite.MakePlot()

    return

    are = ArrivalRateEstimator(observed_passenger_data)
    are.MakePlot()


    wtc = WaitTimeCalculator(observed_interarrival_times)
    wtc.MakePlot()


    elapsed = ElapsedTimeEstimator(wtc,
                                   lam=0.0333,
                                   num_passengers=10)
    elapsed.MakePlot()



    return
    
    for y in [100, 200, 300, 1100]:
        inter = Interarrivals(y)
        inter.Update(y)
        thinkplot.Pmf(inter)
        #thinkplot.Cdf(inter.MakeCdf())
    thinkplot.Show()

    return

    

    # arrival rate of passengers in passengers / second
    lam_pass = 0.05
    
    
    
    return



    n = 300
    lam = 12
    k = 0.9
    sample = WeibullSample(n, lam, k)
    cdf = thinkbayes.MakeCdfFromList(sample)
    thinkplot.Cdf(cdf)
    thinkplot.Show()

    print 'mean', cdf.Mean()
    print 'median', cdf.Percentile(50)
    print '25th', cdf.Percentile(25)
    print '75th', cdf.Percentile(75)

if __name__ == '__main__':
    main(*sys.argv)
