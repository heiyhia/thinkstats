"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import numpy
import cPickle
import random

import brfss
import correlation
import myplot

import thinkbayes
import thinkstats

import matplotlib.pyplot as pyplot


class Height(thinkbayes.Suite):

    def __init__(self, mus, sigmas, name=''):
        """Makes a prior distribution for mu and sigma based on a sample.

        mus: sequence of possible mus
        sigmas: sequence of possible sigmas
        name: string name for the Suite
        """
        thinkbayes.Suite.__init__(self, name=name)

        self.mus = mus
        self.sigmas = sigmas

        # populate the Suite
        for mu in self.mus:
            for sigma in self.sigmas:
                self.Set((mu, sigma), 1)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        Args:
            hypo: tuple of hypothetical mu and sigma
            data: float sample

        Returns:
            likelihood of the sample given mu and sigma
        """
        x = data
        mu, sigma = hypo
        like = thinkbayes.EvalGaussianPdf(mu, sigma, x)
        return like

    def LogLikelihood(self, hypo, data):
        """Computes the log likelihood of the data under the hypothesis.

        Args:
            data: a list of values
            hypo: tuple of hypothetical mu and sigma

        Returns:
            log likelihood of the sample given mu and sigma (unnormalized)
        """
        x = data
        mu, sigma = hypo
        loglike = GaussianLogLikelihood(mu, sigma, x)
        return loglike

    def LogUpdateSetFast(self, data):
        """Computes the log likelihood of the data under the hypothesis.

        Args:
            data: sequence of values
        """
        xs = tuple(data)
        n = len(xs)

        for hypo in self.Values():
            mu, sigma = hypo
            total = Summation(xs, mu)
            loglike = -n * math.log(sigma) - total / 2 / sigma**2
            self.Incr(hypo, loglike)

    def LogUpdateSetMeanVar(self, data):
        """Computes the log likelihood of the data under the hypothesis.

        Estimates log likelihood using Approximate Bayesian Computation (ABC).

        Args:
            data: sequence of values
        """
        xs = data
        n = len(xs)

        # compute summary stats
        xbar, S2 = thinkstats.MeanVar(xs)
        sighat = math.sqrt(S2)

        self.LogUpdateSetABC(n, xbar, sighat)

    def LogUpdateSetMedianIQR(self, data):
        """Computes the log likelihood of the data under the hypothesis.

        Estimates log likelihood using Approximate Bayesian Computation (ABC).

        Args:
            data: sequence of values
        """
        xs = data
        n = len(xs)

        # compute summary stats
        median, iqr = MedianIQR(xs)
        sighat = iqr / 1.349
        print 'median, sighat', median, sighat

        self.LogUpdateSetABC(n, median, sighat)

    def LogUpdateSetABC(self, n, xbar, sighat):
        for hypo in sorted(self.Values()):
            mu, sigma = hypo

            # compute log likelihood of xbar, given hypo
            sample_mu = mu
            sample_sigma = sigma / math.sqrt(n)
            loglike = GaussianLogLikelihood(sample_mu, sample_sigma, xbar)

            #compute log likelihood of sighat, given hypo
            sample_mu = sigma
            sample_sigma = sigma / math.sqrt(2 * (n-1))
            loglike += GaussianLogLikelihood(sample_mu, sample_sigma, sighat)

            self.Incr(hypo, loglike)


def GaussianLogLikelihood(mu, sigma, x):
    """Computes the log-likelihood of x given mu and sigma.

    mu, sigma: paramemters of Gaussian
    x: float values

    returns: float log-likelihood
    """
    z = (x-mu) / sigma
    loglike = -math.log(sigma) - z**2 / 2
    return loglike


def FindPriorRanges(xs, num_points, num_stderrs=3.0):
    """Find ranges for mu and sigma with non-negligible likelihood.

    xs: sample
    num_points: number of values in each dimension
    num_stderrs: number of standard errors to include on either side
    
    Returns: sequence of mus, sequence of sigmas    
    """
    def MakeRange(estimate, stderr):
        """Makes a linear range around the estimate.

        estimate: central value
        stderr: standard error of the estimate

        returns: numpy array of float
        """
        spread = stderr * num_stderrs
        array = numpy.linspace(estimate-spread, estimate+spread, num_points)
        return array

    # estimate mean and stddev of xs
    n = len(xs)
    xbar, S2 = thinkstats.MeanVar(xs)
    sighat = math.sqrt(S2)

    print 'classical estimators', xbar, sighat

    # compute ranges for xbar and sighat
    stderr_xbar = sighat / math.sqrt(n)
    mus = MakeRange(xbar, stderr_xbar)

    stderr_sighat = sighat / math.sqrt(2 * (n-1))
    sigmas = MakeRange(sighat, stderr_sighat)

    return mus, sigmas


def Summation(xs, mu, cache={}):
    """Computes the sum of (x-mu)**2 for x in t.

    Caches previous results.

    xs: tuple of values
    mu: hypothetical mean
    cache: cache of previous results
    """
    try:
        return cache[xs, mu]
    except KeyError:
        ds = [(x-mu)**2 for x in xs]
        total = sum(ds)
        cache[xs, mu] = total
        return total


def ComputeMarginals(suite):
    """Computes the marginal distributions for mu and sigma.

    suite: Pmf that maps (x, y) to z

    Returns: Pmf objects for mu and sigma
    """
    pmf_m = thinkbayes.Pmf()
    pmf_s = thinkbayes.Pmf()
    for (m, s), p in suite.Items():
        pmf_m.Incr(m, p)
        pmf_s.Incr(s, p)
    return pmf_m, pmf_s


def ComputeCoefVariation(suite):
    """Computes the distribution of CV.

    suite: Pmf that maps (x, y) to z

    Returns: Pmf object for CV.
    """
    pmf = thinkbayes.Pmf()
    for (m, s), p in suite.Items():
        pmf.Incr(s/m, p)
    return pmf


def PlotPosterior(suite, pcolor=False, contour=True):
    """Makes a contour plot.
    
    suite: Suite that maps (mu, sigma) to probability
    """
    X, Y = numpy.meshgrid(suite.mus, suite.sigmas)
    func = lambda x, y: suite.Prob((x, y))
    prob = numpy.vectorize(func)
    Z = prob(X, Y)

    myplot.Clf()
    if pcolor:
        pyplot.pcolor(X, Y, Z)
    if contour:
        pyplot.contour(X, Y, Z)

    myplot.Save(root='bayes_height_posterior_%s' % suite.name,
                title='Posterior joint distribution',
                xlabel='Mean height (cm)',
                ylabel='Stddev (cm)')


def PlotCoefVariation(suites):
    """Plot the posterior distributions for CV.

    suites: map from label to Pmf of CVs.
    """
    myplot.Clf()

    pmfs = {}
    for label, suite in suites.iteritems():
        pmf = ComputeCoefVariation(suite)
        cdf = thinkbayes.MakeCdfFromPmf(pmf, label)
        myplot.Cdf(cdf)
    
        pmfs[label] = pmf

    myplot.Save(root='bayes_height_cv',
                xlabel='Coefficient of variation',
                ylabel='Probability')

    print 'female bigger', thinkbayes.PmfProbGreater(pmfs['female'],
                                                     pmfs['male'])
    print 'male bigger', thinkbayes.PmfProbGreater(pmfs['male'],
                                                   pmfs['female'])


def PlotCdfs(samples):
    """Make CDFs showing the distribution of outliers."""
    cdfs = []
    for label, sample in samples.iteritems():
        outliers = [x for x in sample if x < 150]

        cdf = thinkbayes.MakeCdfFromList(outliers, label)
        cdfs.append(cdf)

    myplot.Clf()
    myplot.Cdfs(cdfs)
    myplot.Save(root='bayes_height_cdfs',
                title='CDF of height',
                xlabel='Reported height (cm)',
                ylabel='CDF')


def DumpHeights(data_dir='.', n=10000):
    """Read the BRFSS dataset, extract the heights and pickle them."""
    resp = brfss.Respondents()
    resp.ReadRecords(data_dir, n)

    d = {1:[], 2:[]}
    [d[r.sex].append(r.htm3) for r in resp.records if r.htm3 != 'NA']

    fp = open('bayes_height_data.pkl', 'wb')
    cPickle.dump(d, fp)
    fp.close()


def LoadHeights():
    """Read the pickled height data.

    returns: map from sex code to list of heights.
    """
    fp = open('bayes_height_data.pkl', 'r')
    d = cPickle.load(fp)
    fp.close()
    return d


def UpdateSuite1(suite, xs):
    """Computes the posterior distibution of mu and sigma.

    Computes untransformed likelihoods.

    suite: Suite that maps from (mu, sigma) to prob
    xs: sequence
    """
    suite.UpdateSet(xs)


def UpdateSuite2(suite, xs):
    """Computes the posterior distibution of mu and sigma.

    Computes log likelihoods.

    suite: Suite that maps from (mu, sigma) to prob
    xs: sequence
    """
    suite.Log()
    suite.LogUpdateSet(xs)
    suite.Exp()
    suite.Normalize()


def UpdateSuite3(suite, xs):
    """Computes the posterior distibution of mu and sigma.

    Computes log likelihoods efficiently.

    suite: Suite that maps from (mu, sigma) to prob
    t: sequence
    """
    suite.Log()
    suite.LogUpdateSetFast(xs)
    suite.Exp()
    suite.Normalize()


def UpdateSuite4(suite, xs):
    """Computes the posterior distibution of mu and sigma.

    Computes log likelihoods efficiently.

    suite: Suite that maps from (mu, sigma) to prob
    t: sequence
    """
    suite.Log()
    suite.LogUpdateSetMeanVar(xs)
    suite.Exp()
    suite.Normalize()


def MedianIQR(xs):
    """Computes the median and interquartile range.

    xs: sequence of values

    returns: tuple of float (median, IQR)
    """
    cdf = thinkbayes.MakeCdfFromList(xs)
    median = cdf.Percentile(50)
    iqr = cdf.Percentile(75) - cdf.Percentile(25)
    return median, iqr


def UpdateSuite5(suite, xs):
    """Computes the posterior distibution of mu and sigma.

    Computes log likelihoods efficiently.

    suite: Suite that maps from (mu, sigma) to prob
    t: sequence
    """
    suite.Log()
    suite.LogUpdateSetMedianIQR(xs)
    suite.Exp()
    suite.Normalize()


def RunEstimate(update_func, num_points=31):
    """Runs the whole analysis.

    update_func: which of the update functions to use
    num_points: number of points in the Suite (in each dimension)
    """
    #DumpHeights(n=1000000)
    d = LoadHeights()

    labels = {1:'male', 2:'female'}

    samples = {}
    suites = {}

    for key, xs in d.iteritems():
        name = labels[key]
        print name, len(xs)

        mus, sigmas = FindPriorRanges(xs, num_points)
        suite = Height(mus, sigmas, name)
        suites[name] = suite
        update_func(suite, xs)

        PlotPosterior(suite)

        pmf_m, pmf_s = ComputeMarginals(suite)
        print 'marginal mu', pmf_m.Mean(), pmf_m.Var()
        print 'marginal sigma', pmf_s.Mean(), pmf_s.Var()

    PlotCoefVariation(suites)


def main():
    func = UpdateSuite5
    RunEstimate(func)


if __name__ == '__main__':
    main()


"""
UpdateSuite1 (100):
marginal mu 162.816901408 0.55779791443
marginal sigma 6.36966103214 0.277026082819

UpdateSuite2 (100):
marginal mu 162.816901408 0.55779791443
marginal sigma 6.36966103214 0.277026082819

UpdateSuite3 (100):
marginal mu 162.816901408 0.55779791443
marginal sigma 6.36966103214 0.277026082819

UpdateSuite4 (100):
marginal mu 162.816901408 0.547456009605
marginal sigma 6.30305516111 0.27544106054

UpdateSuite3 (1000):
marginal mu 163.722137405 0.0660294386397
marginal sigma 6.64453251495 0.0329935312671

UpdateSuite4 (1000):
marginal mu 163.722137405 0.0658920503302
marginal sigma 6.63692197049 0.0329689887609

UpdateSuite3 (all):
marginal mu 163.223475005 0.000203282582659
marginal sigma 7.26918836916 0.000101641131229

UpdateSuite4 (all):
marginal mu 163.223475004 0.000203281499857
marginal sigma 7.26916693422 0.000101640932082

UpdateSuite5 (all):
marginal mu 163.1805214 7.9399898468e-07
marginal sigma 7.29969524118 3.26257030869e-14

"""

