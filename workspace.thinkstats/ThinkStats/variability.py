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

NUM_SIGMAS = 1

class Height(thinkbayes.Suite, thinkbayes.Joint):

    def __init__(self, mus, sigmas, name=''):
        """Makes a prior distribution for mu and sigma based on a sample.

        mus: sequence of possible mus
        sigmas: sequence of possible sigmas
        name: string name for the Suite
        """
        self.mus = mus
        self.sigmas = sigmas

        pairs = [(mu, sigma) 
                 for mu in mus
                 for sigma in sigmas]

        thinkbayes.Suite.__init__(self, pairs, name=name)

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
        loglike = EvalGaussianLogPdf(mu, sigma, x)
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
        S = math.sqrt(S2)

        self.LogUpdateSetABC(n, xbar, S)

    def LogUpdateSetMedianIPR(self, data):
        """Computes the log likelihood of the data under the hypothesis.

        Estimates log likelihood using Approximate Bayesian Computation (ABC).

        Args:
            data: sequence of values
        """
        xs = data
        n = len(xs)

        # compute summary stats
        median, sighat = MedianSighat(xs, num_sigmas=NUM_SIGMAS)
        print 'median, sighat', median, sighat

        self.LogUpdateSetABC(n, median, sighat)

    def LogUpdateSetABC(self, n, xbar, sighat):
        for hypo in sorted(self.Values()):
            mu, sigma = hypo

            # compute log likelihood of xbar, given hypo
            sample_mu = mu
            sample_sigma = sigma / math.sqrt(n)
            loglike = EvalGaussianLogPdf(sample_mu, sample_sigma, xbar)

            #compute log likelihood of sighat, given hypo
            sample_mu = sigma
            sample_sigma = sigma / math.sqrt(2 * (n-1))
            loglike += EvalGaussianLogPdf(sample_mu, sample_sigma, sighat)

            self.Incr(hypo, loglike)


def EvalGaussianLogPdf(mu, sigma, x):
    """Computes the log PDF of x given mu and sigma.

    mu, sigma: paramemters of Gaussian
    x: float values

    returns: float log-likelihood
    """
    z = (x-mu) / sigma
    logpdf = -math.log(sigma) - z**2 / 2
    return logpdf


def FindPriorRanges(xs, num_points, num_stderrs=3.0, median_flag=False):
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
    if median_flag:
        xbar, S = MedianSighat(xs, num_sigmas=NUM_SIGMAS)
    else:
        xbar, S2 = thinkstats.MeanVar(xs)
        S = math.sqrt(S2)

    print 'classical estimators', xbar, S

    # compute ranges for xbar and S
    stderr_xbar = S / math.sqrt(n)
    mus = MakeRange(xbar, stderr_xbar)

    stderr_S = S / math.sqrt(2 * (n-1))
    sigmas = MakeRange(S, stderr_S)

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


def ComputeCoefVariation(suite):
    """Computes the distribution of CV.

    suite: Pmf that maps (x, y) to z

    Returns: Pmf object for CV.
    """
    pmf = thinkbayes.Pmf()
    for (m, s), p in suite.Items():
        pmf.Incr(s/m, p)
    return pmf


def PlotCdfs(d, labels):

    myplot.Clf()
    for key, xs in d.iteritems():
        mu = thinkstats.Mean(xs)
        xs = thinkstats.Jitter(xs, 1.3)
        xs = [x-mu for x in xs]
        cdf = thinkbayes.MakeCdfFromList(xs)
        myplot.Cdf(cdf, label=labels[key])
    myplot.Show()
                  

def PlotPosterior(suite, pcolor=False, contour=True):
    """Makes a contour plot.
    
    suite: Suite that maps (mu, sigma) to probability
    """
    myplot.Clf()
    myplot.Contour(suite.GetDict(), pcolor=pcolor, contour=contour)

    myplot.Save(root='variability_posterior_%s' % suite.name,
                title='Posterior joint distribution',
                xlabel='Mean height (cm)',
                ylabel='Stddev (cm)')


def PlotCoefVariation(suites):
    """Plot the posterior distributions for CV.

    suites: map from label to Pmf of CVs.
    """
    myplot.Clf()
    myplot.PrePlot(num=2)

    pmfs = {}
    for label, suite in suites.iteritems():
        pmf = ComputeCoefVariation(suite)
        print 'CV posterior mean', pmf.Mean()
        cdf = thinkbayes.MakeCdfFromPmf(pmf, label)
        myplot.Cdf(cdf)
    
        pmfs[label] = pmf

    myplot.Save(root='variability_cv',
                xlabel='Coefficient of variation',
                ylabel='Probability')

    print 'female bigger', thinkbayes.PmfProbGreater(pmfs['female'],
                                                     pmfs['male'])
    print 'male bigger', thinkbayes.PmfProbGreater(pmfs['male'],
                                                   pmfs['female'])


def PlotOutliers(samples):
    """Make CDFs showing the distribution of outliers."""
    cdfs = []
    for label, sample in samples.iteritems():
        outliers = [x for x in sample if x < 150]

        cdf = thinkbayes.MakeCdfFromList(outliers, label)
        cdfs.append(cdf)

    myplot.Clf()
    myplot.Cdfs(cdfs)
    myplot.Save(root='variability_cdfs',
                title='CDF of height',
                xlabel='Reported height (cm)',
                ylabel='CDF')

def PlotMarginals(suite):
    myplot.Clf()

    pyplot.subplot(1, 2, 1)
    pmf_m = suite.Marginal(0)
    cdf_m = thinkbayes.MakeCdfFromPmf(pmf_m)
    myplot.Cdf(cdf_m)

    pyplot.subplot(1, 2, 2)
    pmf_s = suite.Marginal(1)
    cdf_s = thinkbayes.MakeCdfFromPmf(pmf_s)
    myplot.Cdf(cdf_s)

    myplot.Show()


def DumpHeights(data_dir='.', n=10000):
    """Read the BRFSS dataset, extract the heights and pickle them."""
    resp = brfss.Respondents()
    resp.ReadRecords(data_dir, n)

    d = {1:[], 2:[]}
    [d[r.sex].append(r.htm3) for r in resp.records if r.htm3 != 'NA']

    fp = open('variability_data.pkl', 'wb')
    cPickle.dump(d, fp)
    fp.close()


def LoadHeights():
    """Read the pickled height data.

    returns: map from sex code to list of heights.
    """
    fp = open('variability_data.pkl', 'r')
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


def UpdateSuite5(suite, xs):
    """Computes the posterior distibution of mu and sigma.

    Computes log likelihoods efficiently.

    suite: Suite that maps from (mu, sigma) to prob
    t: sequence
    """
    suite.Log()
    suite.LogUpdateSetMedianIPR(xs)
    suite.Exp()
    suite.Normalize()


def MedianIPR(xs, p):
    """Computes the median and interpercentile range.

    xs: sequence of values

    returns: tuple of float (median, IPR)
    """
    cdf = thinkbayes.MakeCdfFromList(xs)
    median = cdf.Percentile(50)

    alpha = (1-p) / 2
    ipr = cdf.Value(1-alpha) - cdf.Value(alpha)
    return median, ipr


def MedianSighat(xs, num_sigmas):
    """Computes the median and an estimate of sigma.

    Based on an interpercentile range (IPR).

    factor: number of standard deviations spanned by the IPR
    """
    half_p = thinkbayes.StandardGaussianCdf(num_sigmas) - 0.5
    median, ipr = MedianIPR(xs, half_p * 2)
    sighat = ipr / 2 / num_sigmas

    return median, sighat

def Summarize(xs):
    # print outliers
    xs.sort()
    print 'smallest', xs[:10]
    print 'largest', xs[-10:]

    cdf = thinkbayes.MakeCdfFromList(xs)
    print cdf.Percentile(25), cdf.Percentile(50), cdf.Percentile(75)


def RunEstimate(update_func, num_points=31, median_flag=False):
    """Runs the whole analysis.

    update_func: which of the update functions to use
    num_points: number of points in the Suite (in each dimension)
    """
    # DumpHeights(n=10000000)
    d = LoadHeights()
    labels = {1:'male', 2:'female'}

    # PlotCdfs(d, labels)

    suites = {}
    for key, xs in d.iteritems():
        name = labels[key]
        print name, len(xs)
        Summarize(xs)

        xs = thinkstats.Jitter(xs, 1.3)

        mus, sigmas = FindPriorRanges(xs, num_points, median_flag=median_flag)
        suite = Height(mus, sigmas, name)
        suites[name] = suite
        update_func(suite, xs)
        print 'MLE', suite.MaximumLikelihood()

        PlotPosterior(suite)

        pmf_m = suite.Marginal(0)
        pmf_s = suite.Marginal(1)
        print 'marginal mu', pmf_m.Mean(), pmf_m.Var()
        print 'marginal sigma', pmf_s.Mean(), pmf_s.Var()

        # PlotMarginals(suite)

    PlotCoefVariation(suites)


def main():
    random.seed(17)

    func = UpdateSuite5
    median_flag = (func == UpdateSuite5)
    RunEstimate(func, median_flag=median_flag)


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

