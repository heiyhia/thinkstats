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
            data: float sample
            hypo: tuple of hypothetical mu and sigma

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

        z = (x-mu) / sigma
        loglike = -math.log(sigma) - z**2 / 2
        return loglike

    def LogUpdateSet2(self, data):
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


def FindPriorRanges(xs, num_points, num_stderrs=3.0):
    """Find ranges for mu and sigma with non-negligible likelihood.

    xs: sample
    num_points: number of values in each dimension
    num_stderrs: number of standard errors to include on either side
    
    Returns: sequence of mus, sequence of sigmas    
    """
    # estimate mean and stddev of xs
    n = len(xs)
    xbar, S2 = thinkstats.MeanVar(xs)
    sighat = math.sqrt(S2)

    print 'classical estimators', xbar, sighat

    # compute standard error for mu and the range of mus
    stderr_xbar = sighat / math.sqrt(n)
    mspread = num_stderrs * stderr_xbar
    mus = numpy.linspace(xbar-mspread, xbar+mspread, num_points)

    # compute standard error for sigma and the range of ss
    stderr_sighat = sighat / math.sqrt(2 * (n-1))
    sspread = num_stderrs * stderr_sighat
    sigmas = numpy.linspace(sighat-sspread, sighat+sspread, num_points)

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


def NormalProbPlot(samples):
    """Makes a normal probability plot for each sample in samples."""
    myplot.Clf()

    markers = dict(male='b', female='g')

    for label, sample in samples.iteritems():
        NormalPlot(sample, label, markers[label], jitter=0.0)
    
    myplot.Save(show=True,
                #root='bayes_height_normal',
                title='Normal probability plot',
                xlabel='Standard normal',
                ylabel='Reported height (cm)')


def NormalPlot(ys, label, color='b', jitter=0.0, **line_options):
    """Makes a normal probability plot.
    
    Args:
        ys: sequence of values
        label: string label for the plotted line
        color: color string passed along to pyplot.plot
        jitter: float magnitude of jitter added to the ys 
        line_options: dictionary of options for pyplot.plot        
    """
    n = len(ys)
    xs = [random.gauss(0.0, 1.0) for i in range(n)]
    xs.sort()
    ys = [y + random.uniform(-jitter, +jitter) for y in ys]
    ys.sort()

    inter, slope = correlation.LeastSquares(xs, ys)
    fit = correlation.FitLine(xs, inter, slope)
    pyplot.plot(*fit, color=color, linewidth=0.5, alpha=0.5)

    pyplot.plot(sorted(xs), sorted(ys),
                color=color,
                marker='.',
                label=label,
                markersize=3,
                alpha=0.1,
                **line_options)
 

def PlotMarginals(suite):
    """Plot the marginal distributions for a 2-D joint distribution."""
    pmf_m, pmf_s = ComputeMarginals(suite)

    myplot.Clf()
    pyplot.figure(1, figsize=(7, 4))

    pyplot.subplot(1, 2, 1)
    cdf_m = thinkbayes.MakeCdfFromPmf(pmf_m, 'mu')
    myplot.Cdf(cdf_m)
    pyplot.xlabel('Mean height (cm)')
    pyplot.ylabel('CDF')

    pyplot.subplot(1, 2, 2)
    cdf_s = thinkbayes.MakeCdfFromPmf(pmf_s, 'sigma')
    myplot.Cdf(cdf_s)
    pyplot.xlabel('Std Dev height (cm)')
    pyplot.ylabel('CDF')

    myplot.Save(root='bayes_height_marginals_%s' % suite.name)


def PlotAges(resp):
    """Plot the distribution of ages."""
    ages = [r.age for r in resp.records]
    cdf = thinkbayes.MakeCdfFromList(ages)
    myplot.Clf()
    myplot.Cdf(cdf)
    myplot.Show()


def DumpHeights(data_dir='.', n=10000):
    """Read the BRFSS dataset, extract the heights and pickle them."""
    resp = brfss.Respondents()
    resp.ReadRecords(data_dir, n)

    #PlotAges(resp)

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


def Winsorize(xs, p=0.01):
    """Compresses outliers."""
    cdf = thinkbayes.MakeCdfFromList(xs)
    low, high = cdf.Value(p), cdf.Value(1-p)
    print low, high

    outliers = [x for x in xs if x < low or x > high]
    outliers.sort()
    print outliers

    wxs = [min(max(low, x), high) for x in xs]
    return wxs


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
    suite.LogUpdateSet2(xs)
    suite.Exp()
    suite.Normalize()


def RunEstimate(update_func, num_points=31):
    #DumpHeights(n=1000000)
    d = LoadHeights()

    labels = {1:'male', 2:'female'}

    samples = {}
    suites = {}

    for key, xs in d.iteritems():
        name = labels[key]
        print name, len(xs)

        #xs = Winsorize(xs, 0.0001)
        #samples[name] = xs

        mus, sigmas = FindPriorRanges(xs, num_points)
        suite = Height(mus, sigmas, name)
        suites[name] = suite
        update_func(suite, xs)
        print 'MLE', thinkbayes.MaximumLikelihood(suite)

        PlotPosterior(suite)
        PlotMarginals(suite)

    #PlotCdfs(samples)
    #NormalProbPlot(samples)
    PlotCoefVariation(suites)


def main():
    func = UpdateSuite3
    RunEstimate(func)
    return


if __name__ == '__main__':
    main()
