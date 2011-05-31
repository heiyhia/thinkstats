"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import bayes
import copy
import string
import matplotlib.pyplot as pyplot
import myplot
import Pmf
import Cdf
import math
import random
import sys
"""

"""

lowercase = string.ascii_lowercase


class Beta(object):
    """Represents a Beta distribution.

    See http://en.wikipedia.org/wiki/Beta_distribution

    alpha and beta are the parameters 

    """
    def __init__(self, yes, no, name=''):
        """Initializes a Beta distribution.

        yes and no are the number of successful/unsuccessful trials.
        """
        self.alpha = yes+1
        self.beta = no+1
        self.name = name

    def Update(self, yes, no):
        """Updates a Beta distribution."""
        self.alpha += yes
        self.beta += no
        
    def Mean(self):
        """Computes the mean of this distribution."""
        return float(self.alpha) / (self.alpha + self.beta)

    def Random(self):
        """Generates a random variate from this distribution."""

        if self.beta == 0:
            # handle the special case where there is one taxon, so
            # its prevalence is necessarily 1
            return 1

        return random.betavariate(self.alpha, self.beta)

    def Pdf(self, p):
        """Computes the PDF at p."""
        return math.pow(p, self.alpha-1) * math.pow(1-p, self.beta-1)
        
    def Pmf(self, steps=101):
        """Returns a curve that represents the PDF of this distribution."""
        ps = [i / (steps-1.0) for i in xrange(steps)]
        probs = [self.Pdf(p) for p in ps]
        pmf = Pmf.MakePmfFromDict(dict(zip(ps, probs)))
        return pmf

    def ConditionalCdf(self, fraction, steps=101):
        """Generates a CDF conditioned on p <= fraction."""
        ps = [fraction * i / (steps-1.0) for i in xrange(steps)]
        probs = [self.Pdf(p) for p in ps]
        cdf = Cdf.MakeCdfFromItems(zip(ps, probs))
        return cdf                                                    


class Census(object):
    """Represents a belief about the population of a sample.

    n is the hypothetical number of taxa

    multiplier is the number of taxa that have not been observed yet.

    taxa is a mapping from taxon to a distribution of p,
    the fractional population of the taxon.
    
    taxa['other'] is a special entry that represents the prevalence
    of unobserved taxa.

    When all taxa are accounted for, taxa['other'] is removed.
    """

    def __init__(self, n):
        self.n = n
        self.multiplier = n
        self.taxa = dict(other=Beta(yes=0, no=n-2))

    def GetTaxa(self):
        return self.taxa

    def Get(self, taxon):
        """Looks up the distribution for a given taxon (or None)."""
        return self.taxa.get(taxon)

    def Other(self):
        """Returns the distribution for 'other', or None."""
        return self.taxa.get('other')

    def GetMean(self, taxon):
        """Computes the mean value of p for a given taxon, or 0."""
        dist = self.Get(taxon)
        return dist.Mean() if dist else 0

    def Likelihood(self, evidence):
        """Computes the likelihood of seeing a given taxon."""
        taxon = evidence
        if taxon in self.taxa:
            return self.GetMean(taxon)
        else:
            # if we haven't see this taxon before, we have to consider
            # multiple copies of the 'other' dist
            return self.multiplier * self.GetMean('other')

    def AddTaxon(self, taxon):
        """The first time we see a new taxon, add it to the census."""
        
        # make a copy of the 'other' dist
        dist = copy.deepcopy(self.Other())
        if dist is None:
            raise Exception('All taxa are accounted for.')
        self.taxa[taxon] = dist 

        # reduce the multiplier for the 'other' dist;
        # if we just accounted for the last unseen taxon, delete 'other'
        self.multiplier -= 1
        if self.multiplier == 0:
            del self.taxa['other']

    def GenerateTaxon(self):
        """Choose a random taxon from this census."""
        pmf = self.GeneratePrevalence()
        return pmf.Random()

    def GeneratePrevalence(self):
        """Generates a sample census.

        Result is a PMF that maps from taxon to prevalence.

        Gets the probability for 'a' exactly right; the rest are
        approximate.

        """
        dist = self.Get('a')
        prob = dist.Random()

        rest = self.taxa.keys()
        rest.remove('a')

        pmf = self.GenerateConditionalPrevalence(rest, 1-prob)
        pmf.Incr('a', prob)
        return pmf

    def GenerateConditionalPrevalence(self, rest, fraction):
        """Generates a sample census.

        rest is a list of taxa.
        fraction is the total probability for all taxa in rest.

        Result is a PMF that maps from taxon to prevalence.

        Note: instead of generating each prevalance and then choosing
        the next from the conditional distribution, we choose the
        prevalences independently and then normalize.

        """
        pmf = Pmf.Pmf(name='sample census')
        for taxon in rest:
            dist = self.Get(taxon)
            if taxon == 'other':
                prob = sum(dist.Random() for i in xrange(self.multiplier))
            else:
                prob = dist.Random()
            pmf.Incr(taxon, prob)
        pmf.Normalize(fraction)
        return pmf

    def Update(self, evidence):
        """Updates each dist with new evidence.

        evidence is a taxon
        """
        seen_taxon = evidence

        # if we haven't seen this taxon before, add it
        if seen_taxon not in self.taxa:
            self.AddTaxon(seen_taxon)

        for taxon, dist in self.taxa.iteritems():
            if taxon == seen_taxon:
                dist.Update(1, 0)
            else:
                dist.Update(0, 1)

    def Plot(self):
        """
        """
        for taxon, dist in sorted(self.taxa.iteritems()):
            print taxon, bayes.CredibleInterval(dist, 90)

    def PlotPmfs(self):
        """Plots the PMFs for each distribution."""
        # TODO: add labels
        pmfs = self.taxa.values()
        PlotPmfs(pmfs)


class MetaHypo(object):
    """Represents a belief about the number of taxa in a population.

    Maps from a hypothesis to a probability.
    """
    def __init__(self, prior):
        """Makes a meta hypothesis with the given prior.

        prior maps from number of taxa to probability.
        """
        self.hypos = Pmf.Pmf()
        for n, prob in prior.Items():
            hypo = Census(n)
            self.hypos.Incr(hypo, prob)

    def GetHypos(self):
        return self.hypos

    def PlotHypos(self, taxon):
        """Looks up the PMFs for a given taxon and plots them."""
        pmfs = Pmf.Pmf()
        for hypo, prob in self.hypos.Items():
            if hypo.n > 1:
                dist = hypo.Get(taxon)
                pmfs.Incr(dist.Pmf(), prob)
        PlotMixture(pmfs)

    def Update(self, evidence):
        """Updates based on observing a given taxon."""
        for hypo, prob in self.hypos.Items():
            likelihood = hypo.Likelihood(evidence)
            if likelihood:
                self.hypos.Mult(hypo, likelihood)
            else:
                # if a hypothesis has been ruled out, remove it
                self.hypos.Remove(hypo)

        self.hypos.Normalize()

        # update the hypotheses
        for hypo, prob in self.hypos.Items():
            hypo.Update(evidence)

    def Print(self):
        """Prints the PMF for number of taxa."""
        self.Pmf().Print()

    def Pmf(self, name=''):
        """Returns the PMF for number of taxa."""
        t = [(hypo.n, prob) for hypo, prob in self.hypos.Items()]
        return Pmf.MakePmfFromDict(dict(t), name=name)

    def GenerateCurves(self, sample, m=15, iters=10):
        curves = []
        for i in range(iters):
            taxons, curve = self.GenerateRarefactionCurve(sample, m=m)
            print taxons
            curves.append(curve)
        return curves

    def GenerateRarefactionCurve(self, sample, m):
        """Generates a possible rarefaction curve.

        Starts with the given sample and generates m additional taxa.

        Returns the complete sample and the rarefaction curve.
        """
        def NextNameGenerator(skip):
            """Generates names for unobserved taxa."""
            for name in lowercase[skip:]:
                yield name

        names = NextNameGenerator(skip=len(set(sample)))

        meta = copy.deepcopy(self)
        taxons = list(sample)
        random.shuffle(taxons)

        for i in xrange(m):
            taxon = self.GenerateTaxon()
            if taxon == 'other':
                taxon = names.next()
            taxons.append(taxon)
            meta.Update(taxon)

        curve = MakeCurve(taxons)
        return ''.join(taxons), curve

    def GenerateTaxon(self):
        """Chooses a random taxon."""
        pmf = Pmf.Pmf()
        for hypo, prob in sorted(self.hypos.Items()):
            taxon = hypo.GenerateTaxon()
            pmf.Incr(taxon, prob)
        return pmf.Random()


def PlotPmfs(pmfs):
    myplot.Pmfs(pmfs, show=True,
                xlabel='Prevalence',
                ylabel='Probability',
                )

def PlotMixture(pmfs, show=False):
    for dist in pmfs.Values():
        xs, ys = dist.Render()
        pyplot.plot(xs, ys, color='blue', alpha=0.2)

    mix = Pmf.MakeMixture(pmfs)
    xs, ys = mix.Render()
    pyplot.plot(xs, ys, color='blue', alpha=0.9, linewidth=2)

    myplot.Plot(show=show, clf=False,
                xlabel='Prevalence',
                ylabel='Probability',
                legend=False)


def MakeCurve(sample):
    """Makes a rarefaction curve for the given sample."""
    s = set()
    curve = []
    for i, taxon in enumerate(sample):
        s.add(taxon)
        curve.append((i+1, len(s)))
    return curve


def JitterCurve(curve, dx=0.2, dy=0.3):
    """Adds random noise to the pairs in a curve.

    dx and dy control the amplitude of the noise in each dimension.
    """
    curve = [(x+random.uniform(-dx, dx), 
              y+random.uniform(-dy, dy)) for x, y in curve]
    return curve


def PlotCurves(curves, show=False):
    """Plots a set of curves.

    curves is a list of curves; each curve is a list of (x, y) pairs.
    """
    for curve in curves:
        curve = JitterCurve(curve)
        xs, ys = zip(*curve)
        pyplot.plot(xs, ys, color='blue', alpha=0.2)

    myplot.Plot(show=show,
                clf=False,
                #title='Rarefaction curve',
                xlabel='# samples',
                ylabel='# taxa',
                legend=False)


def ProbCurve(curves, n, m):
    ms = range(1, m+1)
    ps = []
    for m in ms:
        p = MoreTaxons(curves, n, m)
        ps.append(p)
    return ms, ps


def MoreTaxons(curves, start, m):
    """Estimates probability of finding at >=1 new taxon after m tries."""
    count = 0
    for curve in curves:
        _, already_found = curve[start-1]
        _, finally_found = curve[start+m-1]
        if finally_found > already_found:
            count += 1
    return float(count) / len(curves)


def MakeSample(freqs):
    """Make a sequence of letters that with the given frequencies.

    In descending order.
    """
    t = []
    freqs.sort(reverse=True)
    for taxon, freq in zip(lowercase, freqs):
        t.extend([taxon] * freq)
    return ''.join(t)


def MakeMeta(sample):
    prior = Pmf.MakePmfFromList(range(1, 21), name='prior')
    meta = MetaHypo(prior)
    for taxon in sample:
        meta.Update(taxon)

    return meta


def main(script, flag=1, *args):
    flag = int(flag)

    pyplot.rc('figure', figsize=(12, 8))

    d = {
        1: [9, 6],
        2: [6, 4, 1, 1, 1, 1],
        }
    freqs = d[flag]

    sample = MakeSample(freqs)
    meta = MakeMeta(sample)

    meta.Print()
    posterior = meta.Pmf(name='posterior')

    # plot posterior on # taxa
    pyplot.subplot(2, 2, 1)
    myplot.Pmfs([posterior],
                clf=False,
                xlabel='# of taxa',
                ylabel='Probability')

    # plot prevalence of 'a'
    pyplot.subplot(2, 2, 2)
    meta.PlotHypos(taxon='a')

    # generate curves
    pyplot.subplot(2, 2, 3)
    curves = meta.GenerateCurves(sample, m=15, iters=100)
    PlotCurves(curves)

    # plot prob of finding more taxa
    pyplot.subplot(2, 2, 4)

    ms, ps = ProbCurve(curves, n=len(sample), m=15)
    pyplot.plot(ms, ps)
    myplot.Plot(xlabel='# additional samples', ylabel='Prob of more taxa')

    #pyplot.show()

    myplot.Plot(root='rare%d'%flag)

if __name__ == '__main__':
    profile = False
    if profile:
        import cProfile
        cProfile.run('main(*sys.argv)')
    else:
        main(*sys.argv)
