"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import bayes
import copy
import matplotlib.pyplot as pyplot
import myplot
import Pmf
import Cdf
import math
import random
"""

"""


def MakeUniformSuite(low, high, steps, name=None):
    """Makes a PMF that represents a suite of hypotheses with equal p.
    
    Args:
        low: low end of range
        high: high end of range
        steps: number of values

    Returns:
        Pmf object
    """
    hypos = [low + (high-low) * i / (steps-1.0) for i in range(steps)]
    pmf = Pmf.MakePmfFromList(hypos, name=name)
    return pmf


def MakeBinomialSuite(yes, no, steps=11):
    """Makes a PMF that represents a suite of hypotheses with equal p.
    
    Args:
        low: low end of range
        high: high end of range
        steps: number of values

    Returns:
        Pmf object
    """
    if no == -1:
        return Pmf.MakePmfFromList([1])

    ps = [i / (steps-1.0) for i in range(steps)]
    probs = [math.pow(p, yes) * math.pow(1-p, no) for p in ps]
    pmf = Pmf.MakePmfFromDict(dict(zip(ps, probs)))
    pmf.Normalize()
    return pmf


def MakeExponentialSuite(low, high, steps, lam=10.0):
    """Makes a PMF that represents a suite of hypotheses with equal p.
    
    Args:
        low: low end of range
        high: high end of range
        steps: number of values

    Returns:
        Pmf object
    """
    ps = [low + (high-low) * i / (steps-1.0) for i in range(steps)]
    probs = [math.exp(-lam*p) for p in ps]
    pmf = Pmf.MakePmfFromDict(dict(zip(ps, probs)))
    pmf.Normalize()
    return pmf


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

    def PlotHypos(self, taxon):
        """Looks up the PMFs for a given taxon and plots them."""
        pmfs = []
        for hypo in self.hypos.Values():
            if hypo.n > 1 and hypo.Get(taxon):
                pmfs.append(hypo.Get(taxon))
        PlotPmfs(pmfs)

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
        for hypo, prob in sorted(self.hypos.Items()):
            print hypo.n, prob


class Census(object):
    """Represents a belief about the population of a sample.

    n is the hypothetical number of taxa

    multiplier is the number of taxa that have not been observed.

    suites is a mapping from taxon to a suite of hypotheses
    about p, the fractional population of the taxon.
    
    suites['other'] is a special entry that represents the prevalence
    of unobserved taxa.

    When all taxa are accounted for, suites['other'] is removed.
    """

    def __init__(self, n):
        self.updater = bayes.Binomial()
        self.n = n
        self.multiplier = n
        self.suites = {}
        self.suites['other'] = MakeBinomialSuite(yes=0, no=n-2, steps=101)

    def Get(self, taxon):
        """Looks up the suite for a given taxon (or None)."""
        return self.suites.get(taxon)

    def Other(self):
        """Returns the suite for 'other', or None."""
        return self.suites.get('other')

    def GetMean(self, taxon):
        """Computes the mean value of p for a given taxon, or 0."""
        suite = self.Get(taxon)
        return suite.Mean() if suite else 0

    def Likelihood(self, evidence):
        """Computes the likelihood of seeing a given taxon."""
        taxon = evidence
        if taxon in self.suites:
            return self.GetMean(taxon)
        else:
            # if we haven't see this taxon before, we have to consider
            # multiple copies of the 'other' suite
            return self.multiplier * self.GetMean('other')

    def AddTaxon(self, taxon):
        """The first time we see a new taxon, add it to the census."""
        
        # make a copy of the 'other' suite
        suite = copy.deepcopy(self.Other())
        if suite is None:
            raise Exception('All taxa are accounted for.')
        self.suites[taxon] = suite 

        # reduce the multiplier for the 'other' suite
        self.multiplier -= 1

        # if we just accounted for the last unseen taxon, delete 'other'
        if self.multiplier == 0:
            del self.suites['other']

    def MakeCdfs(self):
        """Generates a CDF for each suite.

        Makes Generate() faster.
        """
        self.cdfs = {}
        for taxon, suite in self.suites.iteritems():
            self.cdfs[taxon] = Cdf.MakeCdfFromPmf(suite)

    def Generate(self):
        """Given the current CDFs, generate a possible census.

        Result is a PMF that maps from taxon to prevalence.
        """
        pmf = Pmf.Pmf(name='sample census')
        # TODO: go in descending order by mean
        for taxon, cdf in self.cdfs.iteritems():
            prob = cdf.Random()
            if taxon == 'other':
                prob *= self.multiplier
            pmf.Incr(taxon, prob)
        pmf.Normalize()
        return pmf

    def GenerateTaxon(self):
        """Choose a random taxon from this census."""
        pmf = self.Generate()
        cdf = Cdf.MakeCdfFromPmf(pmf)
        return cdf.Random()

    def Update(self, evidence):
        """Updates each suite with new evidence.

        evidence is a taxon
        """
        seen_taxon = evidence

        # if we haven't seen this taxon before, add it
        if seen_taxon not in self.suites:
            self.AddTaxon(seen_taxon)

        for taxon, suite in self.suites.iteritems():
            yesno = [1, 0] if taxon==seen_taxon else [0, 1]
            self.updater.Update(suite, yesno)

    def Plot(self):
        """
        """
        for taxon, suite in sorted(self.suites.iteritems()):
            print taxon, bayes.CredibleInterval(suite, 90)

    def PlotPmfs(self):
        """Plots the PMFs for each suite."""
        # TODO: add labels
        pmfs = self.suites.values()
        PlotPmfs(pmfs)

    def GenerateRarefactionCurves(self, n, m, s, iters=50):
        self.MakeCdfs()
        curves = [self.GenerateCurve(n, m, s) for i in range(iters)]
        return curves

    def GenerateCurve(self, n, m, s):
        pmf = self.Generate()
        cdf = Cdf.MakeCdfFromPmf(pmf)
        sample = cdf.Sample(m)
        return self.MakeCurve(sample, n, s)

    def PlotBeforeAndAfter(self, sample):
        curves1 = self.MakePossibleCurves(sample)
        
        n = len(sample)
        s = set(sample)
        curves2 = self.GenerateRarefactionCurves(n, 15, s)

        pyplot.clf()
        PlotCurves(curves1, 'b')
        PlotCurves(curves2, 'g')
        myplot.Plot(show=True,
                    title='Rarefaction curve',
                    xlabel='# samples',
                    ylabel='# taxa',
                    )


def PlotPmfs(pmfs):
    myplot.Pmfs(pmfs, show=True,
                xlabel='Prevalence',
                ylabel='Probability',
                )


def MakePossibleCurves(sample, iters=5):
    t = list(sample)
    curves = []
    for i in range(iters):
        random.shuffle(t)
        curve = MakeCurve(t)
        curves.append(curve)
    return curves


def MakeCurve(sample):
    s = set()
    curve = []
    for i, taxon in enumerate(sample):
        s.add(taxon)
        curve.append((i+1, len(s)))
    return curve


def PlotCurves(curves, color='b'):
    for curve in curves:
        curve = JitterCurve(curve)
        xs, ys = zip(*curve)
        pyplot.plot(xs, ys, color)


def JitterCurve(curve, dx=0.2, dy=0.3):
    curve = [(x+random.uniform(-dx, dx), 
              y+random.uniform(-dy, dy)) for x, y in curve]
    return curve


def main():
    sample = 'aaabbc'
    curves = MakePossibleCurves(sample)
    PlotCurves(curves)
    myplot.Plot(show=True,
                title='Rarefaction curve',
                xlabel='# samples',
                ylabel='# taxa',
                )

    prior = Pmf.MakePmfFromList(range(1,20))
    meta = MetaHypo(prior)
    for taxon in sample:
        meta.Update(taxon)

    meta.Print()
    #meta.PlotHypos('other')
    return

    MakeCurve('aab')
    return

    inventory = list('abcdefghijk')

    census = Census(inventory)

    sample = 'aaaaaaaaabbbbb'
    sample = 'aaaaaabbbbcdef'
    evidence = Pmf.MakeHistFromList(sample)
    census.Update(evidence)
    census.Plot()
    census.PlotPmfs()

    census.PlotBeforeAndAfter(sample)

if __name__ == '__main__':
    main()
