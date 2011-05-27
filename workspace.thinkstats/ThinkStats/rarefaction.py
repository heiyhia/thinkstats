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
"""

"""


class Beta(object):
    def __init__(self, yes, no):
        self.alpha = yes+1
        self.beta = no+1

    def Update(self, yes, no):
        self.alpha += yes
        self.beta += no
        
    def Mean(self):
        """Computes the mean of this distribution."""
        return float(self.alpha) / (self.alpha + self.beta)

    def Random(self):
        """Generate a random variate from this distribution."""
        # handle the special case where there is one taxon, so
        # its prevalence is necessarily 1
        if self.beta == 0:
            return 1

        return random.betavariate(self.alpha, self.beta)

    def Render(self, steps=101):
        """Returns a curve that represents the PDF of this distribution."""
        ps = [i / (steps-1.0) for i in range(steps)]
        probs = [math.pow(p, self.alpha-1) * math.pow(1-p, self.beta-1) 
                 for p in ps]
        pmf = Pmf.MakePmfFromDict(dict(zip(ps, probs)))
        pmf.Normalize()


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
        self.n = n
        self.multiplier = n
        self.suites = {}
        self.suites['other'] = Beta(yes=0, no=n-2)

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

    def GeneratePrevalence(self):
        """Generates a possible census.

        Result is a PMF that maps from taxon to prevalence.
        """
        pmf = Pmf.Pmf(name='sample census')
        for taxon, suite in self.suites.iteritems():
            prob = suite.Random()
            if taxon == 'other':
                prob *= self.multiplier
            pmf.Incr(taxon, prob)
        pmf.Normalize()
        return pmf

    def GenerateTaxon(self):
        """Choose a random taxon from this census."""
        pmf = self.GeneratePrevalence()
        return pmf.Random()

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
            suite.Update(*yesno)

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

    def GenerateTaxon(self):
        """Chooses a random taxon."""
        pmf = Pmf.Pmf()
        for hypo, prob in sorted(self.hypos.Items()):
            taxon = hypo.GenerateTaxon()
            pmf.Incr(taxon, prob)
        return pmf.Random()

    def GenerateRarefactionCurve(self, sample, m):
        """Generates a possible rarefaction curve.

        Starts with the given sample and generates m additional taxa.

        Returns the complete sample and the rarefaction curve.
        """
        def NextNameGenerator(skip):
            """Generates names for unobserved taxa."""
            for name in string.ascii_lowercase[skip:]:
                yield name

        names = NextNameGenerator(skip=len(set(sample)))

        meta = copy.deepcopy(self)
        taxons = list(sample)
        random.shuffle(taxons)

        for i in range(m):
            taxon = self.GenerateTaxon()
            if taxon == 'other':
                taxon = names.next()
            taxons.append(taxon)
            meta.Update(taxon)

        curve = MakeCurve(taxons)
        return ''.join(taxons), curve


def PlotPmfs(pmfs):
    myplot.Pmfs(pmfs, show=True,
                xlabel='Prevalence',
                ylabel='Probability',
                )


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


def PlotCurves(curves, color='b'):
    """Plots a set of curves.

    curves is a list of curves; each curve is a list of (x, y) pairs.
    """
    for curve in curves:
        curve = JitterCurve(curve)
        xs, ys = zip(*curve)
        pyplot.plot(xs, ys, color)

    myplot.Plot(show=True,
                title='Rarefaction curve',
                xlabel='# samples',
                ylabel='# taxa',
                )


def main():
    sample = 'aaaaaaaaabbbbb'
    sample = 'aaaaaabbbbcdef'

    prior = Pmf.MakePmfFromList(range(1, 21))
    meta = MetaHypo(prior)
    for taxon in sample:
        meta.Update(taxon)

    curves = []
    iters = 50
    for i in range(iters):
        taxons, curve = meta.GenerateRarefactionCurve(sample, 15)
        print taxons
        curves.append(curve)

    PlotCurves(curves)
    meta.Print()


if __name__ == '__main__':
    main()
