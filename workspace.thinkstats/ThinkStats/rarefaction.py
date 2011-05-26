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


def MakeLogisticSuite(low, high, steps):
    """Makes a PMF that represents a suite of hypotheses with equal p.
    
    Args:
        low: low end of range
        high: high end of range
        steps: number of values

    Returns:
        Pmf object
    """
    lps = [low + (high-low) * i / (steps-1.0) for i in range(steps)]
    ps = [bayes.Logistic(lp) for lp in lps]
    pmf = Pmf.MakePmfFromList(ps)
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


class Census:
    """Represents a belief about the population of a sample.

    A Census is a mapping from species to a suite of hypotheses
    about p, the fractional population of the species."""

    def __init__(self, inventory):
        self.updater = bayes.Binomial()

        self.suites = {}
        for name in inventory:
            self.suites[name] = MakeUniformSuite(0, 1, 101, name=name)

    def AddName(self, name):
        """The first time we see a new species, add it to the census."""
        if name not in self.suites:
            # if we haven't seen this new species in n tries,
            # we should initialize it with the oh-for-n distribution,
            # which is what the 'other' distribution is.
            suite = copy.deepcopy(self.suites['other'])
            suite.name = name
            self.suites[name] = suite 

    def MakeCdfs(self):
        self.cdfs = {}
        for name, suite in self.suites.iteritems():
            self.cdfs[name] = Cdf.MakeCdfFromPmf(suite)

    def Generate(self):
        pmf = Pmf.Pmf(name='sample census')
        # TODO: go in descending order by mean
        for name, cdf in self.cdfs.iteritems():
            prob = cdf.Random()
            pmf.Incr(name, prob)
        pmf.Normalize()
        return pmf

    def Update(self, evidence):
        
        total = evidence.Total()
        for name, suite in self.suites.iteritems():
            yes = evidence.Freq(name)
            no = total-yes
            self.updater.Update(suite, [yes, no])

    def Plot(self):
        for name, suite in sorted(self.suites.iteritems()):
            print name, bayes.CredibleInterval(suite, 90)

    def PlotPmfs(self):
        pmfs = self.suites.values()
        myplot.Pmfs(pmfs, show=True,
                    #root='rarefaction1',
                    title='Posterior pmfs',
                    xlabel='Prevalence',
                    ylabel='Probability',
                    )

    def GenerateRarefactionCurves(self, n, m, s, iters=50):
        self.MakeCdfs()
        curves = [self.GenerateCurve(n, m, s) for i in range(iters)]
        return curves

    def GenerateCurve(self, n, m, s):
        pmf = self.Generate()
        cdf = Cdf.MakeCdfFromPmf(pmf)
        sample = cdf.Sample(m)
        return self.MakeCurve(sample, n, s)

    def MakePossibleCurves(self, sample, iters=5):
        t = list(sample)
        curves = []
        for i in range(iters):
            random.shuffle(t)
            curve = self.MakeCurve(t)
            curves.append(curve)
        return curves

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
                    ylabel='# species',
                    )


def MakeCurve(sample):
    inventory = [sample[0], 'other']
    census = Census(inventory)

    for name in sample:
        census.AddName(name)

        pmf = Pmf.MakeHistFromList([name])
        census.Update(pmf)
        
    census.PlotPmfs()


def PlotCurves(curves, color):
    for curve in curves:
        curve = JitterCurve(curve)
        xs, ys = zip(*curve)
        pyplot.plot(xs, ys, color)


def JitterCurve(curve, dx=0.2, dy=0.3):
    curve = [(x+random.uniform(-dx, dx), 
              y+random.uniform(-dy, dy)) for x, y in curve]
    return curve


def main():
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
