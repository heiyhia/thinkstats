"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import thinkbayes

import matplotlib.pyplot as pyplot
import myplot
import numpy

import math
import random
import sys




def PlotMixture(pmfs, show=False):
    for dist in pmfs.Values():
        xs, ys = dist.Render()
        pyplot.plot(xs, ys, color='blue', alpha=0.2)

    mix = Pmf.MakeMixture(pmfs)
    xs, ys = mix.Render()
    pyplot.plot(xs, ys, color='blue', alpha=0.9, linewidth=2)

    myplot.Save(show=show, clf=False,
                xlabel='prevalence',
                ylabel='prob',
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


def OffsetCurve(curve, i, n, dx=0.3, dy=0.3):
    """Adds random noise to the pairs in a curve.

    i is the index of the curve
    n is the number of curves

    dx and dy control the amplitude of the noise in each dimension.
    """
    xoff = -dx + 2 * dx * i / (n-1)
    yoff = -dy + 2 * dy * i / (n-1)
    curve = [(x+xoff, y+yoff) for x, y in curve]
    return curve


def PlotCurves(curves, root=None, clf=False):
    """Plots a set of curves.

    curves is a list of curves; each curve is a list of (x, y) pairs.
    """
    if root: 
        pyplot.clf()

    n = len(curves)
    for i, curve in enumerate(curves):
        curve = OffsetCurve(curve, i, n)
        xs, ys = zip(*curve)
        pyplot.plot(xs, ys, color='blue', alpha=0.2)

    myplot.Save(root=root,
                clf=clf,
                xlabel='# samples',
                ylabel='# taxa',
                legend=False)


def ProbCurve(curves, n, m):
    """Computes the probability of finding more taxons.

    n is the number you've already seen
    m is the number of additional samples

    returns a list of ms and a list of probabilities
    """
    ms = range(1, m+1)
    ps = []
    for m in ms:
        p = MoreTaxons(curves, n, m)
        ps.append(p)

        print m, p
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


def MakeMeta(sample, n=31):
    prior = Pmf.MakePmfFromList(range(2, n), name='prior')
    meta = MetaHypo(prior)
    for taxon in sample:
        meta.Update(taxon)

    return meta


def MakeSubplots(root, sample, meta, m=15, iters=100):
    pyplot.figure(1, figsize=(12, 8))

    # plot posterior on # taxa
    pyplot.subplot(2, 2, 1)
    meta.PlotPosterior()

    # plot prevalence of 'a'
    pyplot.subplot(2, 2, 2)
    meta.PlotPrevalence()

    # generate curves
    pyplot.subplot(2, 2, 3)
    curves = meta.GenerateCurves(sample, m=m, iters=iters)
    PlotCurves(curves)

    # plot prob of finding more taxa
    pyplot.subplot(2, 2, 4)

    ms, ps = ProbCurve(curves, n=len(sample), m=m)
    pyplot.plot(ms, ps)
    myplot.Save(xlabel='# samples', 
                ylabel='prob of more taxa')

    #pyplot.show()

    myplot.Save(root=root)


def MakeFigures(root, sample, meta, m=15, iters=100):
    # plot posterior on # taxa
    meta.PlotPosteriorPmf(root + '.1')
    return
    # plot prevalence of 'a'
    meta.PlotPrevalence(root + '.2')

    # generate curves
    curves = meta.GenerateCurves(sample, m=m, iters=iters)
    PlotCurves(curves, root + '.3')

    # plot prob of finding more taxa
    ms, ps = ProbCurve(curves, n=len(sample), m=m)
    myplot.Plot(ms, ps,
                root=root + '.4',
                xlabel='# samples', 
                ylabel='prob of more taxa')


class Species(thinkbayes.Suite):
    
    def Update(self, data):
        thinkbayes.Suite.Update(self)
        for hypo in self.Values():
            hypo.Update(data)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under this hypothesis.

        hypo: Dirichlet object
        data: sequence of frequencies
        """
        dirichlet = hypo
        like = 0
        for i in range(1000):
            like += dirichlet.Likelihood(data)

        return like

    def DistOfN(self):
        """Computes the distribution of n."""
        pmf = thinkbayes.Pmf()
        for hypo, prob in self.Items():
            pmf.Set(hypo.n, prob)
        return pmf
        

def SimpleDirichletExample():
    beta = thinkbayes.Beta()
    beta.Update((3, 3))
    print beta.Mean()

    data = [3, 2, 1]
    dirichlet = thinkbayes.Dirichlet(3)
    dirichlet.Update(data)

    names = ['lions',  'tigers', 'bears']

    for i in range(3):
        beta = dirichlet.MarginalBeta(i)
        print names[i], beta.Mean()

        pmf = beta.MakePmf(name=names[i])
        print names[i], pmf.MaximumLikelihood()
        myplot.Pmf(pmf)

    myplot.Save(root='species1',
                xlabel='Prevalence',
                ylabel='Prob')





def main(script, flag=1, *args):
    SimpleDirichletExample()
    return

    data = [1, 1]

    hypos = [thinkbayes.Dirichlet(i) for i in range(len(data), 10)]
    suite = Species(hypos)

    suite.Update(data)

    pmf = suite.DistOfN()
    myplot.Pmf(pmf)
    myplot.Show()
    return
    
    random.seed(17)
    p = dirichlet.Likelihood(data)
    print p

    random.seed(17)
    lp = dirichlet.LogLikelihood(data)
    print lp, math.exp(lp)

    return


    xs = []
    for i in range(10000):
        x = dirichlet.Random()
        xs.append(x)

    cols = zip(*xs)
    for col in cols:
        cdf = thinkbayes.MakeCdfFromList(col)
        myplot.Cdf(cdf)

    myplot.Show()
    return


    beta = thinkbayes.Beta()
    beta.Update((140, 110))
    print beta.Mean()

    pmf = beta.MakePmf(101)
    myplot.Pmf(pmf)
    myplot.Show()
    return


    flag = int(flag)

    d = {
        0: [1],
        1: [9, 6],
        2: [8, 3],
        3: [11, 4, 2],
        4: [3, 2, 2, 1],
        5: [7, 3, 1, 1, 1, 1],
        6: [7, 2, 1],
        7: [4, 3, 9, 1, 1],
        8: [14, 5, 2, 2, 1,1,1,1,1,1],
        9: [11, 2, 1,1],
        }
    freqs = d[flag]

    sample = MakeSample(freqs)
    meta = MakeMeta(sample, n=31)

    meta.Print()

    root='rare%d' % flag
    MakeFigures(root, sample, meta)


if __name__ == '__main__':
    random.seed(17)

    profile = False
    if profile:
        import cProfile
        cProfile.run('main(*sys.argv)')
    else:
        main(*sys.argv)
