"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

"""
How many students would you have to see to have an 80% chance of
seeing a stat.sig. difference between first-year students and juniors?

Race   Before   After   Diff    Sigma (before)  Sample size
White  1170     1211            179              118
Black   995     1001            167             4790


"""

import math

import Cdf
import Pmf
import myplot

import matplotlib.pyplot as pyplot

def NormalPdf(x):
    """Computes the PDF of x in the standard normal distribution."""
    return math.exp(-x**2/2) / math.sqrt(2 * math.pi)

def Frange(low, high, n):
    return [low + (high-low) * float(i)/(n-1) for i in range(n)]

def RenderPdf(mu, sigma, n=1001):
    xs = Frange(mu-4*sigma, mu+4*sigma, n)
    ys = [NormalPdf((x-mu) / sigma) for x in xs]
    return xs, ys

def MakeFigure():
    frac1 = 0.8
    frac2 = 1 - frac1

    xs, ys = RenderPdf(1170, 179)
    pmf1 = Pmf.MakePmfFromDict(dict(zip(xs, ys)), name='blue')

    xs, ys = RenderPdf(995, 167)
    pmf2 = Pmf.MakePmfFromDict(dict(zip(xs, ys)), name='green')

    myplot.Pmfs([pmf1, pmf2],
                root='normal1',
                xlabel='CLA score',
                ylabel='PDF',
                )

    pmf1.Normalize(frac1)
    pmf2.Normalize(frac2)

    ymax = max(pmf1.MaxLike(), pmf2.MaxLike())
    ymax = 0.003

    pyplot.clf()

    threshes = [1200, 1300, 1400, 1500, 1570]
    for thresh in threshes:
        myplot.Plot([thresh, thresh], [0, ymax],
                    clf=False,
                    line_options=dict(color='gray', alpha=0.5, linewidth=1))

    plot_options=[dict(color='blue', linewidth=2),
                  dict(color='green', linewidth=2)]

    myplot.Pmfs([pmf1, pmf2],
                plot_options=plot_options,
                clf=False,
                root='normal2',
                xlabel='CLA score',
                ylabel='PDF',
                )

    cdf1 = Cdf.MakeCdfFromPmf(pmf1)
    cdf2 = Cdf.MakeCdfFromPmf(pmf2)

    for thresh in threshes:
        p1 = frac1 * (1 - cdf1.Prob(thresh))
        p2 = frac2 * (1 - cdf2.Prob(thresh))

        den = p1 + p2
        rep1 = p1 / den
        rep2 = p2 / den
        print thresh, den, rep1, rep2

    return


    myplot.Cdfs([cdf1, cdf2],
                root='normal2',
                xlabel='',
                ylabel='',
                title='')

MakeFigure()
