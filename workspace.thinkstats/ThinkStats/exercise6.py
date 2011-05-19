"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import erf
import Cdf
import cumulative
import math
import myplot
import random
import rankit
import thinkstats
import matplotlib.pyplot as pyplot


def RenderNormalCdf(mu, sigma, max, n=50):
    """Generates sequences of xs and ps for a normal CDF."""
    xs = [max * i / n for i in range(n)]    
    ps = [erf.NormalCdf(x, mu, sigma) for x in xs]
    return xs, ps


def MakeNormalModel(values):
    """Plot the CDF of birthweights with a normal model."""
    
    # estimate parameters: trimming outliers yields a better fit
    mu, var = thinkstats.TrimmedMeanVar(values, p=0.01)
    print 'Mean, Var', mu, var
    
    # plot the model
    sigma = math.sqrt(var)
    print 'Sigma', sigma
    xs, ps = RenderNormalCdf(mu, sigma, 200)

    pyplot.clf()
    pyplot.plot(xs, ps, label='model', linewidth=4, color='0.8')

    # plot the data
    cdf = Cdf.MakeCdfFromList(values)
    xs, ps = cdf.Render()
    pyplot.plot(xs, ps, label='data', linewidth=2, color='red')
 
    myplot.Plot(show=True,
                ylabel = 'CDF')


def GenerateNormalVariates(mu, sigma, n):
    """Generates a list of normal variates."""
    xs = [random.normalvariate(mu, sigma) for i in range(n)]
    return xs


def PlotSimulatedData(mu, sigma, n):
    """Generates a sample with the given parameters and plots it
    versus a list of normal variates."""
    ys = GenerateNormalVariates(mu, sigma, n)
    xs = GenerateNormalVariates(0, 1, n)
    pyplot.plot(sorted(xs), sorted(ys), '-', markersize=3, color='0.5')


def MakeNormalPlot(values):
    """Makes a normal probability plot.
    
    Args:
        values: sequence of values
        lineoptions: dictionary of options for pyplot.plot        
        options: dictionary of options for myplot.Plot
    """
    # TODO: when n is small, generate a larger sample and desample
    pyplot.clf()

    # compute parameters
    mu, var = thinkstats.TrimmedMeanVar(values, p=0.01)
    sigma = math.sqrt(var)
    n = len(values)

    # plot resampled data
    PlotSimulatedData(mu, sigma, n)

    # plot real data
    xs = GenerateNormalVariates(0, 1, n)
    pyplot.plot(sorted(xs), sorted(values), 'r.', markersize=3)
 
    myplot.Plot(show=True,
                xlabel='Standard normal values',
                legend=False)


def GetWeights(table):
    """Gets birth weight for each record in the table."""
    values = [p.totalwgt_oz for p in table.records
              if p.totalwgt_oz != 'NA']
    return values


def Jitter(values, jitter=0.5):
    """Jitters the values by adding a uniform variate in (-jitter, jitter)."""
    return [x + random.uniform(-jitter, jitter) for x in values]


def main():
    # test the distribution of birth weights for normality
    live_births, _, _ = cumulative.MakeTables()
    
    t = GetWeights(live_births)

    MakeNormalModel(t)

    
if __name__ == "__main__":
    main()
