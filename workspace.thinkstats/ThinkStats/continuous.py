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

def ExpoCdf(x, lam):
    return 1 - math.exp(-lam * x)

def ParetoCdf(x, alpha, xmin):
    if x < xmin:
        return 0
    return 1 - pow(x / xmin, -alpha)

def ParetoMedian(xmin, alpha):
    return xmin * pow(2, 1/alpha)

def MakeExpoCdf():
    n = 40
    max = 2.5
    xs = [max*i/n for i in range(n)]
    
    lam = 2.0
    ps = [ExpoCdf(x, lam) for x in xs]
    
    percentile = -math.log(0.05) / lam
    print 'Fraction <= ', percentile, ExpoCdf(lam, percentile)

    pyplot.clf()
    pyplot.plot(xs, ps)
    myplot.Plot('expo_cdf',
              title = 'Exponential CDF',
              xlabel = 'x',
              ylabel = 'CDF',
              legend=False)
    
def MakeParetoCdf():
    n = 50
    max = 10.0
    xs = [max*i/n for i in range(n)]
    
    xmin = 0.5
    alpha = 1.0
    ps = [ParetoCdf(x, alpha, xmin) for x in xs]
    print 'Fraction <= 10', ParetoCdf(xmin, alpha, 10)
    
    pyplot.clf()
    pyplot.plot(xs, ps)
    myplot.Plot('pareto_cdf',
              title = 'Pareto CDF',
              xlabel = 'x',
              ylabel = 'CDF',
              legend=False)
    
def MakeParetoCdf2():
    n = 50
    max = 1000.0
    xs = [max*i/n for i in range(n)]
    
    xmin = 100
    alpha = 1.7
    ps = [ParetoCdf(x, alpha, xmin) for x in xs]
    print 'Median', ParetoMedian(xmin, alpha)
    
    pyplot.clf()
    pyplot.plot(xs, ps)
    myplot.Plot('pareto_height',
              title = 'Pareto CDF',
              xlabel = 'height (cm)',
              ylabel = 'CDF',
              legend=False)
    


def RenderNormalCdf(mu, sigma, max, n=50):
    xs = [max * i / n for i in range(n)]    
    ps = [erf.NormalCdf(x, mu, sigma) for x in xs]
    return xs, ps


def MakeNormalCdf():
    xs, ps = RenderNormalCdf(2.0, 0.5, 4.0)
    
    pyplot.clf()
    pyplot.plot(xs, ps)
    myplot.Plot('normal_cdf',
              title = 'Normal CDF',
              xlabel = 'x',
              ylabel = 'CDF',
              legend=False)
    
    
def MakeNormalModel():
    pool, _, _ = cumulative.MakeTables()
    
    t = pool.weights[:]
    t.sort()
    rankit.MakeNormalPlot(t, 'nsfg_birthwgt_normal',
                          ylabel='Birth weights (oz)',)
    
    mu, var = thinkstats.TrimmedMeanVar(t)
    print 'Mean, Var', mu, var
    
    pyplot.clf()

    sigma = math.sqrt(var)
    print 'Sigma', sigma
    xs, ps = RenderNormalCdf(mu, sigma, 200)
    pyplot.plot(xs, ps, label='model', linewidth=3, color='0.7')

    xs, ps = pool.weight_cdf.Render()
    pyplot.plot(xs, ps, label='data', color='0.0')
 
    myplot.Plot('nsfg_birthwgt_model',
              title = 'Birth weights',
              xlabel = 'birth weight (oz)',
              ylabel = 'CDF')
    


 
def main():
    MakeNormalModel()
    return

    MakeExpoCdf()
    MakeParetoCdf()
    MakeParetoCdf2()
    MakeNormalCdf()
    
if __name__ == "__main__":
    main()
