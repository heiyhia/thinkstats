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

def RunTests(iters=1000):
    pool, firsts, others = cumulative.MakeTables()
    RunTest('length', pool.lengths, firsts.lengths, others.lengths, iters)
    

def RunTest(root, pool, actual1, actual2, iters=1000, trim=True):
    if trim:
        pool = thinkstats.Trim(pool)
        actual1 = thinkstats.Trim(actual1)
        actual2 = thinkstats.Trim(actual2)

    Test(root + 'deltas_cdf',
         actual1, actual2, pool, pool,
         iters)

    Test(root + 'deltas_ha_cdf',
         actual1, actual2, actual1, actual2,
         iters)

    
def Test(root, actual1, actual2, model1, model2, iters=1000, plot=False):
    n = len(actual1)
    m = len(actual2)
    
    mu1 = thinkstats.Mean(actual1)
    mu2 = thinkstats.Mean(actual2)
    
    delta = mu1 - mu2
    print n, m, mu1, mu2, delta
    
    delta = abs(delta)
    deltas = [Resample(model1, model2, n, m) for i in range(iters)]
    center = thinkstats.Mean(deltas)
    print 'Center', center

    cdf = Cdf.MakeCdfFromList(deltas)

    left = cdf.Prob(-delta)
    right = 1.0 - cdf.Prob(delta)
    
    print 'Tails', left, right, left+right

    if plot:
        PlotCdf(cdf, delta)
    
def PlotCdf(cdf, delta):
    def VertLine(x):
        xs = [x, x]
        ys = [0, 1]
        pyplot.plot(xs, ys, linewidth=2, color='0.7')
        
    VertLine(-delta)
    VertLine(delta)

    xs, ys = cdf.Render()    
    pyplot.plot(xs, ys)
    
    myplot.Plot(root,
                title='Resampled differences',
                xlabel='difference in weeks',
                ylabel='CDF(x)',
                legend=False) 
    

    
def Resample(t1, t2, n, m):
    sample1 = SampleWithReplacement(t1, n)
    sample2 = SampleWithReplacement(t2, m)

    mu1 = thinkstats.Mean(sample1)
    mu2 = thinkstats.Mean(sample2)
    
    delta = mu1 - mu2
    return delta

def SampleWithReplacement(t, n):    
    return [random.choice(t) for i in range(n)]

def SampleWithoutReplacement(t, n):
    return random.sample(t, n)
 
def main():
    RunTests(1000)
    
if __name__ == "__main__":
    main()
