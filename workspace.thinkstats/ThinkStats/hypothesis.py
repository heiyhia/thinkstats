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
    mean_var = thinkstats.MeanVar(pool.lengths)
    print 'Pooled mean var', mean_var
    
    RunTest('length', pool.lengths, firsts.lengths, others.lengths, 
            iters,
            partition=False)
    

def RunTest(root, pool, actual1, actual2, 
            iters=1000,
            trim=True,
            partition=False):
    """Computes the distributions of delta under H0 and HA.
    
    Args:
        root: string filename root for the plots
        
        pool: sequence of values from the pooled distribution
        
        actual1: sequence of values in group 1
        
        actual2: sequence of values in group 2
        
        iters: how many resamples
        
        trim: whether to trim the sequences
        
        partition: whether to cross-validate by partitioning the data
    """
    if trim:
        pool.sort()
        actual1.sort()
        actual2.sort()
        pool = thinkstats.Trim(pool)
        actual1 = thinkstats.Trim(actual1)
        actual2 = thinkstats.Trim(actual2)

    if partition:
        n = len(actual1)
        m = len(actual2)
        actual1, model1 = Partition(actual1, n/2)
        actual2, model2 = Partition(actual2, m/2)
        pool = model1 + model2
    else:
        model1 = actual1
        model2 = actual2
        
    ph0 = Test(root + '_deltas_cdf',
               actual1, actual2, pool, pool,
               iters, plot=True)

    pha = Test(root + '_deltas_ha_cdf',
               actual1, actual2, model1, model2,
               iters)

    prior = 0.5
    pe = prior*pha + (1-prior)*ph0
    posterior = prior * pha / pe
    print 'Posterior', posterior


def DifferenceInMean(actual1, actual2):
    mu1 = thinkstats.Mean(actual1)
    mu2 = thinkstats.Mean(actual2)
    delta = mu1 - mu2
    return mu1, mu2, delta



def Test(root, actual1, actual2, model1, model2, iters=1000, plot=False):
    """Estimates p-values based on differences in the mean.
    
    Args:
        root: string filename base for plots
        
        actual1:
        actual2: sequences of observed values for groups 1 and 2
        
        model1: 
        model2: sequences of values from the hypothetical distributions
        
        iters: how many resamples
        
        plot: whether to plot the distribution of differences in the mean
    """
    n = len(actual1)
    m = len(actual2)
    
    mu1, mu2, delta = DifferenceInMean(actual1, actual2)
    delta = abs(delta)

    cdf, pvalue = PValue(model1, model2, n, m, delta, iters)
    print n, m, mu1, mu2, delta, pvalue

    var_pooled = 7.3018637881954334
    f = 1.0 / n + 1.0 / m
    mu, var = (0, f * var_pooled)
    sigma = math.sqrt(var)
    print 'Expected mean var deltas', mu, var

    if plot:
        PlotCdf(root, cdf, delta)
        # PlotModel(root, cdf, mu, sigma)
        
    return pvalue
    
    
def PValue(model1, model2, n, m, delta, iters=1000, plot=False):
    deltas = [Resample(model1, model2, n, m) for i in range(iters)]
    mean_var = thinkstats.MeanVar(deltas)
    print 'Actual mean var deltas', mean_var

    cdf = Cdf.MakeCdfFromList(deltas)

    left = cdf.Prob(-delta)
    right = 1.0 - cdf.Prob(delta)
    
    pvalue = left + right
    print 'Tails', left, right
    print 'Pvalue', pvalue

    return cdf, pvalue


def PlotModel(root, cdf, mu, sigma):
    low, high = -3 * sigma, 3 * sigma
    steps = 100
    xs = [low + (high - low) * i / steps for i in range(steps)]
    ps = [erf.NormalCdf(x, mu, sigma) for x in xs]  

    pyplot.plot(xs, ps, linewidth=2, color='0.7')
     
    xs, ys = cdf.Render()    
    pyplot.plot(xs, ys)
    
    myplot.Plot(root,
                title='Resampled differences',
                xlabel='difference in weeks',
                ylabel='CDF(x)',
                legend=False) 

def PlotCdf(root, cdf, delta):
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
    """Computes the difference in mean of two samples.
    
    Args:
        t1: sequence of values
        t2: sequence of values
        
        n: size of the sample to draw from t1
        m: size of the sample to draw from t2
    """
    sample1 = SampleWithReplacement(t1, n)
    sample2 = SampleWithReplacement(t2, m)

    mu1 = thinkstats.Mean(sample1)
    mu2 = thinkstats.Mean(sample2)
    
    delta = mu1 - mu2
    return delta


def Partition(t, n):
    """Splits a sequences into two random partitions.
    
    Side effect: shuffles t
    
    Args:
        t: sequence of values
        
        n: size of the first partition

    Returns:
        two lists of values
    """
    random.shuffle(t)
    return t[:n], t[n:]


def SampleWithReplacement(t, n):
    """Generates a sample with replacement.
    
    Args:
        t: sequence of values
        
        n: size of the sample
        
    Returns:
        list of values
    """    
    return [random.choice(t) for i in range(n)]


def SampleWithoutReplacement(t, n):
    """Generates a sample without replacement.
    
    Args:
        t: sequence of values
        
        n: size of the sample
        
    Returns:
        list of values
    """    
    return random.sample(t, n)
 
 
def main():
    RunTests(1000)
    
if __name__ == "__main__":
    main()
