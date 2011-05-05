"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import Cdf
import myplot
import random
import thinkstats


def MakeFigure():
    fp = open('babyboom.dat')
    
    # skip to the beginning of the data
    for line in fp:
        if line.find('START DATA') != -1:
            break
    
    # read a list of times
    times = []
    for line in fp:
        t = line.split()
        time = int(t[-1])
        times.append(time)
    
    # compute interarrival times
    diffs = [times[0]]
    for i in range(len(times)-1):
        diff = times[i+1] - times[i]
        diffs.append(diff)
    
    n = len(diffs)
    mu = thinkstats.Mean(diffs)
        
    print 'mean interarrival time', mu
    
    cdf = Cdf.MakeCdfFromList(diffs, 'actual')

    sample = [random.expovariate(1/mu) for i in range(n)]
    model = Cdf.MakeCdfFromList(sample, 'model')
    
    myplot.Cdf(cdf, root='interarrivals',
              title='Time between births',
              xlabel='minutes',
              ylabel='CDF',
              legend=False)

    myplot.Cdf(cdf, root='interarrivals_logy',
              complement=True,
              title='Time between births',
              xlabel='minutes',
              ylabel='Complementary CDF',
              yscale='log',
              legend=False)

    myplot.Cdfs([cdf, model], root='interarrivals_model',
              complement=True,
              title='Time between births',
              xlabel='minutes',
              ylabel='Complementary CDF',
              yscale='log')

def main():
    MakeFigure()
    
if __name__ == "__main__":
    main()
