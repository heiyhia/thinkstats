"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import Cdf
import myplot
import random
import thinkstats
import matplotlib.pyplot as pyplot

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
    
    myplot.Cdf(cdf)
    myplot.Save(root='interarrivals',
              title='Time between births',
              xlabel='minutes',
              ylabel='CDF',
              legend=False,
                formats=['eps', 'png', 'pdf'])

    myplot.Cdfs([cdf, model], complement=True)
    myplot.Save(root='interarrivals_model',
                title='Time between births',
                xlabel='minutes',
                ylabel='Complementary CDF',
                yscale='log',
                formats=['eps', 'png', 'pdf'])

    pyplot.subplots_adjust(bottom=0.11)
    myplot.Cdf(cdf, complement=True)
    myplot.Save(root='interarrivals_logy',
                title='Time between births',
                xlabel='minutes',
                ylabel='Complementary CDF',
                yscale='log',
                legend=False,
                formats=['eps', 'png', 'pdf'])

def main():
    MakeFigure()
    
if __name__ == "__main__":
    main()
