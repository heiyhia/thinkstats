"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import matplotlib
import matplotlib.pyplot as pyplot
import myplot
import numpy
import random
import sys
import thinkstats
    
def ReadCsv(filename='nlsy79.csv'):
    fp = open(filename)
    reader = csv.reader(fp)
    header = reader.next()
    
    rows = []
    for t in reader:
        rows.append(t)
    
    fp.close()
    return rows

def SummarizeWeight(rows, input_limit=None):
    years = [1981, 1982, 1985, 1986, 1988, 1989, 1990, 1992, 
             1993, 1994, 1996, 1998, 2000, 2002, 2004, 2006, 2008]
    
    all_diffs = []
    
    for i, row in enumerate(rows):
        if i == input_limit:
            break
        
        id, race, sex = row[:3]
        weights = row[3:]
        print id
        diffs = Differences(years, weights, jitter=3)
        all_diffs.extend(diffs)
        
    weights, changes = zip(*all_diffs)
    
    print 'Mean weight', thinkstats.Mean(weights)
    print 'Mean change', thinkstats.Mean(changes)
    print numpy.corrcoef(weights, changes)

    pyplot.hexbin(weights, changes, cmap=matplotlib.cm.gray_r)
    myplot.Plot('nlsy_scatter',
                title = 'Weight change vs. weight',
                xlabel = 'Current weight (pounds)',
                ylabel = 'Weight change (pounds)',
                axis = [70, 270, -25, 25],
                legend=False,
                show=True,
                )
    
def Differences(years, weights, jitter=0.0):
    pairs = [(year, float(weight)) for year, weight in zip(years, weights)
             if float(weight) > 0]
    
    diffs = []
    for i in range(len(pairs)-1):
        y1, w1 = pairs[i]
        if w1 > 400:
            continue

        y2, w2 = pairs[i+1]
        years = y2-y1
        if years > 5:
            continue
        
        weight = w1 + random.uniform(-jitter, jitter)
        pounds = (w2 - w1) + random.uniform(-jitter, jitter)
        if abs(pounds) > 100:
            continue
        
        change = pounds / years
        diffs.append((weight, change))
        
    return diffs

def main(name):
    rows = ReadCsv()

    SummarizeWeight(rows)
    
if __name__ == '__main__':
    main(*sys.argv)
