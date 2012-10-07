"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib.pyplot as pyplot
import myplot
import csv

import Cdf
import correlation
import heri
import math
import random
import thinkstats


def ReadData(filename='heri.csv'):
    """Reads a CSV file of data from HERI's CIRP survey.

    Args:
      filename: string filename

    Returns:
      list of (score, number) pairs
    """
    fp = open(filename)
    reader = csv.reader(fp)
    res = []

    for t in reader:
        try:
            year = int(t[0])
            res.append(t)
        except ValueError:
            pass
    return res


def GetColumn(data, index):
    res = {}
    for row in data:
        try:
            year = int(row[0])
            res[year] = float(row[index]) / 10.0
        except ValueError:
            pass
    return res


def RenderColumn(col):
    return zip(*sorted(col.items()))


def DiffColumns(col1, col2):
    years1 = set(col1)
    years2 = set(col2)
    res = [(year, col1[year] - col2[year])for year in sorted(years1 & years2)]
    return zip(*res)


def MakePlot(filename='heri.csv'):
    """Generates a plot with the data, a fitted model, and error bars."""
    pyplot.clf()

    data = ReadData(filename)

    attended = GetColumn(data, 4)
    del attended[1966]
    ts, ys = RenderColumn(attended)
    ys = [100-y for y in ys]
    pyplot.plot(ts, ys, 'go-', linewidth=3, markersize=0, alpha=0.7,
                label='No attendance')

    nones = GetColumn(data, 1)
    ts, ys = RenderColumn(nones)
    pyplot.plot(ts, ys, 'bs-', linewidth=3, markersize=0, alpha=0.7,
                label='No religion')

    myplot.Save(root='heri2.3',
                title='',
                xlabel='',
                ylabel='Percent',
                axis=[1966, 2013, 0, 30])


def MakeGenderPlot(filename='heri.csv'):
    """Generates a plot with the data, a fitted model, and error bars."""
    pyplot.clf()

    data = ReadData(filename)

    men = GetColumn(data, 6)
    ts, ys = RenderColumn(men)
    pyplot.plot(ts, ys, 'b-', linewidth=3, alpha=0.7, label='men')

    women = GetColumn(data, 11)
    ts, ys = RenderColumn(women)
    pyplot.plot(ts, ys, 'g-', linewidth=3, alpha=0.7, label='women')

    myplot.Save(root='heri2.1',
                title='',
                xlabel='',
                ylabel='Preferred religion None (%)',
                axis=[1967, 2013, 0, 28])

    del men[1969]
    del women[1969]
    ts, ds = DiffColumns(men, women)

    heri.MakePlot(ts, ds,
                  model='ys ~ ts')

    pyplot.plot(ts, ds, color='purple', linewidth=3, alpha=0.7,
                label='Gender gap')

    myplot.Save(root='heri2.2',
                title='',
                xlabel='',
                ylabel='Percentage points',
                axis=[1967, 2013, 0, 6])

def main(script):
    MakePlot()
    MakeGenderPlot()


if __name__ == '__main__':
    import sys
    main(*sys.argv)

