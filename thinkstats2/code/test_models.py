"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2014 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import thinkstats2
import thinkplot


def read_file(filename):
    fp = open(filename)
    data = []
    for line in fp:
        x = float(line.strip())
        data.append(x)
    return data


def main():
    filename = 'mystery0.dat'
    data = read_file(filename)
    cdf = thinkstats2.MakeCdfFromList(data)

    thinkplot.SubPlot(2, 2, 1)
    thinkplot.Cdf(cdf)

    thinkplot.SubPlot(2, 2, 2)
    scale = thinkplot.Cdf(cdf, xscale='log')
    

    thinkplot.Show()


if __name__ == '__main__':
    main()
