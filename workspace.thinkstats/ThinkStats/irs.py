"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import sys

import myplot
import Pmf
import Cdf


def ReadIncomeFile(filename='08in11si.csv'):
    """Reads a data file from the IRS and returns the first two columns.

    Skips the header and returns only the first table (non-cumulative).

    Args:
      filename: string data file

    Returns:
      list of string pairs
    """
    reader = csv.reader(open(filename))
    for line in reader:
        if line[0] == 'All returns':
            break

    t = []
    for line in reader:
        if line[0].startswith('Accumulated'):
            break
        t.append(line[0:2])

    return t


def MakeIncomePdf(data):
    """Converts the strings from the IRS file to a PMF.

    Args:
      data: list of (dollar range, number of returns) string pairs

    Returns:
      tuple of (Hist, Pmf, Cdf) representing the number of returns in each bin
    """
    def clean(s):
        """Converts dollar amounts to integers."""
        try:
            return int(s.lstrip('$'))
        except ValueError:
            if s in ['No', 'income']:
                return 0
            if s == 'more':
                return -1
            return None

    def midpoint(low, high):
        """Finds the midpoint of a range."""
        if high == -1:
            return low * 3 / 2
        else:
            return (low + high) / 2

    hist = Pmf.Hist()

    for column, number in data:
        # convert the number of returns
        number = number.replace(',', '')
        number = int(number)

        # convert the income range
        column = column.replace(',', '')
        t = column.split()
        low, high = t[0], t[-1]
        low, high = clean(low), clean(high)

        # add to the histogram
        x = midpoint(low, high)
        hist.Incr(x, number)
        print x, number

    pmf = Pmf.MakePmfFromHist(hist)
    cdf = Cdf.MakeCdfFromDict(pmf.GetDict())
    return hist, pmf, cdf


def main(script, *args):
    data = ReadIncomeFile()
    hist, pmf, cdf = MakeIncomePdf(data)

    # plot the complementary CDF on a log-log scale
    myplot.Cdf(cdf, 
               complement=True,
               xscale='log',
               yscale='log',
               show=True)


if __name__ == "__main__":
    main(*sys.argv)
