"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import survey
import Pmf
import Cdf
import myplot


def ReadPregnancyRecords():
    """Reads survey data and returns a table of records."""
    table = survey.Pregnancies()
    table.ReadRecords()
    return table


def Validate(table):
    """Runs analysis on the given table.
    
    Args:
        table: table object
    """
    print len(table.records)
    count = 0
    for record in table.records:
        if record.prglength >= 27:
            count += 1
    print count


def MakeHistogram(table):
    lengths = [record.prglength for record in table.records]
    hist = Pmf.MakeHistFromList(lengths, name='pregnancy length')
    myplot.Hist(hist, show=True, xlabel='weeks', ylabel='count')


def MakePmf(table):
    lengths = [record.prglength for record in table.records]
    pmf = Pmf.MakePmfFromList(lengths, name='pregnancy length')
    myplot.Pmf(pmf, show=True, xlabel='weeks', ylabel='probability')


def MakeCdf(table):
    lengths = [record.prglength for record in table.records]
    cdf = Cdf.MakeCdfFromList(lengths, name='pregnancy length')
    myplot.Cdf(cdf, show=True, xlabel='weeks', ylabel='probability')


def main():
    table = ReadPregnancyRecords()
    Validate(table)
    

if __name__ == '__main__':
    main()
