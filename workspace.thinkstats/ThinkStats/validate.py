"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import survey
import Pmf
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


def main():
    table = ReadPregnancyRecords()
    Validate(table)
    MakeHistogram(table)
    

if __name__ == '__main__':
    main()
