"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import survey
import thinkstats

import Cdf
from Pmf import Pmf
import myplot


def Summarize(data_dir):
    """Prints summary statistics for first babies and others.
    
    Returns:
        tuple of Tables
    """
    table = survey.Pregnancies()
    table.ReadRecords(data_dir)

    # make a map from caseid to list of pairs
    d = {}
    for record in table.records:
        # skip non-live births
        if record.outcome != 1:
            continue

        # skip multiple births
        if record.nbrnaliv > 1:
            continue

        pair = record.birthord, record.prglength
        d.setdefault(record.caseid, []).append(pair)

    print len(d)

    # find all caseids with more than one live birth
    pmf = Pmf()
    for caseid, t in d.iteritems():
        if len(t) <= 1:
            continue

        t.sort()
        _, prglength1 = t[0] 
        _, prglength2 = t[1] 

        if prglength1 < 15 or prglength2 < 15:
            continue

        diff = prglength2 - prglength1
        if abs(diff) > 15:
            print caseid, prglength1, prglength2

        pmf.Incr(diff)

    pmf.Normalize()
    return pmf


def main(name, data_dir='.'):
    pmf = Summarize(data_dir)
    cdf = Cdf.MakeCdfFromPmf(pmf)

    myplot.Cdf(cdf)
    myplot.Save(root='first_matched',
                xlabel='Difference in weeks',
                ylabel='Cumulative probability',
                )

    mean = pmf.Mean()
    print mean, 'weeks'
    print mean * 7 * 24, 'hours'
    

if __name__ == '__main__':
    import sys
    main(*sys.argv)
