"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import first
import descriptive
import math
import Cdf
import Pmf
import survey
import thinkstats

import matplotlib.pyplot as pyplot
import plot


def Process(table, name):
    """Runs various analyses on this table."""
    descriptive.Process(table, name)
    GetWeights(table)
    print table.weight_cdf.name, len(table.weights)


def GetWeights(table):
    weights = []
    for record in table.records:
        lbs = record.birthwgt_lb
        oz = record.birthwgt_oz

        if lbs > 20 or oz > 15:
            continue

        if lbs == 'NA' or oz == 'NA':
            continue

        weight_oz = lbs * 16 + oz


        weights.append(weight_oz)

    table.weights = weights
    table.weight_pmf = Pmf.MakePmfFromList(weights, table.name)
    table.weight_cdf = Cdf.MakeCdfFromList(weights, table.name)


def MakeTables():
    """Reads survey data and returns a tuple of Tables"""
    table, firsts, others = first.MakeTables()
    pool = descriptive.PoolRecords(firsts, others)

    Process(pool, 'live births')
    Process(firsts, 'first babies')
    Process(others, 'others')
        
    return pool, firsts, others


def Summarize():
    pool, firsts, others = MakeTables()
    MakeFigures(firsts, others)


def MakeFigures(firsts, others):

    axis = None
    plot.Hists([firsts.weight_pmf, others.weight_pmf], 
               'nsfg_birthwgt_pmf', 
               title='Birth weight PMF',
               xlabel='ounces',
               ylabel='probability',
               axis=axis)

    plot.Cdfs([firsts.weight_cdf, others.weight_cdf], 
               'nsfg_birthwgt_cdf', 
               title='Birth weight CDF',
               xlabel='ounces',
               ylabel='probability',
               axis=axis)

def main():
    Summarize()
    

if __name__ == '__main__':
    main()
