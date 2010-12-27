"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import random

import Cdf
import myplot
import Pmf
import survey
import thinkstats

def MakeTable():
    """Reads survey data and returns table."""
    table = survey.Pregnancies()
    table.ReadRecords()

    sequences = {}
    d={1:'B', 2:'G'}

    for record in table.records:
        num, gender = record.nbrnaliv, record.babysex

        if num==1 and gender in [1,2]:
            sequences.setdefault(record.caseid, []).append(d[gender])

    return table, sequences


def MakeSequence(num=2, p=0.508, d={True:'B', False:'G'}):
    t = [d[random.random()<p] for i in range(num)]
    return t


def MakeSequences(seqs=1621):
    t = [MakeSequence() for i in range(seqs)]
    return t


def DiagDiff(hist):
    """Difference between the sum of the diagonal and off-diagonal elements."""
    x = hist.Freq('BB') + hist.Freq('GG') - hist.Freq('BG') - hist.Freq('GB')
    return x


def PrintGrid(hist):
    print hist.Freq('BB'), hist.Freq('BG')
    print hist.Freq('GB'), hist.Freq('GG')


def ProbSequence(sequence, p=0.508):
    d = dict(B=math.log(p), G=math.log(1-p))
    t = [d[x] for x in sequence]
    prob = math.exp(sum(t))
    return prob
    

def SubtractExpected(hist):
    total = hist.Total()
    diff = hist.Copy('diff')

    for val, freq in hist.Items():
        expected = total * ProbSequence(val)
        diff.Incr(val, -expected)

    return diff


def SummarizeSequences(sequences, trim=True):
    pairs = Pmf.Hist()
    for t in sequences:

        # use the first two babies
        if trim:
            t = t[0:2]

        if len(t) != 2:
            continue
        s = ''.join(t)
        pairs.Incr(s)
    return pairs


def SimulateDiags(num=1000):
    diags = []
    for i in range(num):
        fake = MakeSequences()
        pairs = SummarizeSequences(fake)
        diag = DiagDiff(pairs)
        diags.append(diag)

    print thinkstats.Mean(diags)

    cdf = Cdf.MakeCdfFromList(diags, 'diagonal bias')
    myplot.Cdf(cdf, show=True)

    return diags, cdf

def main():
    table, sequences = MakeTable()

    # get the sex ratio (boys per hundred girls)
    genders = Pmf.Hist()
    for t in sequences.values():
        for x in t:
            genders.Incr(x)

    boys, girls = genders.Freq('B'), genders.Freq('G')
    print boys, girls, 100.0 * boys / girls
    print boys, girls, 1.0 * boys / (boys+girls)

    pairs = SummarizeSequences(sequences.values())
    PrintGrid(pairs)

    diff = SubtractExpected(pairs)
    PrintGrid(diff)

    diag = DiagDiff(pairs)
    print pairs.Total(), diag

    return
    diags, cdf = SimulateDiags(1000)


if __name__ == '__main__':
    main()
