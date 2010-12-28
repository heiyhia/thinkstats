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
    """Reads survey data and returns table.

    Returns:
        tuple of (Pregnancies table, list of sequences)
    """
    table = survey.Pregnancies()
    table.ReadRecords()

    sequences = {}
    d={1:'B', 2:'G'}

    for record in table.records:
        num, gender = record.nbrnaliv, record.babysex

        if num==1 and gender in [1,2]:
            sequences.setdefault(record.caseid, []).append(d[gender])

    return table, sequences.values()


def MakeSequence(num=2, p=0.508, d={True:'B', False:'G'}):
    """Makes a random sequence with length num."""
    t = [d[random.random()<p] for i in range(num)]
    return t


def MakeSequences(seqs=1621):
    """Makes a list of random sequences."""
    t = [MakeSequence() for i in range(seqs)]
    return t


def DiagDiff(hist):
    """Difference between the sum of the diagonal and off-diagonal elements."""
    x = hist.Freq('BB') + hist.Freq('GG') - hist.Freq('BG') - hist.Freq('GB')
    return x


def PrintGrid(hist):
    """Prints a 2x2 grid of frequencies for two-baby pairs.

    Args:
        hist: Hist object that maps sequences to frequencies.
    """
    print hist.Freq('BB'), hist.Freq('BG')
    print hist.Freq('GB'), hist.Freq('GG')


def ProbSequence(sequence, p=0.508):
    """Computes the probability of a sequence under the simple model.

    Assumes that the probability of a boy for all births is p.
    """
    d = dict(B=math.log(p), G=math.log(1-p))
    t = [d[x] for x in sequence]
    prob = math.exp(sum(t))
    return prob
    

def ProbSequenceSigma(sequence, probs):
    """Computes the probability of a sequence under the complex model.

    Assumes that the probability of a boy is described by the given
    distribution.
    """
    total = 0.0
    for val, prob in probs.Items():
        total += prob * ProbSequence(sequence, val)
    return total
    

def LikelihoodSequences(hist, sigma):
    """Computes the likelihood of a list of sequences given sigma.

    Args:
        hist: Histogram that maps sequences to their frequencies
        sigma: float std of proclivities

    Returns:
        float likelihood
    """
    probs = ComputeProbs(sigma=sigma)

    total = 1.0
    for sequence, freq in hist.Items():
        prob = ProbSequenceSigma(sequence, probs)
        print sequence, prob
        total *= prob
    return total


def ComputeProbs(p=0.508, sigma=0.5, n=101):
    """Make a Pmf of probabilities with the given parameter.

    The probabilities are the logistic tranform of the normal values.

    Args:
        p: float mean probability.
        sigma: float standard deviation of the normal values.
        n: number of values in the Pmf

    Returns:
        Pmf object
    """
    bias = p - 0.5
    low = -3.0
    high = 3.0
    vs = [low + i * (high-low) / (n-1) for i in range(n)]

    pmf = Pmf.Pmf(name='proclivities')
    for v in vs:
        p = 1 / (1 + math.exp(-v * sigma)) + bias
        prob = NormalPdf(v)
        pmf.Incr(p, prob)

    pmf.Normalize()
    return pmf


def NormalPdf(x):
    """Computes the PDF of x in the standard normal distribution."""
    return math.exp(-x**2/2) / math.sqrt(2 * math.pi)


def SubtractExpected(hist):
    """Subtracts from each cell in a grid the expected number in that cell.

    Args:
        hist: Hist object that maps sequences to frequencies.

    Returns:
        new Hist object
    """
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
    #myplot.Pmf(probs, show=True)

    pairs = SummarizeSequences(sequences)
    
    likelihood = LikelihoodSequences(pairs, sigma=0.1)
    print likelihood
    return

    # get the sex ratio (boys per hundred girls)
    genders = Pmf.Hist()
    for t in sequences:
        for x in t:
            genders.Incr(x)

    boys, girls = genders.Freq('B'), genders.Freq('G')
    print boys, girls, 100.0 * boys / girls
    print boys, girls, 1.0 * boys / (boys+girls)

    pairs = SummarizeSequences(sequences)
    PrintGrid(pairs)

    diff = SubtractExpected(pairs)
    PrintGrid(diff)

    diag = DiagDiff(pairs)
    print pairs.Total(), diag

    return
    diags, cdf = SimulateDiags(1000)


if __name__ == '__main__':
    main()
