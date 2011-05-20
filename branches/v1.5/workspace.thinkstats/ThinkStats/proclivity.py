"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import random
import matplotlib.pyplot as pyplot
import sys

import Cdf
import myplot
import Pmf

import all_survey
import survey
import thinkstats


def MakeUniformSuite(low, high, steps):
    """Makes a PMF that represents a suite of hypotheses with equal prob.
    
    Args:
        low: low end of range
        high: high end of range
        steps: number of values

    Returns:
        Pmf object
    """
    hypos = [low + (high-low) * i / (steps-1.0) for i in range(steps)]
    pmf = Pmf.MakePmfFromList(hypos)
    return pmf


def GetSequences(table, oknums=[1]):
    """Reads a pregnancy table and returns a list of baby gender sequences.

    Args:
        table: Pregnancies object
        oknums: list of nbrnaliv numbers that are considered ok
                (in other words, should we include the first baby
                from multiple births?)

    Returns:
        list of sequences
    """
    sequences = {}
    d={1:'B', 2:'G'}

    for record in table.records:
        num, gender = record.nbrnaliv, record.babysex

        if num in oknums and gender in [1,2]:
            sequences.setdefault(record.caseid, []).append(d[gender])

    return sequences.values()


def RandomSequence(pb, num=2, d={True:'B', False:'G'}):
    """Makes a random sequence with length num."""
    t = [d[random.random()<pb] for i in range(num)]
    return t


def RandomSequences(pb, num=2, seqs=1621):
    """Makes a list of random sequences."""
    t = [RandomSequence(pb, num) for i in range(seqs)]
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


def LogProbSequence(sequence, pb):
    """Computes the probability of a sequence under the simple model.

    Assumes that the probability of a boy for all births is pb.
    """
    d = dict(B=math.log(pb), G=math.log(1-pb))
    t = [d[x] for x in sequence]
    return sum(t)
    

def ProbSequence(sequence, pb):
    """Computes the probability of a sequence under the simple model.

    Assumes that the probability of a boy for all births is pb.
    """
    prob = math.exp(LogProbSequence(sequence, pb))
    return prob
    

def ProbSequenceSigma(sequence, probs):
    """Computes the probability of a sequence under the complex model.

    Assumes that the probability of a boy is described by the given
    distribution.
    """
    total = 0.0
    for mp, prob in probs.Items():
        total += prob * ProbSequence(sequence, mp)
    return total
    

def LogLikelihood(hist, sigma, pb):
    """Computes the likelihood of a list of sequences given sigma.

    Args:
        hist: Histogram that maps sequences to their frequencies
        sigma: float std of proclivities

    Returns:
        float log likelihood
    """
    probs = ComputeProbs(pb, sigma=sigma)

    total = 0.0
    for sequence, freq in hist.Items():
        prob = ProbSequenceSigma(sequence, probs)
        total += math.log(prob) * freq
    return total


def Likelihood(hist, sigma, pb):
    """Computes the likelihood of a list of sequences given sigma.

    Args:
        hist: Histogram that maps sequences to their frequencies
        sigma: float std of proclivities

    Returns:
        float likelihood
    """
    logprob = LogLikelihood(hist, sigma, pb)
    prob = math.exp(logprob)
    return prob


def ComputeProbs(pb, sigma=0.5, n=101, bound=2.5):
    """Make a Pmf of probabilities with the given parameter.

    The probabilities are the logistic tranform of the normal values.

    Args:
        p: float mean probability.
        sigma: float standard deviation of the normal values.
        n: number of values in the Pmf

    Returns:
        Pmf object
    """
    bias = pb - 0.5
    low = -bound
    high = bound
    xs = [low + i * (high-low) / (n-1) for i in range(n)]

    pmf = Pmf.Pmf(name='proclivities')
    for x in xs:
        p = 1 / (1 + math.exp(-x * sigma)) + bias
        prob = NormalPdf(x)
        pmf.Incr(p, prob)

    pmf.Normalize()
    return pmf


def NormalPdf(x):
    """Computes the PDF of x in the standard normal distribution."""
    return math.exp(-x**2/2) / math.sqrt(2 * math.pi)


def SubtractExpected(hist, pb):
    """Subtracts from each cell in a grid the expected number in that cell.

    Args:
        hist: Hist object that maps sequences to frequencies.

    Returns:
        new Hist object
    """
    total = hist.Total()
    diff = hist.Copy('diff')

    for val, freq in hist.Items():
        expected = total * ProbSequence(val, pb)
        diff.Incr(val, -expected)

    return diff


def SummarizeSequences(sequences, num=2):
    """Make a histogram of num-baby sequences.

    Args:
        sequences: list of string baby sequences
        num: sequence length

    Returns:
        Hist object that maps string sequences to frequencies
    """
    pairs = Pmf.Hist()
    for t in sequences:

        # use the first num babies
        t = t[0:num]
        if len(t) != num:
            continue
        s = ''.join(t)
        pairs.Incr(s)
    return pairs


def GenderRatioPairs(pairs):
    """

    Returns:
        tuple of (num boys, num girls, probability of a boy)
    """
    genders = Pmf.Hist()
    for pair, freq in pairs.Items():
        for x in pair:
            genders.Incr(x, freq)

    boys, girls = genders.Freq('B'), genders.Freq('G')
    boysper = boys, girls, 100.0 * boys / girls
    pb = boys, girls, 1.0 * boys / (boys+girls)

    return pb


def GenderRatioSequences(sequences):
    """

    Returns:
        tuple of (num boys, num girls, probability of a boy)
    """
    genders = Pmf.Hist()
    for t in sequences:
        for x in t:
            genders.Incr(x)

    boys, girls = genders.Freq('B'), genders.Freq('G')
    boysper = boys, girls, 100.0 * boys / girls
    pb = boys, girls, 1.0 * boys / (boys+girls)

    return pb


def SimulateDiags(num=1000):
    diags = []
    for i in range(num):
        fake = RandomSequences()
        pairs = SummarizeSequences(fake)
        diag = DiagDiff(pairs)
        diags.append(diag)

    # print thinkstats.Mean(diags)

    cdf = Cdf.MakeCdfFromList(diags, 'diagonal bias')
    myplot.Cdf(cdf, show=True)

    return diags, cdf


def LogUpdate(prior, evidence, pb):
    """Updates a prior distribution based on new evidence.

    Prior and result are expressed in log likelihoods.

    Args:
        prior: Pmf object
        evidence: whatever kind of object Likelihood expects
    """
    posterior = prior.Copy()
    
    for hypo in prior.Values():
        loglikelihood = LogLikelihood(evidence, hypo, pb)
        posterior.Incr(hypo, loglikelihood)

    ShiftPmf(posterior)
    return posterior


def ShiftPmf(pmf):
    items = pmf.Items()
    vals, probs = zip(*items)
    shift = -max(probs)

    for val, prob in items:
        pmf.Incr(val, shift)


def LogPmf(pmf):
    """Applies a log transform to the probs in a PMF."""
    return _MapPmf(pmf, math.log)


def ExpPmf(pmf):
    """Applies an exp transform to the probs in a PMF."""
    return _MapPmf(pmf, math.exp)


def _MapPmf(pmf, func):
    new = Pmf.Pmf()
    for val, prob in pmf.Items():
        new.Incr(val, func(prob))
    return new


def EstimateSigma(pairs, pb):

    # make the prior
    prior = MakeUniformSuite(0.0, 1.0, 1001)
    prior = LogPmf(prior)

    # make the posterior
    posterior = LogUpdate(prior, pairs, pb)
    posterior.name='sigma'
    posterior = ExpPmf(posterior)
    posterior.Normalize()

    return posterior


def PlotPosteriorSigma(posterior):

    ci = CredibleInterval(posterior, 90)
    print 'CI:', ci

    pyplot.clf()
    PlotCredibleInterval(posterior, ci)

    myplot.Pmf(posterior,
               root='sigma',
               clf=False,
               title='Posterior PMF',
               xlabel='sigma',
               ylabel='probability',
               show=True)


def CredibleInterval(pmf, percentage):
    """Computes a credible interval for a given distribution.

    If percentage=90, computes the 90% CI.

    Args:
        pmf: Pmf object representing a posterior distribution
        percentage: float between 0 and 100

    Returns:
        sequence of two floats, low and high
    """
    cdf = Cdf.MakeCdfFromDict(pmf.GetDict())
    prob = (1 - percentage/100.0) / 2
    interval = [cdf.Value(p) for p in [prob, 1-prob]]
    return interval


def PlotCredibleInterval(pmf, ci):
    for x in ci:
        p = pmf.Prob(x)
        xs = [x, x]
        ps = [0, p]
        pyplot.plot(xs, ps, linewidth=2, color='red')


def Summarize(sequences):
    diag = DiagDiff(pairs)
    print pairs.Total(), diag

    return
    diags, cdf = SimulateDiags(1000)


def MakeTable(constructor=survey.Pregnancies):
    """Makes a table and reads a pregnancy data file.

    Args:
        constructor: class object for a subtype of Pregnancies

    Returns:
        new Pregnancies object (of the given subtype)
    """
    table = constructor()
    table.ReadRecords()
    return table


def main(script, num=1):

    num = int(num)

    constructors = [
        survey.Pregnancies,
        all_survey.Pregnancies1995,
        all_survey.Pregnancies1988,
        all_survey.Pregnancies1982,
        all_survey.Pregnancies1976,
        ]

    constructors = constructors[:num]

    all_sequences = []
    for constructor in constructors:

        # table contains one record for each pregnancy
        table = MakeTable(constructor)

        # sequences contains one item for each respondent with at least
        # one child with known gender
        # (multiple births are excluded by default; alternatively, we
        # could only keep track of the first baby in a multiple birth)
        sequences = GetSequences(table)
        ProcessSequences(sequences, constructor.__name__)

        all_sequences.extend(sequences)

    ProcessSequences(all_sequences, 'all data')


def ProcessSequences(sequences, name):
    
    print
    print name

    boys, girls, pb = GenderRatioSequences(sequences)
    print boys, girls, pb

    # pairs contains one entry for each respondent with at least two
    # children with known gender
    pairs = SummarizeSequences(sequences, num=2)

    boys, girls, pb = GenderRatioPairs(pairs)
    print boys, girls, pb

    print 'pairs'
    PrintGrid(pairs)

    diff = SubtractExpected(pairs, pb)
    print 'diff'
    PrintGrid(diff)
    
    #posterior = EstimateSigma(pairs, pb)
    #PlotPosteriorSigma(posterior)

    trips = SummarizeSequences(sequences, num=3)
    ProcessTrips(trips)


def ProcessTrips(trips):
    for prefix in ['BB', 'BG', 'GB', 'GG']:
        boys = trips.Freq(prefix + 'B')
        girls = trips.Freq(prefix + 'G')
        pb = Frac(boys, girls)
        print prefix, boys+girls, pb


def Frac(x, y):
    return 1.0 * x / (x + y)


if __name__ == '__main__':
    main(*sys.argv)
