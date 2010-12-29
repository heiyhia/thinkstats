"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import random
import matplotlib.pyplot as pyplot

import Cdf
import myplot
import Pmf
import survey
import thinkstats


def MakeUniformSuite(low, high, steps):
    """Makes a PMF that represents a suite of hypotheses with equal p.
    
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


def LogProbSequence(sequence, p):
    """Computes the probability of a sequence under the simple model.

    Assumes that the probability of a boy for all births is p.
    """
    d = dict(B=math.log(p), G=math.log(1-p))
    t = [d[x] for x in sequence]
    return sum(t)
    

def ProbSequence(sequence, p=0.508):
    """Computes the probability of a sequence under the simple model.

    Assumes that the probability of a boy for all births is p.
    """
    prob = math.exp(LogProbSequence(sequence, p))
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
    

def LogLikelihood(hist, sigma):
    """Computes the likelihood of a list of sequences given sigma.

    Args:
        hist: Histogram that maps sequences to their frequencies
        sigma: float std of proclivities

    Returns:
        float log likelihood
    """
    probs = ComputeProbs(sigma=sigma)
    #print 'probs'
    #probs.Print()

    total = 0.0
    for sequence, freq in hist.Items():
        prob = ProbSequenceSigma(sequence, probs)
        total += math.log(prob) * freq
    return total


def Likelihood(hist, sigma):
    """Computes the likelihood of a list of sequences given sigma.

    Args:
        hist: Histogram that maps sequences to their frequencies
        sigma: float std of proclivities

    Returns:
        float likelihood
    """
    logprob = LogLikelihood(hist, sigma)
    prob = math.exp(logprob)
    return prob


def ComputeProbs(p=0.508, sigma=0.5, n=101, bound=2.5):
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

    # print thinkstats.Mean(diags)

    cdf = Cdf.MakeCdfFromList(diags, 'diagonal bias')
    myplot.Cdf(cdf, show=True)

    return diags, cdf


def LogUpdate(prior, evidence):
    """Updates a prior distribution based on new evidence.

    Prior and result are expressed in log likelihoods.

    Args:
        prior: Pmf object
        evidence: whatever kind of object Likelihood expects
    """
    posterior = prior.Copy()
    
    for hypo in prior.Values():
        loglikelihood = LogLikelihood(evidence, hypo)
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

def main():
    posterior = EstimateSigma()
    PlotPosteriorSigma(posterior)

def EstimateSigma():

    # read the first-two-baby sequences
    table, sequences = MakeTable()
    pairs = SummarizeSequences(sequences)
    
    # make the prior
    prior = MakeUniformSuite(0.0, 1.0, 1001)
    prior = LogPmf(prior)

    # make the posterior
    posterior = LogUpdate(prior, pairs)
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


def Summarize(sequences):
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
