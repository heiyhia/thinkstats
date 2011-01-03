"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import unittest
import proclivity

import Pmf

class Test(unittest.TestCase):

    def testProbSequence(self):
        p = proclivity.ProbSequence('BB', pb=0.5)
        self.assertAlmostEquals(p, 0.25)
        p = proclivity.ProbSequence('BG', pb=0.508)
        self.assertAlmostEquals(p, 0.249936)

    def testProbSequenceSigma(self):
        probs = proclivity.ComputeProbs(pb=0.508, sigma=0.0, n=101, bound=2.5)
        p = proclivity.ProbSequenceSigma('BG', probs)
        self.assertAlmostEquals(p, 0.249936)

        probs = proclivity.ComputeProbs(pb=0.508, sigma=0.5, n=101, bound=2.5)
        p = proclivity.ProbSequenceSigma('BB', probs)
        self.assertAlmostEquals(p, 0.2710908)

        p = proclivity.ProbSequenceSigma('BG', probs)
        self.assertAlmostEquals(p, 0.23690915)

    def testLikelihood(self):
        sequences = ['BB']
        hist = Pmf.MakeHistFromList(sequences)
        likelihood = proclivity.Likelihood(hist, sigma=0.5, pb=0.508)
        self.assertAlmostEquals(likelihood, 0.2710908)

        sequences = ['BB', 'BG', 'GB', 'GG']
        hist = Pmf.MakeHistFromList(sequences)
        likelihood = proclivity.Likelihood(hist, sigma=0.5, pb=0.508)
        self.assertAlmostEquals(likelihood, 0.003881266018)

        loglikelihood = proclivity.LogLikelihood(hist, sigma=0.5, pb=0.508)
        self.assertAlmostEquals(loglikelihood, -5.5515938852)
        likelihood = math.exp(loglikelihood)
        self.assertAlmostEquals(likelihood, 0.003881266018)

    def testLogExpPmf(self):
        pmf = Pmf.MakePmfFromList([1, 3])
        new = proclivity.LogPmf(pmf)
        self.assertAlmostEquals(new.Prob(1), -0.69314718056)

        old = proclivity.ExpPmf(new)
        self.assertAlmostEquals(old.Prob(1), 0.5)

    def testMakeTable(self):
        table = proclivity.MakeTable()
        self.assertEquals(len(table.records), 13593)

        sequences = proclivity.GetSequences(table, oknums=[1])
        self.assertEquals(len(sequences), 4387)
        
        pairs = proclivity.SummarizeSequences(sequences)
        self.assertEquals(pairs.Freq('GG'), 709)
        self.assertEquals(pairs.Freq('BB'), 740)

if __name__ == "__main__":
    unittest.main()
