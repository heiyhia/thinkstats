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
        p = proclivity.ProbSequence('BB', p=0.5)
        self.assertAlmostEquals(p, 0.25)
        p = proclivity.ProbSequence('BG')
        self.assertAlmostEquals(p, 0.249936)

    def testProbSequenceSigma(self):
        probs = proclivity.ComputeProbs(sigma=0.0)
        p = proclivity.ProbSequenceSigma('BG', probs)
        self.assertAlmostEquals(p, 0.249936)

        probs = proclivity.ComputeProbs(sigma=0.5)
        p = proclivity.ProbSequenceSigma('BB', probs)
        self.assertAlmostEquals(p, 0.2710908)

        p = proclivity.ProbSequenceSigma('BG', probs)
        self.assertAlmostEquals(p, 0.23690915)

    def testLikelihood(self):
        sequences = ['BB']
        hist = Pmf.MakeHistFromList(sequences)
        likelihood = proclivity.Likelihood(hist, sigma=0.5)
        self.assertAlmostEquals(likelihood, 0.2710908)

        sequences = ['BB', 'BG', 'GB', 'GG']
        hist = Pmf.MakeHistFromList(sequences)
        likelihood = proclivity.Likelihood(hist, sigma=0.5)
        self.assertAlmostEquals(likelihood, 0.003881266018)

        loglikelihood = proclivity.LogLikelihood(hist, sigma=0.5)
        self.assertAlmostEquals(loglikelihood, -5.5515938852)
        likelihood = math.exp(loglikelihood)
        self.assertAlmostEquals(likelihood, 0.003881266018)

    def testLogExpPmf(self):
        pmf = Pmf.MakePmfFromList([1, 3])
        new = proclivity.LogPmf(pmf)
        self.assertAlmostEquals(new.Prob(1), -0.69314718056)

        old = proclivity.ExpPmf(new)
        self.assertAlmostEquals(old.Prob(1), 0.5)

if __name__ == "__main__":
    unittest.main()
