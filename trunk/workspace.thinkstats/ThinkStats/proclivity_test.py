"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import proclivity

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
        self.assertAlmostEquals(p, 0.27177189186277162)

        p = proclivity.ProbSequenceSigma('BG', probs)
        self.assertAlmostEquals(p, 0.23622810813722883)

    def testLikelihoodSequences(self):
        sequences = ['BB']
        likelihood = proclivity.LikelihoodSequences(sequences, sigma=0.5)
        self.assertAlmostEquals(likelihood, 0.27177189186277162)

        sequences = ['BB', 'BG', 'GB', 'GG']
        likelihood = proclivity.LikelihoodSequences(sequences, sigma=0.5)
        self.assertAlmostEquals(likelihood, 0.0038790064091085911)

if __name__ == "__main__":
    unittest.main()
