"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import random
import thinkbayes

class Test(unittest.TestCase):

    def testPmfProbLess(self):
        pmf = thinkbayes.MakePmfFromList([1,2,2,3,3,3,4,4,4,4])
        p = pmf.ProbLess(3)
        self.assertAlmostEquals(p, 0.3)

        p = pmf < 4
        self.assertAlmostEquals(p, 0.6)

        p = pmf < pmf
        self.assertAlmostEquals(p, 0.35)

        p = pmf > pmf
        self.assertAlmostEquals(p, 0.35)

        p = pmf <= pmf
        self.assertAlmostEquals(p, 0.65)

        p = pmf >= pmf
        self.assertAlmostEquals(p, 0.65)

        p = pmf == pmf
        self.assertAlmostEquals(p, 0.3)

        p = pmf != pmf
        self.assertAlmostEquals(p, 0.7)


    def testGaussianCdf(self):
        x = 1.0
        p = thinkbayes.GaussianCdf(x, mu=0, sigma=1)
        self.assertAlmostEquals(p, 0.841344746069)

        x = thinkbayes.GaussianCdfInverse(p, mu=0, sigma=1)  
        self.assertAlmostEquals(x, 1.0)


if __name__ == "__main__":
    unittest.main()
