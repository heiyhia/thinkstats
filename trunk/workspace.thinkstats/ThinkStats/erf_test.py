"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest

import erf
import myplot

class Test(unittest.TestCase):

    def testErf(self):
        y = erf.erf(0.0)
        self.assertAlmostEquals(y, 0.0)

        y = erf.erf(0.1)
        self.assertAlmostEquals(y, 0.11246296562219549)

        y = erf.erf(-0.1)
        self.assertAlmostEquals(y, -0.11246296562219549)

    def testErfinv(self):
        x = erf.erfinv(0.0)
        self.assertAlmostEquals(x, 0.0)

        x = erf.erfinv(0.11246296562219549)
        self.assertAlmostEquals(x, 0.1)

        x = erf.erfinv(-0.11246296562219549)
        self.assertAlmostEquals(x, -0.1)

    def testMakeNormalCdf(self):
        cdf = erf.MakeNormalCdf(digits=2)
        self.assertAlmostEquals(cdf.Prob(-3.0), 0.0013498980316301035)
        self.assertAlmostEquals(cdf.Prob(0.0), 0.5)
        self.assertAlmostEquals(cdf.Prob(3.0), 0.9986501019683699)
        self.assertAlmostEquals(cdf.Prob(5.0), 0.99996832875816688)

    def testMakeNormalPmf(self):
        pmf = erf.MakeNormalPmf(digits=1)

        self.assertAlmostEquals(pmf.Prob(-3.0), 0.00038229481841173252)
        self.assertAlmostEquals(pmf.Prob(0.0), 0.039827837277028821)
        self.assertAlmostEquals(pmf.Prob(3.0), 0.00051591526875394145)
       
    def testFixedPointNormalPmf(self):
        pmf = erf.FixedPointNormalPmf(spread=3, digits=1)
        
        self.assertAlmostEquals(pmf.NormalProb(-1.01), 0.9888250336816421)
        self.assertAlmostEquals(pmf.NormalProb(0.0), 1.629618)
        self.assertAlmostEquals(pmf.NormalProb( 1.01), 0.9888250336816421)


if __name__ == "__main__":
    unittest.main()
