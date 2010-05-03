"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import Pmf

class Test(unittest.TestCase):

    def testHist(self):
        t = [1, 2, 2, 3, 5]
        hist = Pmf.MakeHist(t)

        self.assertEquals(hist.Freq(1), 1)
        self.assertEquals(hist.Freq(2), 2)
        self.assertEquals(hist.Freq(3), 1)
        self.assertEquals(hist.Freq(4), 0)
        self.assertEquals(hist.Freq(5), 1)
        
        pmf = Pmf.MakePmfFromHist(hist)
        self.checkPmf(pmf)

    def checkPmf(self, pmf):
        self.assertAlmostEquals(pmf.Prob(1), 0.2)
        self.assertAlmostEquals(pmf.Prob(2), 0.4)
        self.assertAlmostEquals(pmf.Prob(3), 0.2)
        self.assertAlmostEquals(pmf.Prob(4), 0.0)
        self.assertAlmostEquals(pmf.Prob(5), 0.2)

    def testMakePmf(self):
        t = [1, 2, 2, 3, 5]
        pmf = Pmf.MakePmfFromList(t)
        self.checkPmf(pmf)


if __name__ == "__main__":
    unittest.main()
