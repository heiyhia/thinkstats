"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import erf

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


if __name__ == "__main__":
    unittest.main()
