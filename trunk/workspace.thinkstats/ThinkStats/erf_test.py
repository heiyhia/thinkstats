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


if __name__ == "__main__":
    unittest.main()
