"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import descriptive

class Test(unittest.TestCase):

    def testVar(self):
        pool, firsts, others = descriptive.MakeTables()
        descriptive.Process(pool, 'Live births')
        descriptive.Process(firsts, 'First babies')
        descriptive.Process(others, 'Others')

        self.assertAlmostEquals(firsts.var, 7.79294720207)
        self.assertAlmostEquals(others.var, 6.84123839078)

        self.assertAlmostEquals(pool.mu, 38.5605596852)
        self.assertAlmostEquals(pool.var, 7.3018637882)

if __name__ == "__main__":
    unittest.main()
