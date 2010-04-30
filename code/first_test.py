"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import first

class Test(unittest.TestCase):

    def testSummarize(self):
        table, firsts, others = first.MakeTables()
        first.ProcessTables(table, firsts, others)

        self.assertEquals(len(table.records), 13593)
        self.assertEquals(len(firsts.records), 4413)
        self.assertEquals(len(others.records), 4735)

        self.assertAlmostEquals(firsts.mu, 38.6009517335)
        self.assertAlmostEquals(others.mu, 38.5229144667)


if __name__ == "__main__":
    unittest.main()
