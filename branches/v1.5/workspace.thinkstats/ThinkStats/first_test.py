"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import first

import Pmf

class Test(unittest.TestCase):

    def testMakeTables(self):
        table, firsts, others = first.MakeTables()
        self.assertEquals(len(table.records), 13593)
        self.assertEquals(len(firsts.records), 4413)
        self.assertEquals(len(others.records), 4735)

def MakeHist(table, attr):
    t = [getattr(record, attr) for record in table.records]
    hist = Pmf.MakeHistFromList(t, name=attr)
    return hist

if __name__ == "__main__":
    unittest.main()
