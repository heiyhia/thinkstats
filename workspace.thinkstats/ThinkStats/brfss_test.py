"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import brfss

import Pmf

class Test(unittest.TestCase):

    def testRespondents(self):
        resp = brfss.Respondents()
        resp.ReadRecords(n=10000)
        self.assertEquals(len(resp.records), 10000)

        hist = MakeHist(resp, 'wtkg2')
        t = hist.Values()
        low, high = min(t), max(t)
        self.assertAlmostEquals(low, 22.73)
        self.assertEquals(hist.Freq('NA'), 343)

        hist = MakeHist(resp, 'weight2')
        t = hist.Values()
        low, high = min(t), max(t)
        self.assertAlmostEquals(low, 22.727272727)
        self.assertEquals(hist.Freq('NA'), 343)

        hist = MakeHist(resp, 'wtyrago')
        t = hist.Values()
        low, high = min(t), max(t)
        self.assertAlmostEquals(low, 27.27272727)
        self.assertEquals(hist.Freq('NA'), 616)

        hist = MakeHist(resp, 'htm3')
        t = hist.Values()
        low, high = min(t), max(t)
        self.assertAlmostEquals(low, 104)
        self.assertEquals(hist.Freq('NA'), 101)

        hist = MakeHist(resp, 'sex')
        self.assertEquals(hist.Freq(1), 3669)
        self.assertEquals(hist.Freq(2), 6331)


def MakeHist(table, attr):
    t = [getattr(record, attr) for record in table.records]
    hist = Pmf.MakeHistFromList(t, name=attr)
    return hist

if __name__ == "__main__":
    unittest.main()
