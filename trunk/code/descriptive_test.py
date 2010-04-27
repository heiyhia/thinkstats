"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import matplotlib.pyplot as plt
import Pmf
import descriptive
import plot

class Test(unittest.TestCase):

    def testPlot(self):
        t = [1, 2, 2, 3, 5]
        hist = Pmf.MakeHist(t)
        plot.Hist(hist)

    def testDescriptive(self):
        constructor = descriptive.Pregnancies
        table = constructor()
        table.ReadRecords()

        print
        print 'Number of pregnancies', len(table.records)
        self.assertEquals(len(table.records), 13593)

        firsts, others = table.PartitionRecords(constructor)

        print 'Number of first babies', len(firsts.records)
        print 'Number of others', len(others.records)

        self.assertEquals(len(firsts.records), 4413)
        self.assertEquals(len(others.records), 4735)

        firsts.Process()
        others.Process()

        firsts.MakePmf(name='first babies')
        others.MakePmf(name='others')

        axis = [20, 50, 0, 2700]
        plot.Hists([firsts.hist, others.hist], 'nsfg_hist', axis=axis)

        axis = [20, 50, 0, 0.6]
        plot.Hists([firsts.pmf, others.pmf], 'nsfg_pmf', axis=axis)

        axis = [20, 50, 0, 1.0]
        plot.Cdfs([others.cdf, firsts.cdf], 'nsfg_cdf', axis=axis)

if __name__ == "__main__":
    unittest.main()
