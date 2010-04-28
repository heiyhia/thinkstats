"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import Pmf
import Cdf

class Test(unittest.TestCase):

    # TODO: rewrite testHist and testCdf
    def testHist(self):
        t = [1, 2, 2, 3, 5]
        hist = Pmf.MakeHist(t)
        plot.Hist(hist)

    def ReadTables(self):
        constructor = descriptive.Pregnancies
        table = constructor()
        table.ReadRecords()

        firsts, others = table.PartitionRecords(constructor)

        table.Process()
        firsts.Process()
        others.Process()

        table.MakeDist(name='all')        
        firsts.MakeDist(name='first babies')
        others.MakeDist(name='others')

        return firsts, others

    def testConditional(self):
        firsts, others = self.ReadTables()

        self.assertAlmostEquals(firsts.cdf.Prob(40), 0.844550192613)
        self.assertAlmostEquals(others.cdf.Prob(40), 0.906230200634)

        axis = [23, 46, 0, 2700]
        plot.Hists([firsts.hist, others.hist], 
                   'nsfg_hist', 
                   title='Histogram',
                   xlabel='weeks',
                   ylabel='frequency',
                   axis=axis)

        axis = [23, 46, 0, 0.6]
        plot.Hists([firsts.pmf, others.pmf],
                   'nsfg_pmf',
                   title='PMF',
                   xlabel='weeks',
                   ylabel='probability',
                   axis=axis)

        axis = [23, 46, 0, 1.0]
        plot.Cdfs([others.cdf, firsts.cdf],
                  'nsfg_cdf',
                  title='CDF',
                  xlabel='weeks',
                  ylabel='probability',
                  styles=[':', '-'],
                  axis=axis)

    def testTrimmedMean(self):
        firsts, others = self.ReadTables()

        print firsts.mu, others.mu
        print firsts.trim, others.trim

        weeks = range(35, 46)
        ratios = []
        for week in weeks:
            p1 = firsts.pmf.Prob(week)
            p2 = others.pmf.Prob(week)
            diff = 100 * (p1 - p2)
            ratios.append(diff)

        plt.clf()
        plt.bar(weeks, ratios, align='center')
        plot.Plot('nsfg_ratios',
                  title='Difference in PMFs',
                  xlabel='weeks',
                  ylabel='$100 ({PMF}_{first} - {PMF}_{other})$',
                  legend=False
                  )
        
def Log2(x, denom=math.log(2)):
    return math.log(x) / denom

if __name__ == "__main__":
    unittest.main()
