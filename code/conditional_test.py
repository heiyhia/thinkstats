"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import Pmf
import Cdf
import descriptive
import conditional

class Test(unittest.TestCase):

    def testConditional(self):
        constructor = descriptive.Pregnancies
        all_table = constructor()
        all_table.ReadRecords()

        firsts, others = all_table.PartitionRecords(constructor)

        all_table.Process()
        firsts.Process()
        others.Process()

        all_table.MakeDist(name='all')        
        firsts.MakeDist(name='first babies')
        others.MakeDist(name='others')

        weeks = range(30, 45)
        d = {}

        for week in weeks:
            print
            for table in [all_table, firsts, others]:

                cond = ConditionOnWeeks(table.pmf, week)
                prob = cond.Prob(week)
                print week, prob, table.cdf.name
                d[week, table.cdf.name] = prob

        print
        for week in weeks:
            ratio = d[week, 'first babies'] / d[week, 'others']
            print week, ratio

def ConditionOnWeeks(pmf, week=39):
    def filter_func(x):
        return x < week

    cond = conditional.ConditionPmf(pmf, filter_func)
    return cond

if __name__ == "__main__":
    unittest.main()
