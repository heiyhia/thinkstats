"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import unittest
import survey
import first_baby
import descriptive

class Test(unittest.TestCase):

    def testFirstBaby(self):
        constructor = first_baby.Pregnancies
        table, firsts, others = self.runTest(constructor)

    def testDescriptive(self):
        constructor = descriptive.Pregnancies
        table, firsts, others = self.runTest(constructor)

        var1, var2  = firsts.var, others.var

        print
        print 'Variance'
        print 'First babies', var1 
        print 'Others', var2

        self.assertAlmostEquals(var1, 7.79294720207)
        self.assertAlmostEquals(var2, 6.84123839078)

        diff_mu = firsts.mu - others.mu

        print 'Difference in mean', diff_mu

        pool = first_baby.PoolRecords(constructor, firsts, others)
        pool.Process()
        sigma = math.sqrt(pool.var)

        print 'Pooled mean', pool.mu
        print 'Pooled variance', pool.var
        print 'Pooled sigma', sigma

        self.assertAlmostEquals(pool.mu, 38.5605596852)
        self.assertAlmostEquals(pool.var, 7.3018637882)

    def runTest(self, constructor):
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

        table.Process()
        firsts.Process()
        others.Process()

        mu1, mu2 = firsts.mu, others.mu

        self.assertAlmostEquals(mu1, 38.6009517335)
        self.assertAlmostEquals(mu2, 38.5229144667)

        print 'Mean gestation in weeks:' 
        print 'First babies', mu1 
        print 'Others', mu2
        
        print 'Difference in days', (mu1 - mu2) * 7.0

        return table, firsts, others

    def testVar(self):
        t = [1, 1, 1, 3, 3, 591]
        mu = descriptive.Mean(t)
        var1 = descriptive.Var(t)
        var2 = descriptive.Var(t, mu)
        
        print
        print 'Pumpkins'
        print 'mean', mu 
        print 'var', var1
        print 'var', var2

        self.assertAlmostEquals(mu, 100.0)
        self.assertAlmostEquals(var1, 48217.0)
        self.assertAlmostEquals(var2, 48217.0)

if __name__ == "__main__":
    unittest.main()
