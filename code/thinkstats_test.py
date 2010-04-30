import unittest
import thinkstats

class Test(unittest.TestCase):

    def testMean(self):
        t = [1, 1, 1, 3, 3, 591]
        mu = thinkstats.Mean(t)
        self.assertEquals(mu, 100)

    def testVar(self):
        t = [1, 1, 1, 3, 3, 591]
        mu = thinkstats.Mean(t)
        var1 = thinkstats.Var(t)
        var2 = thinkstats.Var(t, mu)
        
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
