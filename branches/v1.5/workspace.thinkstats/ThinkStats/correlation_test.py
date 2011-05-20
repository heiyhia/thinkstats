"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import correlation

import thinkstats


class Test(unittest.TestCase):

    def testCov(self):
        xs = [1, 2, 3]
        ys = [3, 4, 5]
        cov = correlation.Cov(xs, ys)
        self.assertAlmostEquals(cov, 0.666666666)

        var = thinkstats.Var(xs)
        cov = correlation.Cov(xs, xs)
        self.assertAlmostEquals(var, cov)

    def testCorr(self):
        xs = [1, 2, 3]
        ys = [3, 4, 5]
        cor = correlation.Corr(xs, ys)
        self.assertAlmostEquals(cor, 1.0)

        xs = [1, 2, 100]
        ys = [3, 4, 5]
        cor = correlation.Corr(xs, ys)
        self.assertAlmostEquals(cor, 0.8703878312633373)

        cor = correlation.Corr(xs, xs)
        self.assertAlmostEquals(cor, 1.0)

    def testMapToRanks(self):
        t = correlation.MapToRanks([7, 1, 2, 10, 5])
        expected = [4, 1, 2, 5, 3]
        for got, exp in zip(t, expected):
            self.assertEquals(got, exp)

    def testSpearmanCorr(self):
        # the outlier doesn't affect rank correlation
        xs = [1, 2, 100]
        ys = [3, 4, 5]
        cor = correlation.SpearmanCorr(xs, ys)
        self.assertAlmostEquals(cor, 1.0)

    def testLeastSquares(self):
        xs = [1, 2, 3]
        ys = [3, 6, 8]
        inter, slope = correlation.LeastSquares(xs, ys)
        self.assertAlmostEquals(inter, 0.66666666)
        self.assertAlmostEquals(slope, 2.5)

        res = correlation.Residuals(xs, ys, inter, slope)
        for got, exp in zip(res, [-0.166666666, 0.33333333, -0.16666666666]):
            self.assertAlmostEquals(got, exp)

        R2 = correlation.CoefDetermination(ys, res)
        self.assertAlmostEquals(R2, 0.986842105263)

if __name__ == "__main__":
    unittest.main()
