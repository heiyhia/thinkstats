"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import rarefaction

import Cdf
import Pmf
import myplot

class Test(unittest.TestCase):

    def testMeta(self):
        sample = 'a'
        prior = Pmf.MakePmfFromList(range(1, 4))
        meta = rarefaction.MetaHypo(prior, sample)

        hypos = meta.GetHypos()
        for hypo, prob in hypos.Items():
            self.assertAlmostEquals(prob, 1.0/3.0)
            
            if hypo.k != 3:
                continue

            taxa = hypo.GetTaxa()
            for taxon, dist in taxa.iteritems():
                self.assertAlmostEquals(dist.Mean(), 1.0/3.0)

        meta.Update(sample)

        hypos = meta.GetHypos()

        for hypo, prob in hypos.Items():
            self.assertAlmostEquals(prob, 1.0/3.0)
            
            if hypo.k != 2:
                continue

            taxa = hypo.GetTaxa()
            for taxon, dist in taxa.iteritems():
                p = dist.Mean()
                self.assertAlmostEquals(p, 2.0/3.0) if taxon=='a' else 0
                self.assertAlmostEquals(p, 1.0/3.0) if taxon=='other' else 0


    def testGeneratePrevalence(self):
        sample = 'ab'
        prior = Pmf.MakePmfFromList(range(1, 4))
        meta = rarefaction.MetaHypo(prior, sample)

        for taxon in sample:
            meta.Update(taxon)

        k = 3
        taxon = 'other'
        iters = 1000
        plot = False

        hypos = meta.GetHypos()
        for hypo, prob in hypos.Items():
            self.assertAlmostEquals(prob, 0.4) if hypo.k==2 else 0
            self.assertAlmostEquals(prob, 0.6) if hypo.k==3 else 0

            if hypo.k != k or not plot:
                continue

            # check the distribution of generated prevalences
            ps = [hypo.GeneratePrevalence().Prob(taxon) for i in xrange(iters)]
            cdf = Cdf.MakeCdfFromList(ps)

            # compare to what the distribution is supposed to be
            dist = hypo.Get(taxon)
            ps2 = [dist.Random() for i in xrange(iters)]
            cdf2 = Cdf.MakeCdfFromList(ps2)

            myplot.Cdfs([cdf, cdf2], show=True)


if __name__ == "__main__":
    unittest.main()
