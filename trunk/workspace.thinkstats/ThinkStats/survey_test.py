"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import survey

import Pmf

class Test(unittest.TestCase):

    def testMean(self):
        resp = survey.Respondents()
        resp.ReadRecords()
        self.assertEquals(len(resp.records), 7643)

    def testPregnancies(self):
        preg = survey.Pregnancies()
        preg.ReadRecords()
        self.assertEquals(len(preg.records), 13593)

        hist = MakeHist(preg, 'nbrnaliv')
        self.assertEquals(hist.Freq(1), 8981)
        
        hist = MakeHist(preg, 'babysex')
        self.assertEquals(hist.Freq(1), 4641)
        self.assertEquals(hist.Freq(2), 4500)
        
        hist = MakeHist(preg, 'outcome')
        self.assertEquals(hist.Freq(1), 9148)
        
        hist = MakeHist(preg, 'birthord')
        self.assertEquals(hist.Freq(1), 4413)

        hist = MakeHist(preg, 'birthwgt_lb')
        self.assertEquals(hist.Freq(6), 2223)

        hist = MakeHist(preg, 'birthwgt_oz')
        self.assertEquals(hist.Freq(6), 709)

        hist = MakeHist(preg, 'agepreg')
        self.assertEquals(hist.Freq('NA'), 352)
        self.assertEquals(hist.Freq(25.0), 58)

        hist = MakeHist(preg, 'totalwgt_oz')
        self.assertEquals(hist.Freq('NA'), 4509)

        hist = MakeHist(preg, 'finalwgt')
        t = hist.Values()
        low, high = min(t), max(t)
        self.assertAlmostEquals(low, 118.656789706)
        self.assertAlmostEquals(high, 261879.9538641)
        

def MakeHist(table, attr):
    t = [getattr(record, attr) for record in table.records]
    hist = Pmf.MakeHistFromList(t, name=attr)
    return hist

if __name__ == "__main__":
    unittest.main()
