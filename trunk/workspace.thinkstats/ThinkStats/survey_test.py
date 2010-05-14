"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import survey

class Test(unittest.TestCase):

    def testMean(self):
        resp = survey.Respondents()
        resp.ReadRecords()
        self.assertEquals(len(resp.records), 7643)

    def testPregnancies(self):
        preg = survey.Pregnancies()
        preg.ReadRecords()
        self.assertEquals(len(preg.records), 13593)

if __name__ == "__main__":
    unittest.main()
