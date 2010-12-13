"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import cPickle
import math
import sys
import survey
import thinkstats

class Respondent(survey.Respondent): 
    """Represents a respondent."""
    

class Respondents(survey.Table):
    """Represents the respondent table."""

    def ReadRecords(self, data_dir='', filename='CDBRFS08.ASC.gz', n=None):
        self.ReadFile(data_dir, filename, self.GetFields(), Respondent, n)

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.
        
        BRFSS codebook 
        http://www.cdc.gov/brfss/technical_infodata/surveydata/2008.htm

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.
        """
        return [
            ('weight2', 119, 122, int),
            ('wtyrago', 127, 130, int),
            ('wtkg2', 1254, 1258, int),
            ('htm3', 1251, 1253, int),
            ('sex', 143, 143, int),
            ]

    def SummarizeHeight(self):
        d = {1:[], 2:[]}
        [d[r.sex].append(r.htm3) for r in self.records
                                     if r.htm3 != 999]
        
        for key, t in d.iteritems():
            mu, var = thinkstats.TrimmedMeanVar(t)
            sigma = math.sqrt(var)
            cv = sigma / mu
            print key, mu, var, sigma, cv
        
    def SummarizeWeightChange(self):
        
        data = [(r.weight2, r.wtyrago) for r in self.records
                    if r.weight2 <= 999 and r.wtyrago <= 999]
        
        changes = [(curr - prev) for curr, prev in data]
            
        print 'Mean change', thinkstats.Mean(changes)
        
    
def main(name, data_dir=''):
    resp = Respondents()
    resp.ReadRecords(data_dir)
    resp.SummarizeHeight()
    resp.SummarizeWeightChange()
    
if __name__ == '__main__':
    main(*sys.argv)
