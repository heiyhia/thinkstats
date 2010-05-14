"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import sys
import survey
import continuous
import Cdf
import math
import thinkstats
import matplotlib.pyplot as pyplot
import myplot
import cPickle


class Respondent(survey.Respondent): 
    """Represents a respondent."""
    

class Respondents(survey.Table):
    """Represents the respondent table."""

    def ReadRecords(self, filename='CDBRFS08.ASC.gz'):
        self.ReadFile(filename, self.GetFields(), Respondent)

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.
        """
        return [
            ('wtkg2', 1254, 1258, int),
            ]

    def MakeFigure(self, weights, root,
                   xmax=175, 
                   xlabel='adult weight (kg)',
                   axis=None):
        cdf = Cdf.MakeCdfFromList(weights)
                
        pyplot.clf()
        
        t = weights[:]
        t.sort()
        mu, var = thinkstats.TrimmedMeanVar(t)
        print 'n, Mean, Var', len(weights), mu, var
        
        sigma = math.sqrt(var)
        print 'Sigma', sigma

        xs, ps = continuous.RenderNormalCdf(mu, sigma, xmax)
        pyplot.plot(xs, ps, label='model', linewidth=3, color='0.7')
    
        xs, ps = cdf.Render()
        pyplot.plot(xs, ps, label='data', color='0.0')
     
        myplot.Plot(root,
                  title = 'Adult weight',
                  xlabel = xlabel,
                  ylabel = 'CDF',
                  axis=axis or [0, xmax, 0, 1])
    
    def MakeFigures(self):
        weights = [record.wtkg2/100.0 for record in self.records
                   if record.wtkg2 != 99999]
        self.MakeFigure(weights, 'brfss_weight_model')
        
        log_weight = [math.log(weight) for weight in weights]
        xmax = math.log(175.0)
        axis = [3.5, 5.2, 0, 1]
        self.MakeFigure(log_weight, 'brfss_weight_log',
                        xmax=xmax,
                        xlabel='adult weight (log kg)',
                        axis=axis)
    


def WritePickle(filename='brfss.pkl'):
    resp = Respondents()
    resp.ReadRecords()
    print 'Number of respondents', len(resp.records)
    print 'Writing', filename
    fp = open(filename, 'wb')
    cPickle.dump(resp, fp)
    fp.close()

    
def ReadPickle(filename='brfss.pkl'):
    fp = open(filename)
    obj = cPickle.load(fp)
    fp.close()
    return obj


def main(name):
    #WritePickle()
    
    resp = ReadPickle()
    resp.MakeFigures()
    
if __name__ == '__main__':
    main(*sys.argv)
