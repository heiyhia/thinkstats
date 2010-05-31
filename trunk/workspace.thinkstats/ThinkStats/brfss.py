"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import cPickle
import continuous
import Cdf
import math
import matplotlib
import matplotlib.pyplot as pyplot
import myplot
import random
import sys
import survey
import thinkstats

class Respondent(survey.Respondent): 
    """Represents a respondent."""
    

class Respondents(survey.Table):
    """Represents the respondent table."""

    def ReadRecords(self, filename='CDBRFS08.ASC.gz', n=None):
        self.ReadFile(filename, self.GetFields(), Respondent, n)

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
        
    def SummarizeHeight(self):
        d = {1:[], 2:[]}
        [d[r.sex].append(r.htm3) for r in self.records
                                     if r.htm3 != 999]
        
        for key, t in d.iteritems():
            mu, var = thinkstats.TrimmedMeanVar(t)
            sigma = math.sqrt(var)
            cv = sigma / mu
            print key, mu, var, sigma, cv
        
    def SummarizeWeight(self):
        
        data = [(r.weight2, r.wtyrago) for r in self.records
                    if r.weight2 <= 999 and r.wtyrago <= 999]
        
        changes = [(curr - prev) for curr, prev in data]
            
        print 'Mean change', thinkstats.Mean(changes)
        
    def ScatterWeight(self):
        
        weights = []
        changes = []
        for r in self.records:
            if r.weight2 > 999 or r.wtyrago > 999:
                continue
            
            if r.weight2 == r.wtyrago:
                continue
            
            jitter = 3
            fudge = random.uniform(-jitter, jitter)
            change = (r.wtyrago - r.weight2) + jitter
            
            weights.append(r.weight2)
            changes.append(change)
            
        print 'Mean change', thinkstats.Mean(changes)
        
        
        # pyplot.scatter(weights, changes, s=1)
        pyplot.hexbin(weights, changes, cmap=matplotlib.cm.gray_r)
        myplot.Plot('brfss_scatter',
                  title = 'Weight change vs. weight',
                  xlabel = 'Current weight (pounds)',
                  ylabel = 'Weight change (pounds)',
                  axis = [50, 350, -50, 50],
                  legend=False,
                  show=True,
                  )


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
    
    resp = Respondents()
    resp.ReadRecords(n=100000)
    #resp.MakeFigures()
    #resp.SummarizeHeight()
    resp.SummarizeWeight()
    
if __name__ == '__main__':
    main(*sys.argv)
