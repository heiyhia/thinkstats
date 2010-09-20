"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import brfss
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

class Respondent(brfss.Respondent): 
    """Represents a respondent."""
    

class Respondents(brfss.Respondents):
    """Represents the respondent table."""

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



def main(name):
    resp = Respondents()
    resp.ReadRecords(n=100000)
    resp.MakeFigures()

    
if __name__ == '__main__':
    main(*sys.argv)
