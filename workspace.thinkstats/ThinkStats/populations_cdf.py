"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import Cdf
import math
import myplot
import populations
import random
import rankit
import thinkstats
import matplotlib.pyplot as pyplot


def MakeFigure():
    pops = populations.Process()
    print len(pops)
    
    cdf = Cdf.MakeCdfFromList(pops, 'populations')

    myplot.Cdfs([cdf], 'populations',
              title='City/Town Populations',
              xlabel='population',
              ylabel='CDF',
              legend=False)

    myplot.Cdfs([cdf], 'populations_logx',
              title='City/Town Populations',
              xlabel='population',
              ylabel='CDF',
              xscale='log',
              legend=False)

    myplot.Cdfs([cdf], 'populations_loglog',
              complement=True,
              title='City/Town Populations',
              xlabel='population',
              ylabel='Complementary CDF',
              yscale='log',
              xscale='log',
              legend=False)
    
    t = [math.log(x) for x in pops]
    t.sort()
    rankit.MakeNormalPlot(t, 'populations_rankit')

def main():
    MakeFigure()
    
if __name__ == "__main__":
    main()
