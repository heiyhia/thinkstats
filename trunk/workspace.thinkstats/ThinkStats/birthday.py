"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import math
import random
import matplotlib.pyplot as pyplot
import Cdf
import myplot
import Pmf


cdfs = []
allbday = []

for i in range(10):
    n = 30
    t = [random.randrange(365) for i in range(n)]
    t.sort()

    pmf = Pmf.Pmf()
    for i in range(len(t)-1):
        x = t[i+1] - t[i]
        pmf.Incr(x)
        allbday.append(x)

    cdf = Cdf.MakeCdfFromPmf(pmf)
    cdfs.append(cdf)

cdf = Cdf.MakeCdfFromList(allbday, 'all')
cdfs.append(cdf)
myplot.Cdfs(cdfs, root='birthday', transform='exponential')
