import math

import Cdf
import Pmf
import myplot

def NormalPdf(x):
    """Computes the PDF of x in the standard normal distribution."""
    return math.exp(-x**2/2) / math.sqrt(2 * math.pi)

def Frange(low, high, n):
    return [low + (high-low) * float(i)/(n-1) for i in range(n)]

def RenderPdf(mu, sigma, n=1001):
    xs = Frange(mu-4*sigma, mu+4*sigma, n)
    ys = [NormalPdf((x-mu) / sigma) for x in xs]
    return xs, ys

def MakeFigure():
    xs, ys = RenderPdf(1170, 179)
    pmf1 = Pmf.MakePmfFromDict(dict(zip(xs, ys)), name='blue')
    total = pmf1.Total()
    pmf1.Normalize(total / 0.9)
    print pmf1.Total()

    xs, ys = RenderPdf(995, 167)
    pmf2 = Pmf.MakePmfFromDict(dict(zip(xs, ys)), name='green')
    total = pmf2.Total()
    pmf2.Normalize(total / 0.1)
    print pmf2.Total()

    #myplot.Pmfs([pmf1, pmf2], show=True)

    cdf1 = Cdf.MakeCdfFromPmf(pmf1)
    cdf2 = Cdf.MakeCdfFromPmf(pmf2)
    myplot.Cdfs([cdf1, cdf2], show=True)

MakeFigure()
