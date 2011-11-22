"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

Credits:

erf is a modified version of a function contributed to the public
domain by John D. Cook
http://www.johndcook.com/python_erf.html
"""

import math
from scipy.special import erf, erfinv

root2 = math.sqrt(2.0)


def NormalCdf(x, mu=0, sigma=1):
    """Evaluates the CDF of the normal distribution.
    
    Args:
        x: float

        mu: mean parameter
        
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    xx = (x - mu) / sigma
    y = (erf(xx / root2) + 1) / 2
    return y


def NormalCdfInverse(p, mu=0, sigma=1):
    """Evaluates the inverse CDF of the normal distribution.
    
    Args:
        p: float

        mu: mean parameter
        
        sigma: standard deviation parameter
                
    Returns:
        float
    """
    x = root2 * erfinv(2*p - 1)
    return mu + x * sigma

