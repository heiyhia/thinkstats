"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import sys

import Pmf

def RemainingLifetime(pmf, age):
    """Takes a Pmf of lifetimes and an age, and returns a new Pmf 
    that represents the distribution of remaining lifetimes.
    """
    rem = Pmf.Pmf()
    for val, prob in pmf.Items():
        if (val > age):
            rem.Incr(val-age, prob)
    rem.Normalize()
    return rem


def main(script):
    pmf = Pmf.MakePmfFromList(range(10))
    rem = RemainingLifetime(pmf, 5)

    for val, prob in sorted(rem.Items()):
        print val, prob


if __name__ == '__main__':
    main(*sys.argv)
