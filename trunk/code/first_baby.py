"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import sys
import survey

def main(name):
    """Computes summary statistics from survey data."""
    preg = survey.Pregnancies()
    print 'Number of pregnancies', len(preg.records)

    firsts = []
    others = []

    for p in preg.records:
        if p.outcome != 1:
            continue

        if p.birthord == 1:
            firsts.append(p)
        else:
            others.append(p)

    print 'Number of first babies', len(firsts)
    print 'Number of others', len(others)

    first_lens = [p.prglength for p in firsts]
    other_lens = [p.prglength for p in others]

    mu1 = Mean(first_lens)
    mu2 = Mean(other_lens)

    print 'Mean gestation in weeks:' 
    print 'First babies', mu1 
    print 'others', mu2

    print 'Difference in days', (mu1 - mu2) * 7.0


def Mean(t):
    """Computes the mean of a sequence of numbers.

    Args:
        t: sequence of numbers

    Returns:
        float
    """
    return float(sum(t)) / len(t)

if __name__ == '__main__':
    main(*sys.argv)
