"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import descriptive
import itertools
import Pmf
import random
import risk
import thinkstats


def SimulateRolls(sides, num_rolls):
    """Generates a Hist of simulated die rolls.
    
    Args:
      sides: number of sides on the die
      num_rolls: number of times to rolls

    Returns:
      Hist object
    """
    hist = Pmf.Hist()
    for i in range(num_rolls):
        roll = random.randint(1, sides)
        hist.Incr(roll)
    return hist


def MakeExpected(sides, num_rolls):
    """Generates a Hist of expected die rolls.
    
    Args:
      sides: number of sides on the die
      num_rolls: number of times to rolls

    Returns:
      Hist object
    """
    hist = Pmf.Hist()
    for i in range(sides):
        exp = num_rolls / float(sides)
        hist.Incr(i, exp)
    return hist
    

def ChiSquared(expected, observed):
    """Compute the Chi-squared statistic for two tables.
    
    Args:
      expected: Hist of expected values
      observed: Hist of observed values
      
    Returns:
      float chi-squared statistic
    """
    # TODO: check whether there are any observed values with 0 expected
    total = 0.0
    for x, exp in expected.Items():
        obs = observed.Freq(x)
        total += (obs - exp)**2 / exp
    return total


def Test(expected, observed, sides=6, num_trials=1000):
    """Run a simulation to estimate the p-value of the observed values.

    Args:
      expected: Hist of expected values
      observed: Hist of observed values
      sides: number of sides on the die
      num_trials: how many simulations to run

    Returns:
      float p-value
    """

    # compute the chi-squared stat
    threshold = ChiSquared(expected, observed)
    print 'chi-squared', threshold

    print 'simulated %d trials' % num_trials
    chi2s = []
    count = 0
    num_rolls = observed.Total()
    for _ in range(num_trials):
        simulated = SimulateRolls(sides, num_rolls)
        chi2 = ChiSquared(expected, simulated)
        chi2s.append(chi2)
        if chi2 >= threshold:
            count += 1
            
    print 'max chi2', max(chi2s)
    
    pvalue = 1.0 * count / num_trials
    print 'p-value', pvalue

    return pvalue


def main():
    # make a Hist of observed values
    d = {1:8, 2:9, 3:19, 4:6, 5:8, 6:10}
    observed = Pmf.MakeHistFromDict(d)
    print observed.d

    # make a Hist of expected values
    sides = 6
    total = observed.Total()
    expected = MakeExpected(sides, total)
    print expected.d

    # compute the chi-squared statistic
    print ChiSquared(expected, observed)

    # Test(expected, observed, num_trials=1000)


if __name__ == "__main__":
    main()
