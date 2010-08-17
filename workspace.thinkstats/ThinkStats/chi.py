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


def main():
    Test()


def Test():
    pool, firsts, others = descriptive.MakeTables()

    funcs = [risk.ProbEarly, risk.ProbOnTime, risk.ProbLate]
    
    print 'observed'
    observed = ComputeRows(firsts, others, funcs, probs=None)
    print observed

    tables = [firsts, others]
    probs = [func(pool.pmf) for func in funcs]
    print 'expected'
    expected = ComputeRows(firsts, others, funcs, probs=probs)
    print expected
    
    print 'chi-squared'
    threshold = ChiSquared(expected, observed)
    print threshold

    print 'simulated'
    num_trials = 10000
    chi2s = []
    count = 0
    for _ in range(num_trials):
        simulated = ComputeRows(firsts, others, funcs, probs=probs, 
                            row_func=SimulateRow)
        chi2 = ChiSquared(expected, simulated)
        chi2s.append(chi2)
        if chi2 >= threshold:
            count += 1
            
    print 'max chi2'
    print max(chi2s)
    
    print 'p-value'
    print 1.0 * count / num_trials
    

def ComputeRow(n, probs):
    """Multiplies out a row of a table.
    
    Args:
      n: sum of the elements in the row
      prob: sequence of float probabilities
    """
    row = [n * prob for prob in probs]
    return row


def SimulateRow(n, probs):
    """Generates a random row of a table.
    
    Chooses all but the last element at random, then chooses the
    last element to make the sums work.
    
    Args:
      n: sum of the elements in the row
      prob: sequence of float probabilities
    """
    row = [Binomial(n, prob) for prob in probs]
    row[-1] += n - sum(row)
    return row


def Binomial(n, prob):
    """Returns a random sample from a binomial distribution.
    
    Args:
      n: int number of trials
      prob: float probability
      
    Returns:
      int: number of successes
    """
    t = [1 for _ in range(n) if random.random() < prob]
    return sum(t)


def ComputeRows(firsts, others, funcs, probs=None, row_func=ComputeRow):
    """Computes a table suitable for use with chi-squared stats.
    
    There are three uses of this function:
    
    1) To compute observed values, use probs=None and row_func=ComputeRow
    
    2) To compute expected values, provide probs from the pooled data,
        and row_func=ComputeRow
        
    3) To generate random values, provide pooled probs,
        and row_func=SimulateRow
        
    Returns:
      row of rows of float values
    """
    rows = []
    for table in [firsts, others]:
        n = len(table)
        row_probs = probs or [func(table.pmf) for func in funcs]
        row = row_func(n, row_probs)
        rows.append(row)
    
    return rows


def ChiSquared(expected, observed):
    """Compute the Chi-squared statistic for two tables.
    
    Args:
      expected: row of rows of values
      observed: row of rows of values
      
    Returns:
      float chi-squared statistic
    """
    it = zip(itertools.chain(*expected), 
             itertools.chain(*observed))
    t = [(obs - exp)**2 / exp for exp, obs in it]
    return sum(t)


if __name__ == "__main__":
    main()
