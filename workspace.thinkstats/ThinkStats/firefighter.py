"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

Answer to a question here:

http://www.reddit.com/r/statistics/comments/w2t2f/theres_no_raskstatistics_would_you_help_me/


Exercise A: My score 19.2, low 7.7, high 29.6, average 18.9050
Exercise B: My score 18.3, low 8.4, high 26.9, average 18.4325
Exercise C: My score 20.1, low 1.1, high 21.0, average 9.8163
"""

import random
import bisect

from correlation import CorrelatedGenerator


data = [
    [19.2, 7.7, 29.6, 18.9050],
    [18.3, 8.4, 26.9, 18.4325],
    [20.1, 1.1, 21.0, 9.8163]
]


def standard_dev(low, high):
    """Estimates the standard deviation from the range.

    Based on the first 84-rankit, which is about 2.4.
    """
    return (high - low) / 4.8

def z_score(score, mean, sd):
    return (score - mean) / sd

stats = []
for score, low, high, mean in data:
    print score, low, high, mean
    sd = standard_dev(low, high)
    z = z_score(score, mean, sd)
    print sd, z

    stats.append((mean, sd))

real = sum(score for score, _,_,_ in data)
print real

def generate(stats, rho=0.4):
    """Generates triplets of exam scores with the given stats and correlation.

    stats: list of (mean, std) pairs
    rho: corellation coefficient
    """
    it = CorrelatedGenerator(0.4)
    xs = [it.next() for i in range(3)]
    return [x * sd + mean for (mean, sd), x in zip(stats, xs)]

def simulate(stats):
    """Simulate the other 83 test takers and see where the actual score falls.

    stats: list of (mean, std) pairs    
    """
    fake = [sum(generate(stats)) for i in range(83)]
    fake.sort()
    index = bisect.bisect(fake, real)
    return index


t = [simulate(stats) for i in range(100)]
t.sort()
print 84-t[5]
print 84-t[95]
