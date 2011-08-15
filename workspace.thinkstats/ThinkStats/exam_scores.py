"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import random
import sys

import thinkstats

def flip(p):
    """Returns 1 with probability p; 0 otherwise."""
    if random.random() < p:
        return 1
    else:
        return 0

def fake_exam(pc, num_questions):
    """Generates a exam score for a student with probability pc."""
    return sum(flip(pc) for i in range(num_questions))

def fake_diff(pc, num_questions):
    """Generates the difference in two exam scores for a student
    with probability pc."""
    return fake_exam(pc, num_questions) - fake_exam(pc, num_questions)

def fake_diffs(pcs, num_questions):
    """Generates differences in exam scores for students with given
    values of pc."""
    diffs = [fake_diff(pc, num_questions) for pc in pcs]
    return diffs

def p_value(delta, pcs, num_questions, iterations=10000):
    """Computes the probability of seeing a mean difference in exam
    scores >= delta for students with given values of pc."""
    count = 0.0
    for i in range(iterations):
        if thinkstats.Mean(fake_diffs(pcs, num_questions)) >= delta:
            count += 1
    return count / iterations

def generate_pcs(mean, std, num_students):
    """Generates random values for pc."""
    return [random.normalvariate(mean, std) for i in range(num_students)]

def main(script):
    random.seed(17)

    mean = 0.7
    std = 0.15
    num_students = 30
    pcs = generate_pcs(mean, std, num_students)

    delta = 2
    num_questions = 50
    print p_value(delta, pcs, num_questions)


if __name__ == '__main__':
    main(*sys.argv)
