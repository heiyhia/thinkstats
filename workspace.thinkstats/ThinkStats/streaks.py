"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import random

"""
To prepare the input file:

more */random.txt > all_random.txt
python streaks.py


"""

def ReadFile(filename='all_random.txt'):
    """Reads a file containing the output from more."""
    fp = open(filename)
    for line in fp:
        line = line.strip()

        if len(line) == 0 or line.startswith(':'):
            continue
        if line[0] in '01':
            ProcessLine(name, line)
        else:
            name = line.split('/')[0]

def ProcessRandom(n=13):
    """Generate n truly random lines and run the tests on them."""
    for i in range(n):
        line = [random.choice('01') for i in range(100)]
        ProcessLine('random.choice', line)

def ProcessLine(name, line):
    """Run the tests on a string of 0s and 1s."""
    switches = Switches(line)
    run1, run2 = Runs(line)
    print '%s, %d, %d, %d' % (name, switches, run1, run2)

def Switches(line):
    """Count the number of 0-1 or 1-0 transitions."""
    switches = 0
    for x, y in zip(line, line[1:]):
        if x != y:
            switches += 1
    return switches

def Runs(line):
    """Count the lengths of the runs and return the length of the longest."""
    runs = []
    run = 1
    for x, y in zip(line, line[1:]):
        if x != y:
            runs.append(run)
            run = 1
        else:
            run += 1

    runs.append(run)
    runs.sort()
    return runs[-1], runs[-4]

print 'Student-generated sequences'
ReadFile()

print ''
print 'Truly random sequences'
ProcessRandom()

