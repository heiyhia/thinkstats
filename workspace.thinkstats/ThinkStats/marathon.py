"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import urllib

import myplot
import Pmf

results = 'Oct17_BaySta_set1.shtml'

"""
Sample line.

Place Div/Tot  Name                 No.   Age S Club                 City       
        St Time    Net Time  Pace  Qual 
===== ======== ==================== ===== === = ==================== ===========
========== ======= ========= ===== ==== 
    1   1/186  PATRICK MOULTON       2260  28 M BOSTON ATHLETIC ASSO PROVIDENCE 
RI         2:24:41 2:24:38.2  5:32 *

"""

def ConvertPaceToSpeed(pace):
    """Converts pace in HH:MM:SS format to MPH."""
    m, s = [int(x) for x in pace.split(':')]
    secs = m*60 + s
    mph  = 1.0 / secs * 60 * 60 
    return mph


def CleanLine(line):
    """Converts a line from coolrunning results to a tuple of values."""
    t = line.split()
    if len(t) < 6:
        return None
    
    place, divtot = t[:2]

    if not '/' in divtot:
        return None

    if t[-1] == '*':
        qual = True
        t = t[:-1]
    else:
        qual = False

    gun, net, pace = t[-3:]

    for time in [gun, net, pace]:
        if ':' not in time:
            return None

    age = line[42:45]
    age = int(age)
    gender = line[46]

    if gender not in ['M', 'F']:
        return None

    print place, age, gender
    return place, age, gun, net, pace


def ReadResults(filename=results):
    """Read results from coolrunning and return a list of tuples."""
    results = []
    for line in open(filename):
        t = CleanLine(line)
        if t:
            results.append(t)
    return results


def GetSpeeds(results, column=4):
    """Extract the pace column and return a list of speeds in MPH."""
    speeds = []
    for t in results:
        pace = t[column]
        speed = ConvertPaceToSpeed(pace)
        speeds.append(speed)
    return speeds


def main():
    results = ReadResults()
    print len(results)

    speeds = GetSpeeds(results)
    pmf = Pmf.MakePmfFromList(speeds, 'speeds')
    myplot.Pmf(pmf, 
               title='PMF of running speed',
               xlabel='speed (mph)',
               ylabel='probability',
               show=True)


if __name__ == '__main__':
    main()
