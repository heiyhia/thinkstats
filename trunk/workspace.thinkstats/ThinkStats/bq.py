"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import urllib

import Cdf
import myplot
import Pmf

"""
Sample line.

Place Div/Tot  Name                 No.   Age S Club                 City       
        St Time    Net Time  Pace  Qual 
===== ======== ==================== ===== === = ==================== ===========
========== ======= ========= ===== ==== 
    1   1/186  PATRICK MOULTON       2260  28 M BOSTON ATHLETIC ASSO PROVIDENCE 
RI         2:24:41 2:24:38.2  5:32 *

"""

time_table = """
0-34	3:10	3:40
35-39	3:15	3:45
40-44	3:20	3:50
45-49	3:30	4:00
50-54	3:35	4:05
55-59	3:45	4:15
60-64	4:00	4:30
65-69	4:15	4:45
70-74	4:30	5:00
75-79	4:45	5:15
80-99	5:00	5:30
"""

def QualifyingTimes(time_table=time_table):
    times = []
    for line in time_table.split('\n'):
        if not line:
            continue
        ages, male, female = line.split()
        ages = [int(x) for x in ages.split('-')]
        male += ':59'
        female += ':59'
        times.append((ages, male, female))
    return times


def LookupQualifyingTime(gender, age, times):
    for ages, male, female in times:
        if ages[0] <= age <= ages[1]:
            if gender == 'M':
                return male
            else:
                return female
    print gender, age
    return None


def ConvertPaceToSpeed(pace):
    """Converts pace in MM:SS format to MPH."""
    m, s = [int(x) for x in pace.split(':')]
    secs = m*60 + s
    mph  = 1.0 / secs * 60 * 60 
    return mph


def ConvertTimeToMinutes(time):
    """Converts pace in HH:MM:SS format to minutes."""
    t = time.split('.')
    if len(t) == 2:
        time, fraction = time.split('.')
    
    h, m, s = [int(x) for x in time.split(':')]
    mins = h * 60 + m + s / 60.0
    return mins


def CleanLine(line, times, half=False):
    """Converts a line from coolrunning results to a tuple of values."""
    t = line.split()
    if len(t) < 6:
        return None
    
    net, pace, gun = t[-3:]

    for time in [gun, net, pace]:
        if ':' not in time:
            return None

    age = line[37:39]
    age = int(age)
    gender = line[40]

    if gender not in ['M', 'F']:
        return None

    qual_time = LookupQualifyingTime(gender, age, times)
    qual_time = ConvertTimeToMinutes(qual_time)

    net = ConvertTimeToMinutes(net)

    return gender, age, gun, net, qual_time, pace


def ReadResults(times, filename='result-10-all.php', half=False):
    """Read results from coolrunning and return a list of tuples."""
    results = []
    for line in open(filename):
        t = CleanLine(line, times, half)
        if t:
            results.append(t)
    return results


def GetSpeeds(results, column=-1):
    """Extract the pace column and return a list of speeds in MPH."""
    speeds = []
    for t in results:
        pace = t[column]
        speed = ConvertPaceToSpeed(pace)
        speeds.append(speed)
    return speeds


def main():
    times = QualifyingTimes()
    
    filenames = [('08', 'result-08-all.php'),
                 ('09', 'result-09-all.php'),
                 ('10', 'result-10-all.php'),
                 ]

    diff_list = []
    for year, filename in filenames:
        res, diffs = Process(times, filename)
        diff_list.append((year, diffs))

    PlotDiffs(diff_list)

def Process(times, filename):
    res = ReadResults(times, filename=filename)
    diffs = ComputeDiffs(res)

    print len(res), len(diffs)

    return res, diffs


def ComputeDiffs(results):
    diffs = []
    for gender, age, gun, net, qual_time, pace in results:
        diff = net - qual_time
        if abs(diff) < 60:
            diffs.append(diff)
    return diffs


def PlotDiffs(diff_list):
    cdfs = []
    for year, diffs in diff_list:
        cdf = Cdf.MakeCdfFromList(diffs, name=year)
        cdfs.append(cdf)

    options = [dict(linewidth=2) for cdf in cdfs]

    myplot.Cdfs(cdfs, 
                xlabel='time - qualifying time (min)',
                ylabel='CDF',
                plot_options=options,
                root='bq_cdf')


def PlotPmf(results):
    speeds = GetSpeeds(results)
    pmf = Pmf.MakePmfFromList(speeds, 'speeds')
    myplot.Pmf(pmf, 
               title='PMF of running speed',
               xlabel='speed (mph)',
               ylabel='probability',
               show=True)


if __name__ == '__main__':
    main()
