"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import glob
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


class Standard(object):
    time_table = """
    10-34	3:10	3:40
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

    def __init__(self):
        self.times = []
        lines = self.time_table.split('\n')
        for line in lines:
            t = line.split()
            if not t:
                continue
            ages, male, female = t
            ages = [int(x) for x in ages.split('-')]
            male += ':59'
            female += ':59'
            self.times.append((ages, male, female))

    def LookupAgeGroup(self, gender, age):
        for ages, male, female in self.times:
            low, high = ages
            if low <= age <= high:
                return '%s%2.2d%d' % (gender, low, high)
        return None

    def LookupQualifyingTime(self, gender, age):
        for ages, male, female in self.times:
            if ages[0] <= age <= ages[1]:
                if gender == 'M':
                    return male
                else:
                    return female
        return None

standard = Standard()


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


def CleanLine(line, half=False):
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

    if age == 0:
        return None

    gender = line[40]

    if gender not in ['M', 'F']:
        return None

    qual_time = standard.LookupQualifyingTime(gender, age)
    qual_time = ConvertTimeToMinutes(qual_time)

    net = ConvertTimeToMinutes(net)

    return gender, age, gun, net, qual_time, pace


def ReadResults(filename='result-10-all.php', half=False):
    """Read results from coolrunning and return a list of tuples."""
    results = []
    for line in open(filename):
        t = CleanLine(line, half)
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


def PartitionResults(res):
    groups = {}
    for t in res:
        gender, age, gun, net, qual_time, pace = t
        group = standard.LookupAgeGroup(gender, age)
        groups.setdefault(group, []).append(t)
    return groups


def FractionQualified(res, qual_time):
    count = 0
    for gender, age, gun, net, qual_time, pace in res:
        if net <= qual_time:
            count += 1
    return count, len(res)


def ComputeFractions(groups):
    for group, res in sorted(groups.iteritems()):
        qual_time = res[0][4]
        qual, total = FractionQualified(res, qual_time)
        print group, qual_time, qual, total


def main():
    filenames = glob.glob('result-*-all.php')

    diff_list = []
    all_res = []
    for filename in sorted(filenames):
        year = filename.split('-')[1]
        res, diffs = Process(filename)
        diff_list.append((year, diffs))
        all_res.extend(res)

    groups = PartitionResults(all_res)
    ComputeFractions(groups)
    return

    PlotDiffs(diff_list)

def Process(filename):
    res = ReadResults(filename=filename)
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
