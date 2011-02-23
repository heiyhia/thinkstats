"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import glob
import sys

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

    def __init__(self, offset=0):
        self.offset = offset
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
                    return male, self.offset
                else:
                    return female, self.offset
        return None, self.offset


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

    qual_time, offset = standard.LookupQualifyingTime(gender, age)
    qual_time = ConvertTimeToMinutes(qual_time) + offset

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


def ReadAllChicago():
    filenames = ['Chicago2008.csv',
                 'Chicago2009.csv',
                 'Chicago2010.csv',
                 ]

    all_res = []
    for filename in filenames:
        res = ReadChicago(filename)
        all_res.extend(res)

    return all_res


def ReadChicago(filename='Chicago2010.csv'):
    """Read results from coolrunning and return a list of tuples."""
    results = []
    for line in open(filename):
        try:
            age, gender, net = line.split()
        except ValueError:
            continue

        gender = gender.upper()

        age = int(age)
        if age < 10 or age > 34:
            continue

        gun = None
        pace = None

        qual_time, offset = standard.LookupQualifyingTime(gender, age)
        qual_time = ConvertTimeToMinutes(qual_time) + offset

        net = ConvertTimeToMinutes(net)

        t = gender, age, gun, net, qual_time, pace
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


def dappend(d, key, val):
    d.setdefault(key, []).append(val)


def ComputeFractions(groups):
    """
    Args:
        groups: map from string group to list of result tuples
    """
    # map from gender to list of qualifier counts
    qualifiers = {}

    for group, res in sorted(groups.iteritems()):

        # get stats from the first person in the list
        sample = res[0]
        gender = sample[0]
        qual_time = sample[4]

        qual, total = FractionQualified(res, qual_time)
        dappend(qualifiers, gender, qual)
        #print group, qual_time, qual, total

    return qualifiers


def GenderRatio(qualifiers):
    men = sum(qualifiers['M'])
    women = sum(qualifiers['F'])
    ratio = float(women) / (men + women)
    return men, women, ratio


def ComputeField(men, women, field=9602):
    factor = float(field) / (men + women)
    return men*factor, women*factor


def PartitionGenders(res):
    d = {}
    for t in res:
        gender = t[0]
        dappend(d, gender, t)
    return d


def ComputeDiffs(results, low, high):
    diffs = []
    for gender, age, gun, net, qual_time, pace in results:
        diff = net - qual_time
        if low <= diff <= high:
            diffs.append(diff)
    return diffs


def PlotDiffs(groups, low, high, root):
    diff_list = []
    for gender, res in groups.iteritems():
        diffs = ComputeDiffs(res, low=low, high=high)
        diff_list.append((gender, diffs))
        print 'PlotDiffs', gender, len(diffs)

    cdfs = []
    for name, diffs in diff_list:
        cdf = Cdf.MakeCdfFromList(diffs, name=name)
        cdfs.append(cdf)

    options = [dict(linewidth=2) for cdf in cdfs]

    myplot.Cdfs(cdfs, 
                xlabel='time - qualifying time (min)',
                ylabel='P(difference < x)',
                plot_options=options,
                root=root)


def PlotPmf(results):
    speeds = GetSpeeds(results)
    pmf = Pmf.MakePmfFromList(speeds, 'speeds')
    myplot.Pmf(pmf, 
               title='PMF of running speed',
               xlabel='speed (mph)',
               ylabel='probability',
               show=True)


def SummarizeChange(s, old, new):
    print 'Change', s, old, new, new-old, float(new-old) / old


def SummarizeImpact(men, new_men, women, new_women):
    old_field = ComputeField(men, women)
    new_field = ComputeField(new_men, new_women)
    #print 'old field', old_field
    #print 'new field', new_field
    print 'Impact (change in number of men)', new_field[0] - old_field[0]


def ReadCapeCod():
    all_res = []
    filenames = glob.glob('result-*-all.php')

    for filename in sorted(filenames):
        year = filename.split('-')[1]

        res = ReadResults(filename=filename)
        all_res.extend(res)
    return all_res


def GetContenders(res, gender, cutoff, spread):
    contenders = []
    count = 0
    low, high = cutoff-spread, cutoff+spread
    for t in res:
        g, age, gun, net, qual_time, pace = t
        if g == gender and low <= net <= high:
            contenders.append(t)
            if net <= cutoff:
                count += 1

    fraction = float(count) / len(contenders)
    return contenders, count, fraction


def FindFairStandard(res):
    spread = ConvertTimeToMinutes('00:30:00')

    male_cutoff = ConvertTimeToMinutes('3:05:00')
    males, count, fraction = GetContenders(res, 'M', male_cutoff, spread)
    print len(males), count, fraction

    offsets = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
    for offset in offsets:
        cutoff = ConvertTimeToMinutes('3:35:00') + offset
        females, count, fraction = GetContenders(res, 'F', cutoff, spread)
        print offset, len(females), count, fraction


def RunAnalysis(offset=0):
    global standard
    standard = Standard(offset=offset)

    res = ReadAllChicago()
    groups = PartitionResults(res)

    if offset == 0:
        for group, res in groups.iteritems():
            print 'Participants', group, len(res)
        PlotDiffs(groups, low=-190, high=30, root='bq_cdf1')
        PlotDiffs(groups, low=-190, high=-30, root='bq_cdf2')

    qualifiers = ComputeFractions(groups)
    men, women, ratio = GenderRatio(qualifiers)
    return men, women, ratio


def main():
    res = []
    offsets = [0, -1, -6, -11, -21]
    for offset in offsets:
        t = RunAnalysis(offset=offset)
        print 'Qualifiers', t
        res.append(t)

    men, women, ratio = res[0]
    for offset, t in zip(offsets, res):
        print
        print 'offset', offset

        new_men, new_women, new_ratio = t
        SummarizeChange('M', men, new_men)
        SummarizeChange('F', women, new_women)
        SummarizeImpact(men, new_men, women, new_women)

    


if __name__ == '__main__':
    main()
