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
    """Represents a qualifying standard"""

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
        """Initializes the standard.

        offset: how much to add to all qualifying times

        self.times is a list of (agerange, male qualtime, female qualtime)
        """
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
        """Returns the age group string for a given gender and age."""
        for ages, male, female in self.times:
            low, high = ages
            if low <= age <= high:
                return '%s%2.2d%d' % (gender, low, high)
        return None

    def LookupQualifyingTime(self, gender, age):
        """Returns (qualtime, offset) for a given gender and age."""
        for ages, male, female in self.times:
            if ages[0] <= age <= ages[1]:
                if gender == 'M':
                    return male, self.offset
                else:
                    return female, self.offset
        return None, self.offset


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
    """Reads results from coolrunning and returns a list of tuples."""
    results = []
    for line in open(filename):
        t = CleanLine(line, half)
        if t:
            results.append(t)
    return results


def ReadAllCapeCod():
    """Reads data from the Cape Cod Marathon."""
    all_res = []
    filenames = glob.glob('result-*-all.php')

    for filename in sorted(filenames):
        year = filename.split('-')[1]

        res = ReadResults(filename=filename)
        all_res.extend(res)
    return all_res


def ReadAllChicago():
    """Reads all data from Chicago.

    Returns:
        list of res
    """
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


def PartitionResults(res):
    """Returns map from group string to list of res tuples."""
    groups = {}
    for t in res:
        gender, age, gun, net, qual_time, pace = t
        group = standard.LookupAgeGroup(gender, age)
        groups.setdefault(group, []).append(t)
    return groups


def NumberQualified(res, qual_time):
    """Returns (number of qualifiers, number of res)."""
    count = 0
    for t in res:
        net = t[3]
        if net <= qual_time:
            count += 1
    return count, len(res)


def dappend(d, key, val):
    """Appends a value to the list of values in a map."""
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

        qual, total = NumberQualified(res, qual_time)
        dappend(qualifiers, gender, qual)
        #print group, qual_time, qual, total

    return qualifiers


def FindQualifiers(groups, qual_times):
    """Makes a histogram of the number of qualifiers in each group.

    Args:
        groups: map from string group to list of result tuples
        qual_times: map from group to qualifying time in minutes

    Returns:
        Histogram that maps from group to number of qualifiers
    """
    qualifiers = Pmf.Hist()

    for group, res in groups.iteritems():

        qual_time = qual_times[group]
        qual, total = NumberQualified(res, qual_time)
        qualifiers.Incr(group, qual)

    return qualifiers


def GenderRatio(qualifiers):
    men = sum(qualifiers['M'])
    women = sum(qualifiers['F'])
    ratio = float(women) / (men + women)
    return men, women, ratio


def ComputeField(men, women, field=9602):
    """Returns the (men, women) in a field with the given size.

    Args:
        # of men and women in the population from which the field is formed.
    """
    factor = float(field) / (men + women)
    return men*factor, women*factor


def PartitionGenders(res):
    d = {}
    for t in res:
        gender = t[0]
        dappend(d, gender, t)
    return d


def ComputeDiffs(res, low, high):
    """Returns the difference between times and qualifying times.

    Args:
        res: list of results
        low, high: range of diffs to include
    """
    diffs = []
    for gender, age, gun, net, qual_time, pace in res:
        diff = net - qual_time
        if low <= diff <= high:
            diffs.append(diff)
    return diffs


def PlotDiffs(groups, low, high, root):
    """Plots the CDF of diffs for each group.

    Args:
        low, high: range of diffs to include
        root: string filename root
    """
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


def SummarizeChange(s, old, new):
    """Prints old, new, diff, percent diff."""
    print 'Change', s, old, new, new-old, 100.0 * (new-old) / old


def GetContenders(res, gender, cutoff, spread):
    """Returns a list of contenders.

    Args:
        res: list of results
        gender: string M or F
        cutoff: qual time in minutes
        spread: minutes +- to collect results

    Returns:
        tuple of (list of results, number of results, fraction qualified)
    """
    contenders = []
    count = 0
    low, high = cutoff-spread, cutoff+spread
    for t in res:
        g, age, gun, net, qual_time, pace = t
        if g == gender and low <= net <= high:
            contenders.append(t)
            if net <= cutoff:
                count += 1

    return contenders, count, Fraction(count, len(contenders))
                                       

def FindFairStandard(res):
    """Finds the standard that qualifies the same fraction of contenders.
    """
    spread = ConvertTimeToMinutes('00:30:00')

    male_cutoff = ConvertTimeToMinutes('3:05:00')
    males, count, fraction = GetContenders(res, 'M', male_cutoff, spread)
    print len(males), count, fraction

    # try out a range of offsets for the female qualifying time
    offsets = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
    for offset in offsets:
        cutoff = ConvertTimeToMinutes('3:35:00') + offset
        females, count, fraction = GetContenders(res, 'F', cutoff, spread)
        print offset, len(females), count, fraction


def ReadGroups(offset=0):
    """Reads Chigaco results and partitions into groups."""
    global standard
    standard = Standard(offset=offset)

    res = ReadAllChicago()
    groups = PartitionResults(res)
    return groups


def RunAnalysis(offset=0):
    """Computes the field for a given gender gap.

    Args:
        offset: minutes added to the current standard.
    """
    groups = ReadGroups(offset)
    qualifiers = ComputeFractions(groups)
    men, women, ratio = GenderRatio(qualifiers)
    return men, women, ratio


def MakeGraphs():
    """Generate figures showing the distribution of diffs."""
    groups = ReadGroups()
    
    PlotDiffs(groups, low=-190, high=30, root='bq_cdf1')
    PlotDiffs(groups, low=-190, high=-30, root='bq_cdf2')


def EvaluateBurfoot():
    """Evaluate the impact of the Burfoot-proposed standard."""
    groups = ReadGroups()

    qual_times = dict(M1034=ConvertTimeToMinutes('3:10:00'),
                      F1034=ConvertTimeToMinutes('3:40:00'))
    before = FindQualifiers(groups, qual_times)

    print 'before'
    SummarizeQualifiers(before)
    print

    qual_times = dict(M1034=ConvertTimeToMinutes('3:12:00'),
                      F1034=ConvertTimeToMinutes('3:30:00'))
    after = FindQualifiers(groups, qual_times)

    print 'after'
    SummarizeQualifiers(after)
    print

    print 'difference'
    SummarizeDifference(before, after)
    print


def SummarizeQualifiers(qualifiers):
    """Prints a summary of qualifiers.

    qualifiers: histogram that maps groups to number of qualifiers
    """
    print qualifiers.Items()

    pmf = Pmf.MakePmfFromHist(qualifiers)
    print pmf.Items()


def SummarizeDifference(before, after):
    """Summarize the difference between two fields.

    before, after: histograms that map groups to number of qualifiers
    """
    for group in before.Values():
        diff = after.Freq(group) - before.Freq(group)
        print group, before.Freq(group), after.Freq(group), 
        print diff, 100.0 * diff / before.Freq(group)

    t = []
    for group in sorted(before.Values()):
        for pmf in [before, after]:
            t.append(pmf.Freq(group))

    SummarizeImpact(*t)


def SummarizeImpact(women, new_women, men, new_men):
    """Prints a summary of the impact of a change.

    women, new_women: number of women in the sample before and after
    men, new_men: number of men in the sample before and after
    """
    print
    print 'impact'

    # compute the swing in a field the size of Boston
    old_field = ComputeField(men, women)
    new_field = ComputeField(new_men, new_women)
    print 'old field', old_field
    print 'new field', new_field

    swing = new_field[0] - old_field[0]
    print 'Change in number of men', swing

    # apply that swing to the baseline numbers from Boston
    base_men, base_women = 4651, 4951
    print 'before', Fraction(base_men, base_women)

    new_men = base_men + swing
    new_women = base_women - swing 
    print 'after', Fraction(new_men, new_women)


def Fraction(x, y):
    """Returns x, y, percent of x in (x+y)."""
    return x, y, 100.0 * x / (x + y)


def EvaluateImpact():
    """Evaluate the impact of the proposed changes to the standard."""
    res = []
    offsets = [0, -1, -6, -11, -16]
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
        SummarizeImpact(women, new_women, men, new_men)

    
def main():
    #EvaluateBurfoot()

    # MakeGraphs()

    EvaluateImpact()




if __name__ == '__main__':
    main()
