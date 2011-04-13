"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from copy import copy

import sys
import gzip
import math
import os

import matplotlib.pyplot as pyplot

import Cdf
import myplot
import Pmf
import table
import thinkstats


class Records(table.Table):
    """Represents the record table."""

    def ReadRecords(self, data_dir=None, n=None):
        if data_dir == None:
            data_dir = './SEER/incidence'
        filename = self.GetFilename()
        self.ReadFile(data_dir, filename, self.GetFields(), table.Record, n)
        self.Recode()

    def GetFilename(self):
        return 'yr1973_2007.seer9/BREAST.TXT.gz'

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.

        Codes:
        sex: 1=male, 2=female
        age: years or 999
        seqno: 0 = only one primary in lifetime
        diagmo: 1-12
        diagyr: 1973+
        site:
        follow: 2 or 4 = active followup
        behavior: 3 = malignant
        race:
        hisp:
        stage: 
        survival: YYMM or 9999
        cause:
        status: 1=alive, 4=dead
        deathclass: 0=considered "alive"  1=death due to primary tumor

        The Survival Time Recode is calculated using the date of
        diagnosis and one of the following: date of death, date last known to
        be alive, or follow-up cutoff date used for this file (see title page
        for date for this file). Thus a person diagnosed in May 1976 and who
        died in May 1980 has a Survival Time Recode of 04 years and 00 months.

        """
        return [
            ('caseid', 1, 8, int),
            ('sex', 24, 24, int),
            ('age', 25, 27, int),
            ('seqno', 35, 36, int),
            ('diagmo', 37, 38, int),
            ('diagyr', 39, 42, int),
            ('site', 43, 46, int),
            ('follow', 181, 181, int),
            ('behavior', 215, 215, int),
            ('race', 224, 224, int),
            ('hisp', 225, 225, int),
            ('stage', 229, 230, int),
            ('survival', 241, 244, str),
            ('cause', 245, 249, int),
            ('status', 255, 255, int),
            ('deathclass', 264, 264, int),
            ]

    def Recode(self):
        """Computes recoded variables."""
        for r in self.records:
            # convert survival time in YYMM to decimal years
            years, months = int(r.survival[:2]), int(r.survival[2:])
            assert months < 12
            r.interval = years + months / 12.0

            r.diagdate = r.diagyr + r.diagmo / 12.0

    def MakeHists(self):
        """Makes a histogram for each attribute."""
        for field in self.GetFields():
            attr = field[0]
            if attr not in ['caseid', 'cause', 'survival']:
                self.MakeHist(attr)

    def MakeHist(self, attr):
        """Makes a histogram for the given attribute and prints it."""
        vals = [getattr(record, attr) for record in self.records]
        hist = Pmf.MakeHistFromList(vals)

        print attr
        for val, freq in hist.Items():
            print val, freq
        print

    def Filter(self):
        """Makes a new table with a subset of the records.

        Selects malignant tumors, patients with only one primary tumor
        in their lifetimes, and cases with active follow-up.
        """
        table = copy(self)
        table.records = [r for r in table.records if (
                r.behavior == 3 and r.seqno == 0 and r.follow in [2, 4]
                )]
        return table

    def FilterAge(self, low, high):
        """Makes a new table with only records in the given age range."""
        table = copy(self)
        table.records = [r for r in table.records if (
                low <= r.age < high
                )]
        return table

    def FilterDate(self, low, high):
        """Makes a new table with only records in the given age range."""
        table = copy(self)
        table.records = [r for r in table.records if (
                low <= r.diagdate < high
                )]
        return table

    def FilterInterval(self, low):
        """Makes a new table with only records with interval > low."""
        table = copy(self)
        table.records = [r for r in table.records if (
                r.interval >= low
                )]
        return table

    def MeanDiagDate(self):
        xs = [r.diagdate for r in self.records]
        return thinkstats.Mean(xs)

    def ComputeHazard(self, t):
        """Computes the fraction of the cohort that dies at time t.

        If we use deathclass, we get the Net cancer-specific survival.

        If we use status, we get Observed all cause survival

        The right stat for purposes of prognosis is Crude probability of death
        """
        deaths = [r for r in self.records if (
                r.interval == t and r.deathclass == 1
                )]
        d, n = len(deaths), len(self.records)
        return d, n, Fraction(d, n)

    def ComputeSurvivalCurve(self, high=20, factor=3):
        """Estimates the survival curve and event density.

        Uses the Kaplan-Meier product-limit estimator:
        http://www.statsdirect.com/help/survival_analysis/kaplan.htm
        http://tkchen.wordpress.com/2008/09/21/
        kaplan-meier-and-nelson-aalen-estimators/
        """
        intervals = set(r.interval for r in self.records)
        table = copy(self)

        lams = Pmf.Pmf(name='hazard')
        tables = []
        ts, ss = [], []
        prod = 1.0

        for t in sorted(intervals):
            if t > high:
                break

            # find the cohort at risk
            tables.append(table)
            table = table.FilterInterval(t)

            # find the number of deaths and death rate
            d, n, f = table.ComputeHazard(t)

            # find the cumulative probability of survival
            s = prod
            prod *= (1 - f)

            ts.append(t)
            ss.append(s)

            bin = math.floor(t*factor) / factor
            lams.Incr(bin, f)

        return tables, ts, lams, ss


def ComputeProbSurvival(ts, ss, t):
    ps = [1-s for s in ss]
    cdf = Cdf.Cdf(ts, ps)
    s = 1 - cdf.Prob(t)
    return s


def ComputeConditionalSurvival(tables):
    for table in tables:
        _, ts, lams, ss = table.ComputeSurvivalCurve()        
        p5 = ComputeProbSurvival(ts, ss, 5)

def PlotSurvivalCurve(ts, lams, ss):
    # scale lams
    denom = max(lams.Probs())
    lams.Normalize(denom)
    myplot.Pmf(lams, 
               line_options=dict(linewidth=2, linestyle='dashed', color='0.7'))

    pyplot.plot(ts, ss, linewidth=2, color='blue', label='survival')
    myplot.Plot(root='seer1',
                title='',
                xlabel='Survival time (years)',
                ylabel='Probability')


def PlotDiagDates(ts, tables):
    mus = [table.MeanDiagDate() for table in tables]
    pyplot.plot(ts, mus, linewidth=2, color='green')
    myplot.Plot(root='seer2',
                title='',
                xlabel='Survival time (years)',
                ylabel='Mean diagnosis date')


def Fraction(n, m):
    return float(n) / m


def PartitionByAge(table):
    for age in [20, 30, 40, 50, 60, 70]:
        part = table.FilterAge(age, age+10)
        print age, len(part.records)
        tables, ts, lams, ss = part.ComputeSurvivalCurve()
        pyplot.plot(ts, ss, linewidth=2, label=str(age))

    myplot.Plot(root='seer3', show=True,
                xlabel='Survival time (years)',
                ylabel='Probability')


def PartitionByDate(table):
    for date in [1995, 1990, 1985, 1980, 1975, 1970]:
        part = table.FilterDate(date, date+5)
        print date, len(part.records)
        tables, ts, lams, ss = part.ComputeSurvivalCurve()
        pyplot.plot(ts, ss, linewidth=2, label=str(date))

    myplot.Plot(root='seer3', show=True,
                xlabel='Survival time (years)',
                ylabel='Probability')


def main(name, data_dir=None):
    table = Records()
    table.ReadRecords(data_dir=data_dir, n=10000)
    print 'Number of records', len(table.records)

    table = table.Filter()
    print 'Malignant, single primary with follow-up', len(table.records)

    #table = table.FilterAge(30, 39)
    #print 'Age in 30s', len(table.records)

    tables, ts, lams, ss = table.ComputeSurvivalCurve()
    #PlotSurvivalCurve(ts, lams, ss)
    #PlotDiagDates(ts, tables)

    ComputeConditionalSurvival(ts, ss)

    #PartitionByDate(table)

if __name__ == '__main__':
    main(*sys.argv)
