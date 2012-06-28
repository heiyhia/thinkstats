"""This file contains code related to "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""
import copy
import numpy as np
import matplotlib.pyplot as pyplot
import myplot
import csv
import re

import HTML

import Pmf
import Cdf
import columns
import glm

import correlation
import math
import random
import thinkstats

import rpy2.robjects as robjects
# r = robjects.r

ORDER_NA = ['prot', 'cath', 'jew', 'other', 'none', 'NA']
ORDER = ['prot', 'cath', 'jew', 'other', 'none']

PROPER_NAME = dict(prot='Protestant', cath='Catholic', jew='Jewish', 
                   other='Other', none='None' )

USE_PARENT_DATA = True

def clean_var(var, na_codes):
    """Replaces invalid codes with NA."""
    if var in na_codes:
        return 'NA'
    else:
        return var


class Respondent(object):
    """Represents a survey respondent.

    Attributes are set in columns.read_csv.
    """ 
    # map from field name to conversion function
    convert = dict(compwt=float, sei=float, masei=float, pasei=float)

    # divide Protestants into denominations?
    divide_prot = False

    # code table for relig and related attributes
    religions = {
        0:'NA',
        1:'prot',
        2:'cath',
        3:'jew',
        4:'none',
        5:'other',
        6:'other',
        7:'other',
        8:'other',
        9:'other',
        10:'other',
        11:'prot',
        12:'other',
        13:'other',
        98:'NA',
        99:'NA',
        }

    # code table for Protestant denominations
    denoms = dict(
        bap = range(10, 20),
        meth = range(20, 30),
        luth = range(30, 40),
        pres = range(40, 50),
        episc = range(50, 60),
        )

    switch_dict = dict(
        prot = (100000, 200000),
        cath = (200000, 300000),
        jew = (300000, 400000),
        none = (400000, 500000),
        other = (500000, 600000),
        )

    switwhy_dict = {
        10: 'marriage',
        20: 'friends',
        30: 'family',
        40: 'location',
        50: 'theological',
        77: 'positive beliefs',
        0: 'NA',
        99: 'NA',
        }

    # codes for marelkid and parelkid
    relkid_dict = {
        0: 'NA',
        1: 'prot',
        2: 'cath',
        3: 'jew',
        4: 'other',
        5: 'other',
        6: 'other',
        7: 'none',
        8: 'NA',
        9: 'NA',
        }
    
    def clean(self):
        self.clean_relig()

        if self.age > 89:
            self.yrborn = 'NA'
            self.decade = 'NA'
            self.age5 = 'NA'
        else:
            self.yrborn = self.year - self.age
            self.decade = int(self.yrborn / 10) * 10
            self.age5 = int(self.age / 5) * 5

        try:
            self.clean_socioeconomic()
        except AttributeError:
            pass

        try:
            self.clean_children()
        except AttributeError:
            pass

        try:
            self.clean_switch()
        except AttributeError:
            pass

    def clean_children(self):
        # loop through the children to get the years they were born
        for i in range(1, 10):
            attr = 'kdyrbrn%d' % i
            year = getattr(self, attr)
            if year == 0 or year > 9995:
                setattr(self, attr, 'NA')

    def clean_socioeconomic(self):
        self.educ = clean_var(self.educ, [97, 98, 99])
        self.college = 1 if self.educ >= 16 else 0
        self.paeduc = clean_var(self.paeduc, [97, 98, 99])
        self.maeduc = clean_var(self.maeduc, [97, 98, 99])

        self.sei = clean_var(self.sei, [-1, 99.8, 99.9])
        self.pasei = clean_var(self.pasei, [-1, 99.8, 99.9])
        self.masei = clean_var(self.masei, [-1, 99.8, 99.9])

        vars = [self.parelkid_name, self.marelkid_name, 
                self.attendpa, self.attendma,
                self.paeduc, self.maeduc, 
                self.pasei, self.masei]
        self.complete = ('NA' not in vars)

    def clean_relig(self):
        self.relig_name = self.lookup_religion(self.relig)
        self.relig16_name = self.lookup_religion(self.relig16)
        self.sprelig_name = self.lookup_religion(self.sprel)

        # for 1988 and before, use parelig and marelig
        # after 1988, use parelkid and marelkid
        if self.year > 1988:
            self.parelig_name = self.relkid_dict[self.parelkid]
            self.marelig_name = self.relkid_dict[self.marelkid]
        else:
            self.parelig_name = self.lookup_religion(self.parelig)
            self.marelig_name = self.lookup_religion(self.marelig)

        self.has_relig = 0 if self.relig_name=='none' else 1
        self.pa_has = 0 if self.parelig_name=='none' else 1
        self.ma_has = 0 if self.marelig_name=='none' else 1
        self.sp_has = 0 if self.sprelig_name=='none' else 1

        # do the parents have the same religion?
        if self.pa_has and self.parelig_name==self.marelig_name:
            self.par_same = 1
        else:
            self.par_same = 0

        # raised in one of the parents' religions?
        if ((self.pa_has and self.parelig_name==self.relig16_name) or
            (self.ma_has and self.marelig_name==self.relig16_name)):
            self.raised = 1
        else:
            self.raised = 0

        # married in the same religion?
        if 'NA' in [self.relig_name, self.sprelig_name]:
            self.married_in = 'NA'
        else:
            if self.has_relig and self.relig_name == self.sprelig_name:
                self.married_in = 1
            else:
                self.married_in = 0

        # TODO: attendpa is available in 2008; in 1988 it was
        # called paattend and the coding was different
        if self.year == 2008:
            try:
                self.attendpa = clean_var(self.attendpa, [0, 10, 98, 99])
                self.attendma = clean_var(self.attendma, [0, 10, 98, 99])
                self.attendkid = (1 if self.attendpa >= 7 or self.attendma >= 7
                              else 0)
            except AttributeError:
                pass

    def clean_lib(self):
        self.lib = self.code_lib(self.relig_name, self.fund)
        self.pa_lib = self.code_lib(self.parelig_name, self.pafund)
        self.ma_lib = self.code_lib(self.marelig_name, self.mafund)

    def clean_switch(self):
        self.switch1 = self.lookup_switch(self.switch1)
        self.switch2 = self.lookup_switch(self.switch2)
        self.switch3 = self.lookup_switch(self.switch3)
        self.switwhy1 = self.lookup_switwhy(self.switwhy1)
        self.switwhy2 = self.lookup_switwhy(self.switwhy2)

    def code_lib(self, relig_name, fund):
        if relig_name == 'none':
            return 4
        if fund in [1,2,3]:
            return fund
        else:
            return 'NA'
                
    def lookup_religion(self, relig, denom=None):
        """Converts religion codes to string names.

        relig: code from relig and related fields
        denom: code from denom and related fields

        Returns: string
        """
        relname = relig

        if relig in self.religions:
            relname = self.religions[relig]

        if self.divide_prot and relig == 1:
            for denom_name, codes in self.denoms.iteritems():
                if denom in codes:
                    relname = denom_name

        return relname

    def lookup_switch(self, switch):
        """Converts religion codes to string names.

        switch: code from one of the switch fields

        Returns: string
        """
        if switch in [0, 999999]:
            return 'NA'

        for name, (low, high) in self.switch_dict.iteritems():
            if low <= switch < high:
                return name

        return '?'

    def lookup_switwhy(self, switwhy):
        """Converts reason codes to text.

        switch: code from one of the switwhy fields

        Returns: string
        """
        try:
            return self.switwhy_dict[switwhy]
        except KeyError:
            return 'other'

    def ages_when_child_born(self):
        """Returns a list of the respondent's age when children were born."""
        ages = []
        for i in range(1, self.childs+1):
            attr = 'kdyrbrn%d' % i
            child_born = getattr(self, attr)
            if child_born == 'NA':
                return 'NA'
            age_when_born = child_born - self.yrborn
            ages.append(age_when_born)

        return ages

    def simulate_all_children(self, birth_model):
        # collect the actual children
        years_born = []
        for i in range(1, self.childs+1):
            attr = 'kdyrbrn%d' % i
            child_born = getattr(self, attr)
            if child_born == 'NA':
                continue
            years_born.append(child_born)

        # simulate fake children
        for i, age in enumerate(range(self.age, 50)):
            p = birth_model.prob(age)
            if random.random() < p:
                child_born = self.year+i
                years_born.append(child_born)

        children = [self.make_child(yb) for yb in years_born]
        return children

    def make_child(self, year_born):
        child = Respondent()
        child.get_next_id()
        child.compwt = self.compwt
        child.year_born = year_born
        
        # TODO: model the spouse's religion to fill in missing data
        if self.sex == 1:
            child.pa_has = self.has_relig
            child.ma_has = self.sp_has
        else:
            child.ma_has = self.has_relig
            child.pa_has = self.sp_has
            
        # do the parents have the same religion?
        if child.pa_has and self.relig_name==self.sprelig_name:
            child.par_same = 1
        else:
            child.par_same = 0

        return child

    def get_next_id(self, t=[90000]):
        self.caseid = t[0]
        t[0] += 1


class Survey(object):
    """Represents a set of respondents as a map from caseid to Respondent."""

    def __init__(self, rs=None):
        if rs is None:
            self.rs = {}
        else:
            self.rs = rs
        self.cdf = None

    def add_respondent(self, r):
        """Adds a respondent to this survey."""
        self.rs[r.caseid] = r

    def len(self):
        """Number of respondents."""
        return len(self.rs)

    def respondents(self):
        """Returns an iterator over the respondents."""
        return self.rs.itervalues()

    def lookup(self, caseid):
        """Looks up a caseid and returns the Respondent object."""
        return self.rs[caseid]

    def read_csv(self, filename, constructor):
        """Reads a CSV file, return the header line and a list of objects.

        filename: string filename
        """
        objs = columns.read_csv(filename, constructor)
        for obj in objs:
            self.rs[obj.caseid] = obj

    def make_pmf(self, attr):
        """Make a PMF for an attribute.  Uses compwt to weight respondents.

        attr: string attr name

        Returns: normalized PMF
        """
        pmf = Pmf.Pmf()
        for r in self.respondents():
            val = getattr(r, attr)
            wt = r.compwt
            pmf.Incr(val, wt)
        pmf.Normalize()
        return pmf

    def make_cdf(self):
        """Makes a CDF with caseids and weights.

        Cdf.Random() selects from this CDF in proportion to compwt
        """
        items = [(caseid, r.compwt) for caseid, r in self.rs.iteritems()]
        self.cdf = Cdf.MakeCdfFromItems(items)

    def partition_by_attr(self, attr):
        """Makes a map from year to Survey.

        attr: string attribute to be used as a key
        """
        surveys = {}
        for r in self.respondents():
            val = getattr(r, attr)
            if val == 'NA':
                continue

            if val not in surveys:
                surveys[val] = Survey()
            surveys[val].add_respondent(r)
        return surveys

    def make_series(self, attr):
        d = {}
        for r in self.respondents():
            val = getattr(r, attr)
            if r.year not in d:
                d[r.year] = Pmf.Pmf()
            d[r.year].Incr(val, r.compwt)

        for pmf in d.itervalues():
            pmf.Normalize()

        return Series(d)

    def cross_tab(self, attr1, attr2):
        table = {}
        hist = Pmf.Hist()

        for name in ORDER:
            table[name] = Pmf.Pmf()

        for r in self.respondents():
            x = getattr(r, attr1)
            y = getattr(r, attr2)
            if x=='NA' or y=='NA':
                continue

            table[x].Incr(y, r.compwt)
            hist.Incr(x, r.compwt)

        normalize_table(table)
        return Table(table, hist)

    def make_religiosity_curves_by_decade(self, relig_name):
        surveys = self.partition_by_attr('decade')

        curves = []
        labels = []
        for decade, survey in sorted(surveys.iteritems()):
            if decade == 1890 or survey.len() < 300:
                continue
            curve = survey.make_religiosity_curve()
            labels.append(str(decade))
            curves.append(curve)

        root = 'gss.religiosity.%s' % relig_name
        title = 'Religiosity curves, %s' % PROPER_NAME[relig_name]

        plot_curves(curves, labels)
        myplot.Save(root=root,
                    title=title,
                    xlabel='Age when surveyed',
                    ylabel='Fraction with any religion',
                    )

    def make_religiosity_curves(self):
        surveys = self.partition_by_attr('relig16_name')

        curves = []
        for name in ORDER:
            survey = surveys[name]
            curve = survey.make_religiosity_curve()
            curves.append(curve)

        plot_relig_curves(curves)
        myplot.Show()

    def make_religiosity_curve(self):
        d = {}
        for r in self.respondents():
            if r.relig_name == 'NA':
                continue

            age = r.age5
            if age not in d:
                d[age] = Pmf.Hist()
                
            d[age].Incr(r.relig_name != 'none')

        rows = []
        for age, hist in sorted(d.iteritems()):
            if hist.Total() < 30:
                continue
            fraction = fraction_true(hist)
            rows.append((age, fraction))

        curve = zip(*rows)
        return curve

    def resample(self, n=None):
        """Form a new cohort by resampling from this survey.

        n: number of respondents in new sample
        """
        if self.cdf is None:
            self.make_cdf()

        n = n or len(self.rs)
        ids = self.cdf.Sample(n)
        rs = dict((i, self.rs[caseid]) for i, caseid in enumerate(ids))
        return Survey(rs)

    def subsample(self, filter_func):
        """Form a new cohort by filtering respondents

        filter_func: function that takes a respondent and returns boolean
        """
        pairs = [(r.caseid, r) for r in self.respondents() if filter_func(r)]
        rs = dict(pairs)
        return Survey(rs)

    def investigate_conversions(self, old, new):
        switches = []

        for r in self.respondents():
            if r.relig16_name != old or r.relig_name != new:
                continue
            
            print r.switch1, r.switch2, r.switch3

    def investigate_switches(self, old, new):
        switches = []

        for r in self.respondents():
            switch1 = Switch(r.switch1, r.switch2,
                             r.switage1, r.switwhy1)
            switch2 = Switch(r.switch2, r.switch3,
                             r.switage2, r.switwhy2)

            if switch1.match(old, new):
                switches.append(switch1)

            if switch2.match(old, new):
                switches.append(switch2)

        for switch in switches:
            print switch.age, switch.why


    def partition_by_yrborn(self, attr, bin_size=10):
        """Partition the sample by binning birthyear.

        attr: which attribute to collect
        bin_size: number of years in each bin

        Returns: map from decade year to Pmf of values
        """
        d = {}
        for r in self.respondents():
            if r.yrborn == 'NA':
                continue

            decade = r.decade
            if decade not in d:
                d[decade] = Pmf.Pmf()
            val = getattr(r, attr)
            d[decade].Incr(val, r.compwt)

        for pmf in d.itervalues():
            pmf.Normalize()

        return d

    def count_partition(self, d, val):
        """Returns a time series of probabilities for the given value.

        d: map from decade year to Pmf of values
        val: which value to select
        """
        rows = []
        for year, pmf in sorted(d.iteritems()):
            p = pmf.Prob(val)
            rows.append((year, p))
        return rows

    def regress_by_yrborn(self, attr, val):
        """Performs a regression on a variable vs year born.

        Logistic regression of the fraction where the given
        attribute has the given value.

        attr: dependent variable
        val: value of the variable

        Returns a Regression object.
        """
        rows = []

        for r in self.respondents():
            if r.yrborn == 'NA':
                continue

            y = 1 if getattr(r, attr) == val else 0
            x = r.yrborn - 1900

            rows.append((y, x))

        ys, xs = zip(*rows)
        x2s = [x**2 for x in xs]
        col_dict = dict(y=ys, x=xs, x2=x2s)
        glm.inject_col_dict(col_dict)

        return Regression(xs)

    def regress_relig(self, model, print_flag=True):
        """Performs a regression.

        model: string model in r format
        print_flag: boolean, whether to print results

        Returns: LogRegression object
        """
        def clean(attr):
            m = re.match('as.factor\((.*)\)', attr)
            if m:
                return m.group(1)
            return attr
                
        # pull out the attributes in the model
        rows = []
        t = model.split()
        attrs = [clean(attr) for attr in model.split() if len(attr)>1]

        for r in self.respondents():
            row = [getattr(r, attr) for attr in attrs]
            rows.append(row)

        rows = [row for row in rows if 'NA' not in row]

        # inject the data and runs the model
        col_dict = dict(zip(attrs, zip(*rows)))
        glm.inject_col_dict(col_dict)

        res = glm.run_model(model, print_flag=print_flag)
        estimates = glm.get_coeffs(res)

        return LogRegression(res, estimates)

    def iterate_respondent_child_ages(self):
        """Loops through respondents and generates (respondent, ages) pairs.
        
        Where ages is the list of ages at which this parent had children.

        Skips parents with unknown year of birth or any children with
        unknown year of birth.
        """
        for r in self.respondents():
            if r.yrborn == 'NA':
                continue

            ages = r.ages_when_child_born()
            if ages == 'NA':
                continue

            yield r, ages

    def plot_child_curves(self):
        """Makes a plot showing child curves for parent's decade of birth."""
        d = {}
        for r, ages in self.iterate_respondent_child_ages():
            for age in range(13, r.age):
                if (r.decade, age) not in d:
                    d[r.decade, age] = Pmf.Hist()
                # record whether this person had a child at this age
                d[r.decade, age].Incr(age in ages)

        table = np.zeros(shape=(8,90), dtype=np.float)
        for (decade, age), hist in sorted(d.iteritems()):
            index = (decade-1900)/10
            table[index, age] = fraction_true(hist)

        self.child_table = table

        decades, all_ages = zip(*d.iterkeys())
        decades = set(decades)
        ages = [age for age in set(all_ages) if age < 50]
        ages.sort()

        options = dict(lw=3, alpha=0.5)

        for decade in sorted(decades):
            if decade < 1930:
                continue
            label = str(decade)
            index = (decade-1900)/10
            ys = np.cumsum([table[index, age] for age in ages])
            myplot.Plot(ages, ys, label=label, **options)

        myplot.Save(root='gss4',
                    xlabel='Age of parent',
                    ylabel='Cumulative number of children',
                    )

    def plot_child_curve(self):
        model = self.make_birth_model()
        ages = [age for age in model.all_ages if age < 50]
        ages.sort()

        table = model.table
        ps = [table[age] for age in ages]
        ys = np.cumsum(ps)
        myplot.Plot(ages, ys, color='purple', 
                    lw=3, alpha=0.5, linestyle='dashed', label='model')

        myplot.Save(root='gss5',
                    xlabel='Age of parent',
                    ylabel='Cumulative number of children',
                    )

    def make_birth_model(self):
        """Makes a plot showing child curves for parent's decade of birth
        and the aggregated model."""
        d = {}
        for r, ages in self.iterate_respondent_child_ages():
            if r.decade < 1940:
                continue

            # loop through the ages we know about for this respondent
            for age in range(13, r.age):
                if age not in d:
                    d[age] = Pmf.Hist()
                # record whether this person had a child at this age
                d[age].Incr(age in ages)

        table = np.zeros(shape=(90), dtype=np.float)
        for age, hist in sorted(d.iteritems()):
            yes, no = hist.Freq(True), hist.Freq(False)
            table[age] = float(yes) / (yes+no)

        all_ages = set(d.iterkeys())
        return BirthModel(all_ages, table)

    def simulate_generation(self, birth_model):
        all_children = {}
        for r in self.respondents():
            children = r.simulate_all_children(birth_model)
            for child in children:
                print child
                all_children[child.caseid] = child

        next_gen = Survey(all_children)
        return next_gen

    def age_cohort(self, val, start, end):
        """Runs one simulation of the aging cohort.

        val: which religion name to track
        start: low end of the range of year to age by
        end: high end of the range of year to age by

        Returns: a time series of (year, fraction) pairs
        """
        # resample and estimate a linear model
        resampled = self.resample()
        reg = resampled.regress_by_yrborn('relig_name', val)
        fit = reg.linear_model()

        # resample again before aging
        cohort = self.resample()

        # loop through the years and accumulate results
        series = []
        for delta in range(start, end+1):

            total = 0
            count = 0
            for r in cohort.respondents():
                year = r.year + delta
                fake_yrborn = year - r.age
                p = reg.fit_prob(fake_yrborn)

                total += 1
                if random.random() <= p:
                    count += 1

            fraction = float(count) / total
            series.append((year, fraction))

        return series

    def simulate_aging_cohort(self, val, start, end, n=20):
        """Simulates the aging of the cohort for one year
 
        Generates a plot of the results.

        val: which religion name to track
        start: low end of the range of year to age by
        end: high end of the range of year to age by
        n: how many simulations to run
        """
        pyplot.clf()
        random.seed(17)

        # run the simulation
        all_ps = {}
        for i in range(n):
            series = self.age_cohort(val, start, end)
            for x, p in series:
                all_ps.setdefault(x, []).append(p)

        # plot the simulated data
        xs, means = plot_interval(all_ps, color='0.9')
        myplot.Plot(xs, means, color='blue', lw=3, alpha=0.5)

        # plot the real data
        series = Series()
        data = get_series_for_val(series, val)
        xs, ps = zip(*data)
        myplot.Plot(xs, ps, color='red', lw=3, alpha=0.5)

        axes = dict(
            none=[1968, 2011, 0, 0.16],
            prot=[1968, 2011, 0, 1],
            cath=[1968, 2011, 0, 0.5],
            jew=[1968, 2011, 0, 0.2],
            other=[1968, 2011, 0, 0.2],
            )

        myplot.Save(root='gss2',
                    xlabel='Year of survey',
                    ylabel='Fraction with relig=%s' % val,
                    axis=axes[val]
                    )

    def plot_relig_vs_yrborn(self, val):
        """Makes a plot of religious preference by year.

        val: string, which religion name to track.
        """
        random.seed(19)

        # plot some resampled fits
        all_ps = {}
        all_rows = {}
        for i in range(4):
            print i
            resampled = self.resample()

            # collect the partitioned estimates
            d = resampled.partition_by_yrborn('relig_name')
            rows = resampled.count_partition(d, val)
            for x, p in rows:
                all_rows.setdefault(x, []).append(p)

            # collect the resampled values
            reg = resampled.regress_by_yrborn('relig_name', val)
            fit = reg.linear_model()
            for x, p in fit:
                all_ps.setdefault(x, []).append(p)

        plot_interval(all_ps, color='0.9')

        # plot the real fit
        reg = self.regress_by_yrborn('relig_name', val)
        fit = reg.linear_model()
        xs, ps = zip(*fit)
        myplot.Plot(xs, ps, lw=3, color='blue', alpha=0.5)

        # plot the real data with error bars
        d = self.partition_by_yrborn('relig_name')
        rows = self.count_partition(d, val)
        xs, ps = zip(*rows[1:-1])

        plot_errorbars(all_rows, lw=1, color='red', alpha=0.5)
        myplot.Plot(xs, ps, marker='s', markersize=8, 
                    lw=0, color='red', alpha=0.5)

        axes = dict(
            none=[1875, 1995, 0, 0.4],
            prot=[1895, 1965, 0, 1],
            cath=[1895, 1965, 0, 0.5],
            jew=[1895, 1965, 0, 0.2],
            other=[1865, 1995, 0, 0.5],
            )

        # make the figure
        myplot.Save(root='gss1',
                    xlabel='Year born',
                    ylabel='Prob of relig=%s' % val,
                    axis=axes[val])

    def run_transition_model(self, model_factory, n=20):
        before = self.print_state_vector()

        data = []
        for i in range(n):
            trans_model = model_factory()
            next_gen = trans_model.next_gen(self, resample_flag=True)
            after = next_gen.print_state_vector(head_flag=False)
            data.append(after)

        preds = self.extract_predictions(data)
        return preds

    def extract_predictions(self, data):
        """Gets predictions from simulation results."""
        cols = zip(*data)

        preds = []
        for col in cols:
            mean = thinkstats.Mean(col)
            col = list(col)
            col.sort()
            span = col[1], col[-2]
            preds.append((mean, span))

        return preds
            
    def print_state_vector(self, attr='relig_name', head_flag=True):
        """Prints the state of the given attribute.

        attr: string attribute name
        head_flag: boolean, whether to print the header line

        Returns: np state array
        """
        pmf = self.make_pmf(attr)
        vector = pmf_to_vector(pmf)
        print_vector(vector, head_flag)
        return vector


class TransitionModel(object):
    def __init__(self, survey, resample_flag=False):
        """Makes and runs the transition model.

        Returns a vector of predictions.

        survey: Survey object
        resample_flag: boolean, whether to resample
        """
        if resample_flag:
            survey = survey.resample()

        self.spouse_table = SpouseTable(survey)
        self.env_table = EnvironmentTable(survey)
        self.trans_table = TransitionTable(survey, 'Transition table')

    def next_gen(self, survey, resample_flag=False):
        # generate the next generation
        next_gen = self.simulate_transition(survey, resample_flag)

        # compute the generation table (from parent to child religion)
        self.gen_table = TransitionTable(next_gen, 'Generation table',
                                    'fake_parent_relig_name', 'relig_name')

        return next_gen

    def simulate_transition(self, survey, resample_flag):
        """Simulates one generation.

        resample: boolean, whether to resample
        spouse_table: map from r's religion to Pmf of spouse's religion
        env_table: map from (ma, pa) religion to "raised religion"
        trans_table: map from raised religion to r's religion

        Returns: Survey object with next generation
        """
        if resample_flag:
            survey = survey.resample()

        next_gen = {}

        for r in survey.respondents():
            if r.relig_name == 'NA':
                continue

            # choose a random spouse
            sprelig_name = self.spouse_table.generate_spouse(r)

            # choose how to raise the child
            if r.sex == 1:
                raised = self.env_table.generate_raised(sprelig_name, 
                                                        r.relig_name)
            else:
                raised = self.env_table.generate_raised(r.relig_name, 
                                                        sprelig_name)

            # determine the child's religion
            relig_name = self.trans_table.generate_relig(raised)

            # make a Respondent object for the child
            child = copy.copy(r)
            child.fake_parent_relig_name = r.relig_name
            child.relig16_name = raised 
            child.relig_name = relig_name
            next_gen[child.caseid] = child

        return Survey(next_gen)

    def print_model(self):
        print 'prob same', self.spouse_table.prob_same()
        self.spouse_table.print_table()
        self.spouse_table.write_html_table()

        self.env_table.write_combined_table()
        for name in ORDER:
            print 'parent table (%s)' % name
            self.env_table.print_table(name)
            self.env_table.write_html_table(name)

        print 'trans table'
        self.trans_table.print_table()
        self.trans_table.write_html_table()

        print 'gen table'
        self.gen_table.print_table()
        self.gen_table.write_html_table()


class Series(object):
    def __init__(self, d):
        self.d = d

    def print_series(self):
        for year, pmf in sorted(self.d.iteritems()):
            print year,
            for val, prob in sorted(pmf.Items()):
                percent = '%0.0f' % (prob * 100)
                print val, percent, 
            print

    def get_ratios(self, val1, val2):
        rows = []
        for year, pmf in sorted(self.d.iteritems()):
            yes = pmf.Prob(val1)
            no = pmf.Prob(val2)
            if yes + no:
                percent = yes / (yes + no) * 100
                rows.append((year, percent))
            else:
                rows.append((year, None))                

        return rows

    def plot_series(self, rows):
        year, ys = zip(*rows)
        myplot.Plot(year, ys, lw=3, alpha=0.5)
        myplot.Save(root='gss.spouse.series',
            title='Spouses with same religion',
            xlabel='Survey year',
            ylabel='% of respondents')


def make_spouse_series():
    survey = Survey()
    survey.read_csv('gss.series.csv', Respondent)
    
    series = survey.make_series('married_in')
    rows = series.get_ratios(1, 0)
    series.plot_series(rows)


def make_spouse_table(low=1972, high=2010):
    """Makes the spouse tables reading data from the given year.

    year: int

    Returns: SpouseTable
    """
    def filter_func(r):
        return low <= r.year <= high

    survey = Survey()
    survey.read_csv('gss.series.csv', Respondent)
    subsample = survey.subsample(filter_func)
    spouse_table = SpouseTable(subsample)
    return spouse_table


def make_env_table(year):
    """Makes the environment table reading data from the given year.

    year: int

    Returns: EnvironmentTable
    """
    survey = Survey()
    filename = 'gss%d.csv' % year
    survey.read_csv(filename, Respondent)
    env_table = EnvironmentTable(survey)
    return env_table


def make_trans_table(year):
    """Makes the transition table reading data from the given year.

    year: int

    Returns: TransitionTable
    """
    survey = Survey()
    filename = 'gss%d.csv' % year
    survey.read_csv(filename, Respondent)
    trans_table = TransitionTable(survey, 'Transition table')
    return trans_table


class Switch(object):
    def __init__(self, old, new, age, why):
        self.old = old
        self.new = new
        self.age = age
        self.why = why

    def match(self, old, new):
        return self.old==old and self.new==new


class BirthModel(object):
    """Model of the probability of having a child at a given age."""
    def __init__(self, all_ages, table):
        self.all_ages = all_ages
        self.table = table

    def prob(self, age):
        return self.table[age]


class Regression(object):
    """Represents the result of a regression."""
    def __init__(self, xs):
        self.xs = xs

    def linear_model(self, print_flag=False):
        """Runs a linear model and returns fitted values.

        print_flag: boolean, whether to print the R results

        Returns a list of (x, fitted y) pairs
        """
        res = glm.run_model('y ~ x', print_flag=print_flag)
        estimates = glm.get_coeffs(res)

        self.inter = estimates['(Intercept)'][0]
        self.slope = estimates['x'][0]

        xs = np.array(sorted(set(self.xs)))
        log_odds = self.inter + self.slope * xs
        odds = np.exp(log_odds)
        ps = odds / (1 + odds)

        fit = []
        for x, p in zip(xs, ps):
            fit.append((x+1900, p))

        self.fit = fit
        return fit

    def quadratic_model(self, print_flag=False):
        """Runs a quadratic model and returns fitted values.

        print_flag: boolean, whether to print the R results

        Returns a list of (x, fitted y) pairs
        """
        res = glm.run_model('y ~ x + x2', print_flag=print_flag)
        estimates = glm.get_coeffs(res)

        self.inter = estimates['(Intercept)'][0]
        self.slope = estimates['x'][0]
        self.slope2 = estimates['x2'][0]

        xs = np.array(sorted(set(xs)))
        log_odds = self.inter + self.slope * xs + self.slope2 * xs**2
        odds = np.exp(log_odds)
        ps = odds / (1 + odds)

        fit = []
        for x, p in zip(xs, ps):
            fit.append((x+1900, p))

        self.fit = fit
        return fit

    def fit_prob(self, x):
        """Computes the fitted value of y for a given x.

        Only works with the linear model.

        x: float value of x

        Returns: float value of y
        """
        log_odds = self.inter + self.slope * (x-1900)
        odds = math.exp(log_odds)
        p = odds / (1 + odds)
        return p
    

class LogRegression(object):
    def __init__(self, res, estimates):
        self.res = res
        self.estimates = estimates

    def fit_prob(self, r):
        log_odds = 0
        for attr, t in self.estimates.iteritems():
            coef = t[0] 
            if attr == '(Intercept)':
                log_odds += coef
            else:
                x = getattr(r, attr)
                if x == 'NA':
                    return 'NA'
                log_odds += coef * x

        odds = math.exp(log_odds)
        p = odds / (1 + odds)
        return p

    def validate(self, respondents, attr):
        for r in respondents:
            dv = getattr(r, attr)
            p = self.fit_prob(r)
            #print r.caseid, dv, p


def make_regression_model():
    survey = Survey()
    survey.read_csv('gss2008.csv', Respondent)

    complete = survey.subsample(lambda r: r.complete)
    
    attr = 'attendkid'
    print '\nall'
    pmf = survey.make_pmf(attr)
    for val, prob in sorted(pmf.Items()):
        print val, prob

    print '\ncomplete'
    pmf = complete.make_pmf(attr)
    for val, prob in sorted(pmf.Items()):
        print val, prob

    model = 'has_relig ~ pa_has + ma_has + attendkid'
    reg = complete.regress_relig(model, print_flag=True)

    for est, t in reg.estimates.iteritems():
        print est, math.exp(t[0])


def trans_to_matrix(trans):
    """Converts a transition table to a matrix.

    trans: map from explanatory values to Pmf of outcomes.

    Returns: numpy array of float
    """
    n = len(ORDER)
    matrix = np.empty(shape=(n, n), dtype=np.float)

    for i, x in enumerate(ORDER):
        for j, y in enumerate(ORDER):
            percent = trans[x].Prob(y)
            matrix[i][j] = percent

    return np.transpose(matrix)


def fraction_true(hist):
    yes, no = hist.Freq(True), hist.Freq(False)
    return float(yes) / (yes+no)


def print_pmf_sorted(pmf):
    """Prints the values in the Pmf in descending order of prob.

    pmf: Pmf object
    """
    pvs = [(prob, val) for val, prob in pmf.Items()]
    pvs.sort(reverse=True)
    for prob, val in pvs:
        print val, prob * 100


def print_pmf(pmf):
    """Prints the Pmf with values in the given ORDER.

    pmf: Pmf object
    """
    for x in ORDER:
        print '%s\t%0.0f' % (x, pmf.Prob(x) * 100)


def pmf_to_vector(pmf):
    """Converts a Pmf to a vector of probabilites.

    pmf: Pmf object

    Returns: numpy array of float
    """
    t = [pmf.Prob(x) for x in ORDER]
    return np.array(t)


def print_vector(vector, head_flag=True):
    """Prints a 1-D numpy array.

    vector: numpy array
    head_flag: boolean whether to print a header line
    """
    if head_flag:
        for x in ORDER:
            print x, '\t',
        print

    for i, x in enumerate(ORDER):
        percent = vector[i] * 100
        print '%0.1f\t' % percent,
    print


def step(matrix, vector):
    """Advance the simulation by one generation.

    matrix: numpy transition matrix
    vector: numpy state vector

    Returns: new numpy vector
    """
    new = np.dot(matrix, vector)
    return new


def normalize_vector(vector, total=1.0):
    """Normalizes a numpy array to add to total.

    Modifies the array.

    vector: numpy array
    total: float
    """
    vector *= total / np.sum(vector)


class Model(object):
    
    def run_nonlinear(self, matrix, vector):
        """Runs the nonlinear model where the number of conversions depends
        on the prevalance of each category.

        matrix: numpy transition matrix
        vector: numpy state vector
        """
        conversions = matrix / vector

        state = vector
        print_vector(state, True)

        for i in range(10):
            trans = conversions * state
            state = step(trans, state)
            normalize_vector(state)
            print_vector(state, False)

    def run_linear(self, matrix, vector, steps=1):
        """Runs the linear model where the rate of conversions is constant.

        matrix: numpy transition matrix
        vector: numpy state vector
        """
        print_vector(vector, True)

        for i in range(steps):
            vector = step(matrix, vector)
            print_vector(vector, False)

        return vector


class ReligSeries(object):
    def __init__(self, filename='GSS_relig_time_series.csv'):
        """Reads data from CSV file and returns map from year to Pmf of relig.

        filename: string
        """
        fp = open(filename)
        reader = csv.reader(fp)

        header1 = reader.next()[1:]
        header2 = reader.next()[1:]

        self.data = {}
        for t in reader:
            year = int(t[0])
            total = float(t[-1])
            row = [float(x)/total for x in t[1:-1]]
            pmf = combine_row(row, header2)

            # normalizing shouldn't be necessary, but the totals tend
            # to be off in the third decimal place, so I'm cleaning that up
            pmf.Normalize()

            self.data[year] = pmf

        fp.close()

    def plot_time_series_stack(self):
        plot_time_series_stack(self.data)
        myplot.Save(root='gss0',
                    xlabel='Year of survey',
                    ylabel='Fraction of population',
                    legend=True,
                    axis=[1972, 2010, 0, 100.5])
            
    def plot_changes(self, low, high, 
                     change_flag=True,
                     preds=None,
                     spans=None):
        """Makes a plot of changes in religious preference.

        low, high: range of years to plot
        change_flag: boolean: whether to normalize by first year value
        preds: vector of predicted values (should only be used
                     with change_flag=True)
        spans: error ranges for the predictions
        """
        pyplot.clf()
        years = self.get_years(low, high)
        ys = np.zeros(len(years))

        rows = []
        for name in ORDER:
            ys = [self.data[year].Prob(name) * 100  for year in years]
            rows.append(ys)

        colors = ['orange', 'green', 'blue', 'yellow', 'red']
        alphas = [0.5,      0.5,      0.5,    0.8,      0.5]
        stretch = 4 if preds else 1

        for i in range(len(rows)):
            ys = rows[i]
            if change_flag:
                baseline = ys[0]
                ys = [100.0 * y / baseline for y in ys]
                axis = [low-1, high+stretch, 50, 260]
            else:
                axis = [low-1, high+stretch, 0, 90]

            myplot.Plot(years, ys,
                        label=ORDER[i],
                        linewidth=3,
                        color=colors[i],
                        alpha=alphas[i])

            xloc = high + 0.6 * (i+1)

            if spans is not None:
                low_span = 100.0 * 100.0 * spans[i][0] / baseline
                high_span = 100.0 * 100.0 * spans[i][1] / baseline
                myplot.Plot([xloc, xloc], [low_span, high_span],
                            linewidth=3,
                            color=colors[i],
                            alpha=alphas[i])

            if preds is not None:
                pred = 100.0 * 100.0 * preds[i] / baseline
                myplot.Plot(xloc, pred, 
                            marker='s',
                            markersize=10,
                            markeredgewidth=0,
                            color=colors[i],
                            alpha=alphas[i])

        if change_flag:
            if preds is None:
                root = 'gss.change.%d-%d' % (low, high)
            else:
                root = 'gss.pred.%d-%d' % (low, high)
                
            ylabel = '%% change since %d' % low
        else:
            root = 'gss.%d-%d' % (low, high)
            ylabel = '% of respondents'

        myplot.Save(root=root,
                    xlabel='Survey year',
                    ylabel=ylabel,
                    legend=True,
                    axis=axis)

    def get_years(self, low, high):
        years = self.data.keys()
        years.sort()
        years = [year for year in years if low <= year <= high]
        return years

    def test_significance(self, relig_name, low, high):
        """Run a linear regression on market share vs year.

        Prints the results

        series: map from year to Pmf of religions
        relig_name: string, which religion to test
        """
        years = self.get_years(low, high)
        ys = [self.data[year].Prob(relig_name) * 100  for year in years]
        d = dict(years=years, ys=ys)
        glm.inject_col_dict(d)
        glm.linear_model('ys ~ years')

    def get_series_for_val(self, val):
        """Gets the time series for a particular value.

        series: map from year to Pmf of religious preference.
        val: string religion name to track

        Returns: list of (year, fraction) pairs
        """
        res = []
        for year, pmf in sorted(self.data.iteritems()):
            p = pmf.Prob(val)
            res.append((year, p))
        return res


def plot_time_series_stack(data):
    """Makes a plot of the actual data and the model predictions.

    data: map from year to Pmf
    """
    years = data.keys()
    years.sort()

    ys = np.zeros(len(years))

    rows = []
    for name in ORDER:
        for i, year in enumerate(years):
            percent = data[year].Prob(name) * 100
            ys[i] += percent

        rows.append(np.copy(ys))

    colors = ['orange', 'green', 'blue', 'yellow', 'red', '', '', ]

    for i in range(len(rows)-1, -1, -1):
        ys = rows[i]
        if i == 0:
            prev = np.zeros(len(years))
        else:
            prev = rows[i-1]

        pyplot.fill_between(years, prev, ys, 
                            color=colors[i],
                            alpha=0.2)


def combine_row(row, header):
    """Makes a row into a PMF.

    row: list of float data
    header: category each datum should be added to

    Returns: Pmf that maps categories to probs (or fraction of pop)
    """
    pmf = Pmf.Pmf()
    pmf.Incr('NA', 0)
    for name, prob in zip(header, row):
        pmf.Incr(name, prob)
    return pmf


class SpouseTable(object):
    def __init__(self, survey):

        self.sp_female = {}
        self.sp_male = {}
        self.sp_same = {1:Pmf.Hist(), 2:Pmf.Hist()}

        for name in ORDER:
            self.sp_female[name] = Pmf.Pmf()
            self.sp_male[name] = Pmf.Pmf()

        for r in survey.respondents():
            if USE_PARENT_DATA:
                attr1='marelig_name'
                attr2='parelig_name'
            elif r.sex == 1:
                attr1='sprelig_name'
                attr2='relig_name'
            else:
                attr1='relig_name'
                attr2='sprelig_name'

            ma = getattr(r, attr1)
            pa = getattr(r, attr2)
            if ma=='NA' or pa=='NA':
                continue

            self.sp_same[r.sex].Incr(ma==pa)
            self.sp_female[ma].Incr(pa, r.compwt)
            self.sp_male[pa].Incr(ma, r.compwt)

        normalize_table(self.sp_female)
        normalize_table(self.sp_male)

    def print_table(self):
        print 'prob same', self.prob_same()

        for attr in ['sp_male', 'sp_female']:
            print attr
            table = getattr(self, attr)
            print_table(table)

    def prob_same(self):
        t = [fraction_true(self.sp_same[sex]) for sex in [1, 2]]
        return t

    def write_html_table(self):
        for attr in ['sp_male', 'sp_female']:
            table = getattr(self, attr)
            rows, header_row = get_table_rows(table)
            filename = 'gss.%s.html' % attr
            write_html_table(filename, rows, header_row, 'Spouse Table')

    def print_pmf(self, pmf):
        for name in ORDER:
            percent = pmf.Prob(name) * 100
            print '%s %0.0f\t' % (name, percent)

    def generate_spouse(self, r):
        if r.sex == 1:
            pmf = self.sp_male[r.relig_name]
        else:
            pmf = self.sp_female[r.relig_name]
        return pmf.Random()


class EnvironmentTable(object):
    def __init__(self, survey):
        """Makes a 

        Returns map from (marelig, parelig) to normalized Pmf of which
        religion the child is raised in.

        order: string list of relig_names
        """
        self.table = {}
        self.hist = Pmf.Hist()

        for ma in ORDER:
            for pa in ORDER:
                self.table[ma, pa] = Pmf.Pmf()

        for r in survey.respondents():
            ma = r.marelig_name
            pa = r.parelig_name
            raised = r.relig16_name
            if 'NA' in [ma, pa, raised]:
                continue

            self.table[ma, pa].Incr(raised, r.compwt)
            self.hist.Incr((ma, pa))

        normalize_table(self.table)

    def print_table(self, relig_name):
        """Prints a table of mother's religion x father's religion.

        Each row is mother's religion, each column is father's.

        Each entry is the fraction of children raised in relig_name.

        relig_name: string name of religion
        """
        print '\t',
        for y in ORDER:
            print y, '\t',
        print

        for ma in ORDER:
            print ma, '\t', 
            for pa in ORDER:
                pmf = self.table[ma, pa]
                percent = pmf.Prob(relig_name) * 100
                print '%0.0f\t' % percent,
            print

    def write_html_table(self, relig_name):
        """Prints a table of mother's religion x father's religion.

        Each row is mother's religion, each column is father's.

        Each entry is the fraction of children raised in relig_name.

        relig_name: string name of religion
        """
        header_row = ['']
        header_row.extend(ORDER)

        rows = []
        for ma in ORDER:
            row = [ma]
            for pa in ORDER:
                pmf = self.table[ma, pa]
                percent = pmf.Prob(relig_name) * 100
                row.append('%0.0f' % percent)
            rows.append(row)

        filename = 'gss.par.%s.html' % relig_name
        title = 'Parent table (%s)' % relig_name
        write_html_table(filename, rows, header_row, title)

    def print_combined_table(self):
        """Prints a table of mother's religion x father's religion.
        """
        print '\t\t',
        for y in ORDER:
            print y, '\t',
        print

        for ma in ORDER:
            for pa in ORDER:
                print '%4.4s-%4.4s\t' % (ma, pa),
                pmf = self.table[ma, pa]
                for name in ORDER:
                    percent = pmf.Prob(name) * 100
                    print '%0.0f\t' % percent,
                print self.hist.Freq((ma, pa))

    def diff_combined_table(self, other):
        """Prints a table of mother's religion x father's religion.
        """
        print '\t\t',
        for y in ORDER:
            print y, '\t',
        print
        
        total = 0
        for ma in ORDER:
            for pa in ORDER:
                print '%4.4s-%4.4s\t' % (ma, pa),
                pmf = self.table[ma, pa]
                pmf2 = other.table[ma, pa]
                freq = self.hist.Freq((ma, pa))

                for name in ORDER:
                    percent = pmf.Prob(name) * 100
                    percent2 = pmf2.Prob(name) * 100
                    change = percent2 - percent
                    print '%0.0f\t' % percent2,
                    if name == 'none':
                        print '%+0.0f\t' % change,
                        excess = change/100 * freq
                        total += excess
                print '%d\t%0.1f' % (freq, excess)
        print total

    def write_combined_table(self):
        """Prints a table of mother's religion x father's religion.
        """
        header_row = ['parents']
        header_row.extend(ORDER)

        rows = []
        for ma in ORDER:
            for pa in ORDER:
                row = ['%s-%s' % (ma, pa)]
                pmf = self.table[ma, pa]
                for name in ORDER:
                    percent = pmf.Prob(name) * 100
                    row.append('%0.0f' % percent)
                rows.append(row)

        filename = 'gss.parent.html'
        title = 'Parent table'
        write_html_table(filename, rows, header_row, title)

    def generate_raised(self, ma, pa):
        """Chooses a random religion to raise a child in.

        ma: mother's religion
        pa: father's religion

        Returns: string religion name
        """
        pmf = self.table[ma, pa]
        if pmf.Total():
            return pmf.Random()
        else:
            return random.choice([ma, pa])


class TransitionTable(object):
    def __init__(self, survey, title,
                 attr1='relig16_name', attr2='relig_name'):
        """Makes a transition table.

        Returns map from attr1 to normalized Pmf of outcomes.

        attr1: explanatory variable
        attr2: dependent variable
        """
        self.title = title
        self.attr1 = attr1
        self.attr2 = attr2

        self.table = {}
        self.hist = Pmf.Hist()

        for name in ORDER:
            self.table[name] = Pmf.Pmf()

        for r in survey.respondents():
            x = getattr(r, attr1)
            y = getattr(r, attr2)
            if x=='NA' or y=='NA':
                continue

            self.table[x].Incr(y, r.compwt)
            self.hist.Incr(x)

        normalize_table(self.table)

    def print_table(self):
        print_table(self.table)

    def diff_table(self, other):
        """Prints a table of ...
        """
        print '\t',
        for y in ORDER:
            print y, '\t',
        print
        
        total = 0
        for raised in ORDER:
            print raised, '\t',
            pmf = self.table[raised]
            pmf2 = other.table[raised]
            freq = self.hist.Freq(raised)

            for name in ORDER:
                percent = pmf.Prob(name) * 100
                percent2 = pmf2.Prob(name) * 100
                change = percent2 - percent
                print '%0.0f\t' % percent2,
                if name == 'none':
                    print '%+0.0f\t' % change,
                    excess = change/100 * freq
                    total += excess
            print '%d\t%0.1f' % (freq, excess)
        print total

    def write_html_table(self):
        rows, header_row = get_table_rows(self.table)
        filename = 'gss.%s.html' % '.'.join(self.title.lower().split())
        write_html_table(filename, rows, header_row, self.title)

    def generate_relig(self, raised):
        """Chooses a random religious preference.

        raised: string religion raised in

        Returns: string religion name
        """
        pmf = self.table[raised]
        return pmf.Random()


class Table(object):
    def __init__(self, table, hist):
        self.table = table
        self.hist = hist

    def get_pmf(self, key):
        return self.table[key]


def print_table(table):
    """Prints a transition table.

    table: map from explanatory values to Pmf of outcomes.
    """
    print '\t',
    for y in ORDER:
        print y, '\t',
    print

    for x in ORDER:
        print x, '\t', 
        for y in ORDER:
            percent = table[x].Prob(y) * 100
            print '%0.0f\t' % percent,
        print


def get_table_rows(table):
    header_row = ['']
    header_row.extend(ORDER)

    rows = []
    for x in ORDER:
        row = [x]
        for y in ORDER:
            percent = table[x].Prob(y) * 100
            row.append('%0.0f' % percent)
        rows.append(row)

    return rows, header_row


def write_html_table(filename, rows, header_row, title=''):
    """Writes a transition table to a file.

    table: map from explanatory values to Pmf of outcomes.
    """
    print 'Writing', filename
    fp = open(filename, 'w')
    fp.write('<div>\n<h3>%s</h3>\n' % title)

    htmlcode = HTML.table(rows, header_row=header_row)
    fp.write(htmlcode)
    fp.write('</div>\n\n')
 
    fp.close()


def normalize_table(table):
    """Normalize the pmfs in this table."""
    for pmf in table.itervalues():
        if pmf.Total():
            pmf.Normalize()


def make_matrix_model():
    pmf = survey88.make_pmf('relig_name')
    vector = pmf_to_vector(pmf)

    print 'relig16_name'
    trans = survey88.make_trans('relig16_name', 'relig_name')
    print_trans(trans)

    matrix = trans_to_matrix(trans)

    print
    model = Model()
    preds = model.run_linear(matrix, vector, steps=1)

    series = ReligSeries()
    series.plot_changes(1988, 2010, preds=preds)


def plot_time_series():
    """Make time series plots."""
    series = ReligSeries()
    #series.plot_time_series_stack()
    series.plot_changes(1972, 2010, change_flag=False)
    series.plot_changes(1972, 1988)
    series.plot_changes(1988, 2010)


def print_time_series():
    """Print time series data."""
    series = ReligSeries()
    for year, pmf in series.data.iteritems():
        print year
        print_pmf(pmf)


def test_significance():
    series = ReligSeries()
    for name in ORDER:
        print name
        series.test_significance(name, 1972, 2010)
        series.test_significance(name, 1988, 2010)


def plot_interval(all_ps, **options):
    """Plot a 2-standard error interval.

    all_ps: map from x value to list of y values
    options: keyword options passed along to pyplot.fill_between
    """
    xs = all_ps.keys()
    xs.sort()
    columns = [all_ps[x] for x in xs]
    stats = [thinkstats.MeanVar(ys) for ys in columns]
    min_ps = [mu - 2 * math.sqrt(var) for mu, var in stats]
    max_ps = [mu + 2 * math.sqrt(var) for mu, var in stats]
    mean_ps = [mu for mu, var in stats]

    pyplot.fill_between(xs, min_ps, max_ps, linewidth=0, **options)
    return xs, mean_ps


def plot_errorbars(all_ps, n=1, **options):
    """Plot error bars spanning all but n values from the top and bottom.

    all_ps: map from x value to list of y values
    options: keyword options passed along to pyplot.fill_between
    """
    xs = all_ps.keys()
    xs.sort()

    lows = []
    highs = []
    for x in xs:
        col = all_ps[x]
        col.sort()
        low = col[n]
        high = col[-(n+1)]
        lows.append(low)
        highs.append(high)

    for x, low, high in zip(xs, lows, highs):
        myplot.Plot([x, x], [low, high], **options)


def print_transition_model():
    survey88 = Survey()
    survey88.read_csv('gss1988.csv', Respondent)
    trans_model = TransitionModel(survey88)
    next_gen = trans_model.next_gen(survey88)
    trans_model.print_model()


class ModelFactory(object):
    def __init__(self, survey, spouse_flag, env_flag, trans_flag):
        self.survey1 = survey
        self.spouse_flag = spouse_flag
        self.env_flag = env_flag
        self.trans_flag = trans_flag

        if self.any_flag():
            self.survey2 = Survey()
            self.survey2.read_csv('gss2008.csv', Respondent)

    def any_flag(self):
        return self.spouse_flag or self.env_flag or self.trans_flag

    def __call__(self):
        trans_model = TransitionModel(self.survey1,
                                      resample_flag=True)

        if self.any_flag():
            trans_model2 = TransitionModel(self.survey2,
                                           resample_flag=True)

        if self.spouse_flag:
            trans_model.spouse_table = trans_model2.spouse_table

        if self.env_flag:
            trans_model.env_table = trans_model2.env_table

        if self.trans_flag:
            trans_model.trans_table = trans_model2.trans_table

        return trans_model


def run_transition_model(spouse_flag=False, env_flag=False, trans_flag=False):
    """
    
    spouse_flag: boolean, whether to use 2004-2010 spouse tables
    env_flag: boolean, whether to use 2008 parent table
    """
    random.seed(17)

    survey88 = Survey()
    survey88.read_csv('gss1988.csv', Respondent)

    factory = ModelFactory(survey88, spouse_flag, env_flag, trans_flag)
    res = survey88.run_transition_model(factory)
    
    preds, spans = zip(*res)
    for name, pred, span in zip(ORDER, preds, spans):
        print name, '%0.1f (%0.1f %0.1f)' % (pred*100, span[0]*100, span[1]*100)

    series = ReligSeries()
    series.plot_changes(1988, 2010, preds=preds, spans=spans)


def plot_none_vs_yrborn():
    survey = Survey()
    survey.read_csv('gss.series.csv', Respondent)

    def filter_func(r):
        return 1900 <= r.yrborn <= 1990

    subsample = survey.subsample(filter_func)
    subsample.plot_relig_vs_yrborn('none')


def investigate_switches():
    survey = Survey()
    survey.read_csv('gss1988.csv', Respondent)

    surveys = survey.partition_by_attr('relig16_name')
    prot = surveys['prot'].subsample(lambda r: r.relig_name != 'prot')

    rows = []
    for r in prot.respondents():
        row = (r.relig16_name, r.relig_name, 
               r.switch1, r.switwhy1, r.switch2, r.switch3)
        rows.append(row)

    rows.sort()
    for row in rows:
        print row

    return
    for name in ORDER:
        survey = surveys[name]
        survey.make_pmg
        print name, survey.len()

    return
    pmf = survey88.make_pmf('switch1')
    pmf.Set('NA', 0)
    pmf.Normalize()
    for val, prob in sorted(pmf.Items()):
        print val, prob

    survey88.investigate_switches('prot', 'none')


def make_time_series(filename, cutoff=None):
    """Makes a map from decade born to Survey.

    filename: file to read
    cutoff: survey year to cut off results

    Returns: a map from decade born to Survey.
    """
    survey = Survey()
    survey.read_csv(filename, Respondent)

    if cutoff:
        survey = survey.subsample(lambda r: r.year<=cutoff)

    surveys = survey.partition_by_attr('decade')
    return surveys


def make_cross_tabs(surveys, attr1, attr2):
    """Makes a cross tabulation for each year.
    
    surveys: a map from decade born to Survey.
    attr1: string attribute name
    attr2: string attribute name

    Returns: map from year to Table object
    """
    tables = {}
    for year, survey in surveys.iteritems():
        table = survey.cross_tab(attr1, attr2)
        tables[year] = table
    return tables


def plot_upbringing_elements(surveys):
    """Plots elements from upbringing tables as a time series.

    surveys: a map from decade born to Survey.
    """
    tables = make_cross_tabs(surveys, 'parelig_name', 'relig16_name')

    kind = 'upbringing'
    for name in ['prot', 'cath', 'none']:
        plot_element(tables, kind, name)


def plot_transmission_elements(surveys):
    """Plots elements from transmission tables as a time series.

    surveys: a map from decade born to Survey.
    """
    tables = make_cross_tabs(surveys, 'relig16_name', 'relig_name')

    kind = 'transmission'
    for name in ['prot', 'cath', 'none']:
        plot_element(tables, kind, name)


def plot_element(tables, kind, relig_name):
    """Plots an element of a table as a time series.

    tables: map from year to Table object
    kind: string used as part of the output filename
    relig_name: which religion to generate a plot for
    """
    # collect the data
    years = []
    rows = []
    for year, table in sorted(tables.iteritems()):
        if table.hist.Total() < 100:
            continue
        
        row = []
        pmf = table.get_pmf(relig_name)
        for name in ORDER:
            percent = pmf.Prob(name) * 100
            row.append(percent)

        years.append(year)
        rows.append(row)

    # plot it
    cols = zip(*rows)
    plot_simple_series(years, cols)

    root = 'gss.%s.%s' % (kind, relig_name)
    myplot.Save(root=root,
                xlabel='Survey year',
                ylabel='Conversion rate to %s' % relig_name,
                )


def plot_simple_series(years, cols):
    """Plots a set of lines, color coded for religions.

    years: sequence of years
    cols: list of columns, one for each religion, in standard order
    """
    pyplot.clf()

    colors = ['orange', 'green', 'blue', 'yellow', 'red']
    alphas = [0.5,      0.5,      0.5,    0.8,      0.5]

    for i, col in enumerate(cols):
        myplot.Plot(years, col,
                    label=ORDER[i],
                    linewidth=3,
                    color=colors[i],
                    alpha=alphas[i])


def plot_relig_curves(curves):
    """Plots a set of lines, color coded for religions.

    curves: list of (xs, ys) pairs
    """
    pyplot.clf()

    colors = ['orange', 'green', 'blue', 'yellow', 'red']
    alphas = [0.5,      0.5,      0.5,    0.8,      0.5]

    for i, curve in enumerate(curves):
        xs, ys = curve
        myplot.Plot(xs, ys,
                    label=ORDER[i],
                    linewidth=3,
                    color=colors[i],
                    alpha=alphas[i])


def plot_curves(curves, labels):
    """Plots a set of lines.

    curves: list of (xs, ys) pairs
    """
    pyplot.clf()

    for curve, label in zip(curves, labels):
        xs, ys = curve
        myplot.Plot(xs, ys, label=label)


def make_stack_series(tables, name):
    pmfs = {}
    for year, table in tables.iteritems():
        pmf = table.get_pmf(name)
        pmfs[year] = pmf

    plot_time_series_stack(pmfs)
    myplot.Save(root='gss_stack_%s' % name,
                xlabel='Year of survey',
                ylabel='Fraction of population',
                legend=True,
                axis=[1972, 2010, 0, 100.5])


def part_three():
    surveys = make_time_series('gss.series.csv', cutoff=1988)
    plot_transmission_elements(surveys)
    #plot_upbringing_elements(surveys)
    return

    # print the environment tables
    env_table = make_env_table(1988)
    env_table.print_combined_table()

    env_table2 = make_env_table(2008)
    env_table.diff_combined_table(env_table2)

    # print the transmission tables
    trans_table = make_trans_table(1988)
    trans_table.print_table()

    trans_table2 = make_trans_table(2008)
    trans_table2.print_table()
    trans_table.diff_table(trans_table2)

    # plot the model
    run_transition_model(spouse_flag=False, env_flag=False, trans_flag=True)
    return

    # part three
    make_spouse_series()


def part_four():
    survey = Survey()
    survey.read_csv('gss.series.csv', Respondent)
    
    surveys = survey.partition_by_attr('relig16_name')
    for relig_name in ['prot', 'cath', 'none']:
        subsurvey = surveys[relig_name]
        subsurvey.make_religiosity_curves_by_decade(relig_name)

    # sadly, this looks useless
    investigate_switches()


def main(script):
    part_four()
    return

    part_three()
    return

    years = make_time_series()
    tables = make_cross_tabs(years, 'relig_name', 'relig16_name')
    make_stack_series(tables, 'prot')
    return

    plot_none_vs_yrborn()
    return

    make_regression_model()
    return

    spouse_table = make_spouse_table()
    spouse_table.print_table()
    return
    
    print_transition_model()
    return

    print_time_series()
    return

    test_significance()
    return

    plot_time_series()
    return

    series = ReligSeries()
    series.test_significance('none')
    return

    survey94 = Survey()
    survey94.read_csv('gss1994.csv', Respondent)
    #survey94.plot_child_curves()
    #survey94.plot_child_curve()
    birth_model = survey94.make_birth_model()

    survey88 = Survey()
    survey88.read_csv('gss1988.csv', Respondent)
    next_gen = survey88.simulate_generation(birth_model)
    for r in next_gen.respondents():
        print r.caseid, r.compwt
    return
    
    model = 'has_relig ~ pa_has + ma_has + par_same + raised'
    logit = survey88.regress_relig(model)
    logit.validate(survey88.respondents(), 'has_relig')

    pmf = survey88.make_pmf('sprelig_name')
    for val, prob in sorted(pmf.Items()):
        print val, prob

    return

    plot_time_series()
    return
    
    val = 'none'

    reg = survey.regress_by_yrborn('relig_name', val)
    fit = reg.linear_model(True)
    #slope, inter, fit = survey.quadratic_model(xs, True)

    survey.plot_relig_vs_yrborn(val)
    survey.simulate_aging_cohort(val, -16, 23)


if __name__ == '__main__':
    import sys
    main(*sys.argv)

