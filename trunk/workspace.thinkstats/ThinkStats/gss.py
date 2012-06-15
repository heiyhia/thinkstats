"""This file contains code related to "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import numpy as np
import matplotlib.pyplot as pyplot
import myplot
import csv

import Pmf
import Cdf
import columns
import glm

import correlation
import math
import random
import thinkstats

import rpy2.robjects as robjects
r = robjects.r


class Respondent(object):
    """Represents a survey respondent.

    Attributes are set in columns.read_csv.
    """ 
    # map from field name to conversion function
    convert = dict(compwt=float)

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
        8:'NA',
        9:'NA',
        99:'NA'
        }

    # code table for Protestant denominations
    denoms = dict(
        bap = range(10, 20),
        meth = range(20, 30),
        luth = range(30, 40),
        pres = range(40, 50),
        episc = range(50, 60),
        )

    def clean(self):
        """Recodes some attributes.

        Invoked on each object in columns.read_csv
        """
        self.compwt = float(self.compwt)
        self.relig_name = self.lookup_religion(self.relig, self.denom)
        self.parelig_name = self.lookup_religion(self.parelig, self.paden)
        self.marelig_name = self.lookup_religion(self.marelig, self.maden)
        self.relig16_name = self.lookup_religion(self.relig16, self.denom16)

        if self.age > 89:
            self.yrborn = 'NA'
        else:
            self.yrborn = self.year - self.age

    def lookup_religion(self, relig, denom):
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


class Survey(object):
    def __init__(self, rs=None):
        if rs is None:
            self.rs = {}
        else:
            self.rs = rs
        self.cdf = None

    def len(self):
        return len(self.rs)

    def respondents(self):
        return self.rs.itervalues()

    def lookup(self, caseid):
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

    def partition_by_yrborn(self, attr, bin_size=10):
        """Partition the sample by binning birthyear.

        attr: which attribute to collect
        bin_size: number of years in each bin

        Returns: map from index year to Pmf of values
        """
        d = {}
        for r in self.respondents():
            if r.yrborn == 'NA':
                continue

            index = int(r.yrborn / bin_size) * bin_size
            if index not in d:
                d[index] = Pmf.Pmf()
            val = getattr(r, attr)
            d[index].Incr(val, r.compwt)

        for pmf in d.itervalues():
            pmf.Normalize()

        return d

    def count_partition(self, d, val):
        """Returns a time series of probabilities for the given value.

        d: map from index year to Pmf of values
        val: which value to select
        """
        rows = []
        for year, pmf in sorted(d.iteritems()):
            p = pmf.Prob(val)
            rows.append((year, p))
        return rows

    def regress_by_yrborn(self, attr, val):
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

        return xs

    def linear_model(self, xs, print_flag=False):
        res = glm.run_model('y ~ x', print_flag=print_flag)
        estimates = glm.get_coeffs(res)

        inter = estimates['(Intercept)'][0]
        slope = estimates['x'][0]

        xs = np.array(sorted(set(xs)))
        log_odds = inter + slope * xs
        odds = np.exp(log_odds)
        ps = odds / (1 + odds)

        fit = []
        for x, p in zip(xs, ps):
            fit.append((x+1900, p))

        return slope, inter, fit

    def quadratic_model(self, xs, print_flag=False):
        res = glm.run_model('y ~ x + x2', print_flag=print_flag)
        estimates = glm.get_coeffs(res)

        inter = estimates['(Intercept)'][0]
        slope = estimates['x'][0]
        slope2 = estimates['x2'][0]

        xs = np.array(sorted(set(xs)))
        log_odds = inter + slope * xs + slope2 * xs**2
        odds = np.exp(log_odds)
        ps = odds / (1 + odds)

        fit = []
        for x, p in zip(xs, ps):
            fit.append((x+1900, p))

        return slope, inter, fit

    def age_cohort(self, val, start, end):
        resampled = self.resample()
        xs = resampled.regress_by_yrborn('relig_name', val)
        slope, inter, fit = resampled.linear_model(xs)

        cohort = self.resample()

        series = []
        for delta in range(start, end+1):

            total = 0
            count = 0
            for r in cohort.respondents():
                year = r.year + delta
                fake_yrborn = year - r.age
                p = fit_prob(fake_yrborn, slope, inter)

                total += 1
                if random.random() <= p:
                    count += 1

            fraction = float(count) / total
            series.append((year, fraction))

        return series

    def simulate_aging_cohort(self, val, start, end, n=20):
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
        pyplot.plot(xs, means, color='blue', lw=3, alpha=0.5)

        # plot the real data
        series = read_time_series()
        data = get_series_for_val(series, val)
        xs, ps = zip(*data)
        pyplot.plot(xs, ps, color='red', lw=3, alpha=0.5)

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
        random.seed(19)

        # plot some resampled fits
        all_ps = {}
        all_rows = {}
        for i in range(40):
            resampled = self.resample()

            # collect the partitioned estimates
            d = resampled.partition_by_yrborn('relig_name')
            rows = resampled.count_partition(d, val)
            for x, p in rows:
                all_rows.setdefault(x, []).append(p)

            # collect the resampled values
            xs = resampled.regress_by_yrborn('relig_name', val)
            slope, inter, fit = resampled.linear_model(xs)
            for x, p in fit:
                all_ps.setdefault(x, []).append(p)

        plot_interval(all_ps, color='0.9')

        # plot the real fit
        xs = self.regress_by_yrborn('relig_name', val)
        slope, inter, fit = self.linear_model(xs)
        xs, ps = zip(*fit)
        pyplot.plot(xs, ps, lw=3, color='blue', alpha=0.5)

        # plot the real data with error bars
        d = self.partition_by_yrborn('relig_name')
        rows = self.count_partition(d, val)
        xs, ps = zip(*rows[1:-1])

        plot_errorbars(all_rows, lw=1, color='red', alpha=0.5)
        pyplot.plot(xs, ps, marker='s', markersize=8, 
                    lw=0, color='red', alpha=0.5)

        axes = dict(
            none=[1895, 1965, 0, 0.16],
            prot=[1895, 1965, 0, 1],
            cath=[1895, 1965, 0, 0.5],
            jew=[1895, 1965, 0, 0.2],
            other=[1895, 1965, 0, 0.2],
            )

        # make the figure
        myplot.Save(root='gss1',
                    xlabel='Year born',
                    ylabel='Prob of relig=%s' % val,
                    axis=axes[val])


def fit_prob(x, slope, inter):
    log_odds = inter + slope * (x-1900)
    odds = math.exp(log_odds)
    p = odds / (1 + odds)
    return p
    

def make_trans(rs, attr1, attr2):
    """Makes a transition table.

    Returns map from attr1 to normalized Pmf of outcomes.

    rs: list of rects
    attr1: explanatory variable
    attr2: dependent variable
    """
    trans = {}

    for r in rs:
        x = getattr(r, attr1)
        y = getattr(r, attr2)
        wt = r.compwt

        trans.setdefault(x, Pmf.Pmf()).Incr(y, wt)

    for x, pmf in trans.iteritems():
        pmf.Normalize()

    return trans


def print_trans(trans, order):
    """Prints a transition table.

    trans: map from explanatory values to Pmf of outcomes.
    order: string category names in desired order
    """
    print '\t',
    for y in order:
        print y, '\t',
    print

    for x in order:
        print x, '\t', 
        for y in order:
            percent = trans[x].Prob(y) * 100
            print '%0.1f\t' % percent,
        print


def trans_to_matrix(trans, order):
    """Converts a transition table to a matrix.

    trans: map from explanatory values to Pmf of outcomes.
    order: string category names in desired order.

    Returns: numpy array of float
    """
    n = len(order)
    matrix = np.empty(shape=(n, n), dtype=np.float)

    for i, x in enumerate(order):
        for j, y in enumerate(order):
            percent = trans[x].Prob(y)
            matrix[i][j] = percent

    return np.transpose(matrix)


def print_pmf_sorted(pmf):
    """Prints the values in the Pmf in descending order of prob.

    pmf: Pmf object
    """
    pvs = [(prob, val) for val, prob in pmf.Items()]
    pvs.sort(reverse=True)
    for prob, val in pvs:
        print val, prob * 100


def print_pmf(pmf, order):
    """Prints the Pmf with values in the given order.

    pmf: Pmf object
    order: string category names in desired order.
    """
    for x in order:
        print '%s\t%0.1f' % (x, pmf.Prob(x) * 100)


def pmf_to_vector(pmf, order):
    """Converts a Pmf to a vector of probabilites.

    pmf: Pmf object
    order: string category names in desired order.

    Returns: numpy array of float
    """
    t = [pmf.Prob(x) for x in order]
    return np.array(t)


def print_vector(vector, order, head_flag=True):
    """Prints a 1-D numpy array.

    vector: numpy array
    order: string category names
    head_flag: boolean whether to print a header line
    """
    if head_flag:
        for x in order:
            print x, '\t',
        print

    for i, x in enumerate(order):
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
    
    def __init__(self, order):
        self.order = order

    def run_nonlinear(self, matrix, vector):
        """Runs the nonlinear model where the number of conversions depends
        on the prevalance of each category.

        matrix: numpy transition matrix
        vector: numpy state vector
        """
        conversions = matrix / vector

        state = vector
        print_vector(state, self.order, True)

        for i in range(10):
            trans = conversions * state
            state = step(trans, state)
            normalize_vector(state)
            print_vector(state, self.order, False)


    def run_linear(self, matrix, vector, steps=1):
        """Runs the linear model where the rate of conversions is constant.

        matrix: numpy transition matrix
        vector: numpy state vector
        """
        print_vector(vector, self.order, True)

        for i in range(steps):
            vector = step(matrix, vector)
            print_vector(vector, self.order, False)

        return vector


    def display_time_series(self, series):
        years = series.keys()
        years.sort()

        ys = np.zeros(len(years))

        rows = []
        for name in self.order:
            if name == 'NA':
                continue

            for i, year in enumerate(years):
                percent = series[year].Prob(name) * 100
                ys[i] += percent

            print name, ys
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

        myplot.Save(show=True,
                    legend=True,
                    axis=[1972, 2010, 0, 100.5])
            

    def display_changes(self, series,
                        low, high, 
                        change_flag=True,
                        predictions=None):
        years = series.keys()
        years.sort()
        years = [year for year in years if low <= year <= high]

        ys = np.zeros(len(years))

        rows = []
        for name in self.order:
            if name == 'NA':
                continue

            ys = [series[year].Prob(name) * 100  for year in years]
            rows.append(ys)

        colors = ['orange', 'green', 'blue', 'yellow', 'red']
        alphas = [0.5,      0.5,      0.5,    0.8,      0.5]
        markers = ['o', 'o', 'o', 'o', 'o']

        for i in range(len(rows)):
            ys = rows[i]
            if change_flag:
                baseline = ys[0]
                ys = [100.0 * y / baseline for y in ys]
                axis = [low-1, high+1, 50, 260]
            else:
                axis = [low-1, high+1, 0, 90]

            pyplot.plot(years, ys,
                        label=self.order[i],
                        linewidth=3,
                        color=colors[i],
                        alpha=alphas[i])

            if predictions is not None:
                pred = 100.0 * 100.0 * predictions[i] / baseline
                print high, pred
                pyplot.plot(high, pred, 
                            marker=markers[i],
                            linewidth=3,
                            markersize=10,
                            color=colors[i],
                            alpha=alphas[i])

        myplot.Save(show=True,
                    legend=True,
                    axis=axis)


def read_time_series(filename='GSS_relig_time_series.csv'):
    """Reads data from CSV file and returns map from year to Pmf of relig.

    filename: string

    Returns: map from year to Pmf of outcomes
    """
    fp = open(filename)
    reader = csv.reader(fp)

    header1 = reader.next()[1:]
    header2 = reader.next()[1:]

    series = {}
    for t in reader:
        year = int(t[0])
        total = float(t[-1])
        row = [float(x)/total for x in t[1:-1]]
        pmf = combine_row(row, header2)

        # normalizing shouldn't be necessary, but the totals tend
        # to be off in the third decimal place, so I'm cleaning that up
        pmf.Normalize()

        series[year] = pmf

    fp.close()

    return series

def get_series_for_val(series, val):
    res = []
    for year, pmf in sorted(series.iteritems()):
        p = pmf.Prob(val)
        res.append((year, p))
    return res


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


def make_time_series():
    series = read_time_series()
    model = Model(order)

    for year, pmf in sorted(series.iteritems()):
        print year
        for name, prob in pmf.Items():
            print name, prob


def make_model():
    rs = columns.read_csv('gss1988.csv', Respondent)
    pmf = make_pmf(rs, 'relig_name')
    #print_pmf(pmf)

    order = ['prot', 'cath', 'jew', 'other', 'none', 'NA']

    vector = pmf_to_vector(pmf, order)

    trans = make_trans(rs, 'marelig_name', 'relig_name')
    print_trans(trans, order)

    matrix = trans_to_matrix(trans, order)

    print
    model = Model(order)

    predictions = model.run_linear(matrix, vector, steps=1)

    #model.display_time_series(series)
    #model.display_changes(series, 1972, 1988)
    model.display_changes(series, 1988, 2010, predictions=predictions)
    return


def plot_interval(all_ps, **options):
    xs = all_ps.keys()
    xs.sort()
    columns = [all_ps[x] for x in xs]
    stats = [thinkstats.MeanVar(ys) for ys in columns]
    min_ps = [mu - 2 * math.sqrt(var) for mu, var in stats]
    max_ps = [mu + 2 * math.sqrt(var) for mu, var in stats]
    mean_ps = [mu for mu, var in stats]

    pyplot.fill_between(xs, min_ps, max_ps, linewidth=0, **options)
    return xs, mean_ps


def plot_errorbars(all_ps, **options):
    xs = all_ps.keys()
    xs.sort()

    lows = []
    highs = []
    for x in xs:
        col = all_ps[x]
        col.sort()
        low = col[1]
        high = col[-2]
        lows.append(low)
        highs.append(high)

    for x, low, high in zip(xs, lows, highs):
        pyplot.plot([x, x], [low, high], **options)


def main(script):
    survey = Survey()
    survey.read_csv('gss1988.csv', Respondent)

    val = 'none'

    xs = survey.regress_by_yrborn('relig_name', val)
    slope, inter, fit = survey.linear_model(xs, True)
    #slope, inter, fit = survey.quadratic_model(xs, True)

    survey.plot_relig_vs_yrborn(val)
    survey.simulate_aging_cohort(val, -16, 23)

if __name__ == '__main__':
    import sys
    main(*sys.argv)

