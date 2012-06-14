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


def make_trans(objs, attr1, attr2):
    """Makes a transition table.

    Returns map from attr1 to normalized Pmf of outcomes.

    objs: list of objects
    attr1: explanatory variable
    attr2: dependent variable
    """
    trans = {}

    for obj in objs:
        x = getattr(obj, attr1)
        y = getattr(obj, attr2)
        wt = obj.compwt

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


def make_pmf(objs, attr):
    """Make a PMF for an attribute.  Uses compwt to weight respondents.

    objs: list of Respondents
    attr: string attr name

    Returns: normalized PMF
    """
    pmf = Pmf.Pmf()
    for obj in objs:
        val = getattr(obj, attr)
        wt = obj.compwt
        pmf.Incr(val, wt)
    pmf.Normalize()
    return pmf


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
    vector *= total / np.sum(vector)


def run_nonlinear(matrix, vector, order):
    conversions = matrix / vector

    state = vector
    print_vector(state, order, True)

    for i in range(10):
        trans = conversions * state
        state = step(trans, state)
        normalize_vector(state)
        print_vector(state, order, False)


def run_linear(matrix, vector, order):
    print_vector(vector, order, True)

    for i in range(2):
        vector = step(matrix, vector)
        print_vector(vector, order, False)


def read_time_series(filename='GSS_relig_time_series.csv'):
    fp = open(filename)
    reader = csv.reader(fp)

    header1 = reader.next()[1:]
    header2 = reader.next()[1:]
    print header2

    series = {}
    for t in reader:
        year = int(t[0])
        total = float(t[-1])
        row = [float(x)/total for x in t[1:-1]]
        pmf = combine_row(row, header2)
        print year, total
        for name, prob in pmf.Items():
            print name, prob

        # normalizing shouldn't be necessary, but the totals tend
        # to be off in the third decimal place, so I'm cleaning that up
        pmf.Normalize()

        series[year] = pmf

    fp.close()

    return series


def combine_row(row, header2):
    pmf = Pmf.Pmf()
    pmf.Incr('NA', 0)
    for name, prob in zip(header2, row):
        pmf.Incr(name, prob)
    return pmf


def main(script):
    read_time_series()
    return

    order = ['prot', 'cath', 'jew', 'other', 'none', 'NA']
    long_order = []

    objs = columns.read_csv('gss1.csv', Respondent)
    pmf = make_pmf(objs, 'relig_name')
    #print_pmf(pmf)

    vector = pmf_to_vector(pmf, order)

    trans = make_trans(objs, 'marelig_name', 'relig_name')
    print_trans(trans, order)

    matrix = trans_to_matrix(trans, order)

    print
    run_linear(matrix, vector, order)

if __name__ == '__main__':
    import sys
    main(*sys.argv)

