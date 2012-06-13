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

order = ['prot', 'cath', 'jew', 'other', 'none', 'NA']
long_order = []


class Respondent(object):

    divide_prot = False

    convert = dict(compwt=float)

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

    denoms = dict(
        bap = range(10, 20),
        meth = range(20, 30),
        luth = range(30, 40),
        pres = range(40, 50),
        episc = range(50, 60),
        )

    def clean(self):
        self.compwt = float(self.compwt)
        self.relig_name = self.lookup_religion(self.relig, self.denom)
        self.parelig_name = self.lookup_religion(self.parelig, self.paden)
        self.marelig_name = self.lookup_religion(self.marelig, self.maden)
        self.relig16_name = self.lookup_religion(self.relig16, self.denom16)

    def lookup_religion(self, relig, denom):
        relname = relig

        if relig in self.religions:
            relname = self.religions[relig]

        if self.divide_prot and relig == 1:
            for denom_name, codes in self.denoms.iteritems():
                if denom in codes:
                    relname = denom_name

        return relname


def make_trans(objs, attr1, attr2):
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

    n = len(order)
    matrix = np.empty(shape=(n, n), dtype=np.float)

    for i, x in enumerate(order):
        for j, y in enumerate(order):
            percent = trans[x].Prob(y)
            matrix[i][j] = percent

    return np.transpose(matrix)


def make_pmf(objs, attr):
    pmf = Pmf.Pmf()
    for obj in objs:
        val = getattr(obj, attr)
        wt = obj.compwt
        pmf.Incr(val, wt)
    pmf.Normalize()
    return pmf


def print_pmf_sorted(pmf):
    pvs = [(prob, val) for val, prob in pmf.Items()]
    pvs.sort(reverse=True)
    for prob, val in pvs:
        print val, prob * 100


def print_pmf(pmf):
    for x in order:
        print '%s\t%0.1f' % (x, pmf.Prob(x) * 100)


def pmf_to_vector(pmf):
    t = [pmf.Prob(x) for x in order]
    return np.array(t)


def print_vector(vector):
    for x in order:
        print x, '\t',
    
    print
    for i, x in enumerate(order):
        percent = vector[i] * 100
        print '%0.1f\t' % percent,
    print


def step(matrix, vector):
    new = np.dot(matrix, vector)
    return new


def main(script):
    objs = columns.read_csv('gss1.csv', Respondent)
    pmf = make_pmf(objs, 'relig_name')
    #print_pmf(pmf)

    vector = pmf_to_vector(pmf)


    trans = make_trans(objs, 'marelig_name', 'relig_name')
    print_trans(trans, order)

    matrix = trans_to_matrix(trans, order)

    print_vector(vector)

    for i in range(10):
        vector = step(matrix, vector)
        print_vector(vector)


if __name__ == '__main__':
    import sys
    main(*sys.argv)

