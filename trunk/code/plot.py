"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib.pyplot as plt

def Hist(hist):
    xs, fs = hist.Render()
    plt.bar(xs, fs)

    plt.xlabel('x')
    plt.ylabel('frequency')
    plt.title('Histogram')

    format = 'png'
    filename = 'figure.' + format
    plt.savefig(filename, format=format)


def Hists(hists, root, format='png', axis=None, show=False):
    plt.clf()

    width = 0.4
    shifts = [-width, 0.0]

    for i, hist in enumerate(hists):
        xs, fs = hist.Render()
        xs = Shift(xs, shifts[i])
        plt.bar(xs, fs, width=width, color='0.80', label=hist.name)

    Plot(root, format, axis, show)


def Shift(xs, shift):
    return [x+shift for x in xs]


def Plot(root, format='png', axis=None, show=False):
    if axis:
        plt.axis(axis)

    plt.legend()

    if show:
        plt.show()

    filename = '%s.%s' % (root, format)
    plt.savefig(filename, format=format)


def Cdfs(cdfs, root, format='png', axis=None, show=False):
    plt.clf()

    styles = ['r-', 'b--']

    for i, cdf in enumerate(cdfs):
        
        xs, ps = cdf.Render()
        plt.plot(xs, ps, styles[i], label=cdf.name)

    Plot(root, format, axis, show)

