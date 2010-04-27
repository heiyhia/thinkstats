"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib
import matplotlib.pyplot as plt


def Hist(hist):
    xs, fs = hist.Render()
    plt.bar(xs, fs)

    plt.xlabel('x')
    plt.ylabel('frequency')
    plt.title('Histogram')

    Plot('histogram')


def Hists(hists, root, axis=None, show=False):
    plt.clf()

    width = 0.4
    shifts = [-width, 0.0]

    for i, hist in enumerate(hists):
        xs, fs = hist.Render()
        xs = Shift(xs, shifts[i])
        plt.bar(xs, fs, width=width, color='0.80', label=hist.name)

    Plot(root, axis, show)


def Shift(xs, shift):
    return [x+shift for x in xs]


def Plot(root, axis=None, show=False):
    # TODO: make the labelsize work and put a title on these plots
    matplotlib.rc('axes', labelsize='xx-large', titlesize=20.0)
    matplotlib.rc('legend', fontsize=20.0)

    if axis:
        plt.axis(axis)

    plt.legend(loc=2)

    if show:
        plt.show()

    Save(root, 'eps')
    Save(root, 'png')


def Save(root, format='eps'):
    filename = '%s.%s' % (root, format)
    plt.savefig(filename, format=format)


def Cdfs(cdfs, root, axis=None, show=False):
    plt.clf()

    styles = ['r-', 'b--']
    widths = [3, 2]
    colors = ['0.5', '0.0']

    for i, cdf in enumerate(cdfs):
        
        xs, ps = cdf.Render()
        line = plt.plot(xs, ps, linewidth=widths[i], color=colors[i], label=cdf.name)

    Plot(root, axis, show)

