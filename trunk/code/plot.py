"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib
import matplotlib.pyplot as pyplot

# customize some matplotlib attributes
matplotlib.rc('figure', figsize=(5, 4))

#matplotlib.rc('font', size=20.0)
#matplotlib.rc('axes', labelsize=22.0, titlesize=22.0)
#matplotlib.rc('legend', fontsize=20.0)

#matplotlib.rc('xtick.major', size=6.0)
#matplotlib.rc('xtick.minor', size=3.0)

#matplotlib.rc('ytick.major', size=6.0)
#matplotlib.rc('ytick.minor', size=3.0)



def Hist(hist, root, **options):
    """Plots a histogram.

    Args:
      hist: Hist or Pmf object

      root: string filename root

      options: dictionary of options
    """
    pyplot.clf()

    xs, fs = hist.Render()
    pyplot.bar(xs, fs, align='center')

    Plot(root, **options)


def Hists(hists, root, **options):
    """Plots two histograms.

    Args:
      hists: list of two Hist or Pmf objects

      root: string filename root

      options: dictionary of options
    """
    pyplot.clf()

    width = 0.4
    shifts = [-width, 0.0]
    colors = ['0.5', '0.8']

    for i, hist in enumerate(hists):
        xs, fs = hist.Render()
        xs = Shift(xs, shifts[i])
        pyplot.bar(xs, fs, width=width, color=colors[i], label=hist.name)

    Plot(root, **options)


def Shift(xs, shift):
    """Adds a constant to a sequence of values.

    Args:
      xs: sequence of values

      shift: value to add

    Returns:
      sequence of numbers
    """
    return [x+shift for x in xs]


def Plot(root, formats=None, **options):
    """Generate plots in the given formats.

    Args:
      root: string filename root

      formats: list of string formats

      options: dictionary of options
    """
    title = options.get('title', '')
    pyplot.title(title)

    xlabel = options.get('xlabel', '')
    pyplot.xlabel(xlabel)

    ylabel = options.get('ylabel', '')
    pyplot.ylabel(ylabel)

    axis = options.get('axis', None)
    if axis:
        pyplot.axis(axis)

    loc = options.get('loc', 0)
    legend = options.get('legend', True)
    if legend:
        pyplot.legend(loc=loc)

    show = options.get('show', False)
    if show:
        pyplot.show()

    if formats is None:
        formats = ['eps', 'png']

    for format in formats:
        Save(root, format)


def Save(root, format='eps'):
    """Writes the current figure to a file in the given format.

    Args:
      root: string filename root

      format: string format
    """
    filename = '%s.%s' % (root, format)
    print 'Writing', filename
    pyplot.savefig(filename, format=format)


def Cdfs(cdfs, root, **options):
    pyplot.clf()

    styles = options.get('styles', None)
    if styles is None:
        styles = InfiniteList('-')

    widths = [2, 2]
    colors = ['0.0', '0.0']

    for i, cdf in enumerate(cdfs):
        
        xs, ps = cdf.Render()
        line = pyplot.plot(xs, ps,
                        styles[i],
                        linewidth=widths[i], 
                        color=colors[i],
                        label=cdf.name)

    Plot(root, **options)

class InfiniteList(list):
    def __init__(self, val):
        self.val = val

    def __getitem__(self, index):
        return self.val
