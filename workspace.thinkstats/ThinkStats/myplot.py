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


class InfiniteList(list):
    def __init__(self, val):
        self.val = val

    def __getitem__(self, index):
        return self.val


def Hist(hist, root=None, bar_options={}, **options):
    """Plots a histogram.

    Args:
      hist: Hist or Pmf object
      root: string filename root
      bar_options: dictionary of options passed to pylot.bar
      options: dictionary of options
    """
    pyplot.clf()

    xs, fs = hist.Render()
    pyplot.bar(xs, fs, align='center', **bar_options)

    Plot(root, **options)


def Pmf(pmf, root=None, plot_options={}, **options):
    """Plots a PMF as a line.

    Args:
      pmf: Hist or Pmf object
      root: string filename root
      bar_options: dictionary of options passed to pylot.plot
      options: dictionary of options
    """
    pyplot.clf()

    xs, fs = pmf.Render()
    pyplot.plot(xs, fs, **plot_options)

    Plot(root, **options)


def Hists(hists, root=None, bar_options=InfiniteList(dict()), **options):
    """Plots two histograms.

    Args:
      hists: list of two Hist or Pmf objects
      root: string filename root
      bar_options: sequence of option dictionaries
      options: dictionary of options
    """
    pyplot.clf()

    width = 0.4
    shifts = [-width, 0.0]

    for i, hist in enumerate(hists):
        xs, fs = hist.Render()
        xs = Shift(xs, shifts[i])
        pyplot.bar(xs, fs, label=hist.name, width=width, **bar_options[i])

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


def Cdfs(cdfs, 
         root=None, 
         line_options=InfiniteList(dict()), 
         complement=False,
         **options):
    """Plots a sequence of CDFs.
    
    Args:
      cdfs: sequence of CDF objects
      root: string root of the filename to write
      line_options: sequence of option dictionaries
      complement: boolean, whether to plot the complementary CDF
      options: dictionary of keyword options passed along to Plot
    """
    pyplot.clf()

    styles = options.get('styles', None)
    if styles is None:
        styles = InfiniteList('-')

    for i, cdf in enumerate(cdfs):
        
        xs, ps = cdf.Render()
        if complement:
            ps = [1.0-p for p in ps]
            
        line = pyplot.plot(xs, ps,
                           styles[i],
                           label=cdf.name,
                           **line_options[i]
                           )

    Plot(root, **options)


def Plot(root=None, formats=None, **options):
    """Generate plots in the given formats.

    Pulls options out of the option dictionary and passes them to
    title, xlabel, ylabel, xscale, yscale, axis and legend.

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

    yscale = options.get('yscale', 'linear')
    pyplot.yscale(yscale)

    xscale = options.get('xscale', 'linear')
    pyplot.xscale(xscale)

    axis = options.get('axis', None)
    if axis:
        pyplot.axis(axis)

    loc = options.get('loc', 0)
    legend = options.get('legend', True)
    if legend:
        pyplot.legend(loc=loc)

    if formats is None:
        formats = ['eps', 'png']

    if root:
        for format in formats:
            Save(root, format)

    show = options.get('show', False)
    if show:
        pyplot.show()

def Save(root, format='eps'):
    """Writes the current figure to a file in the given format.

    Args:
      root: string filename root

      format: string format
    """
    filename = '%s.%s' % (root, format)
    print 'Writing', filename
    pyplot.savefig(filename, format=format)


