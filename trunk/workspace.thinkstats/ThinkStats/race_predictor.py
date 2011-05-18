import lxml.html
import math
import mechanize
import random
import shelve
import time

import Cdf
import correlation
import matplotlib.pyplot as pyplot
import myplot
import thinkstats


class Browser(object):

    def __init__(self):
        """Initializes the Browser."""
        br = mechanize.Browser()
        br.set_handle_robots(False)

        br.addheaders = [('User-agent', 
                          'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1)'
                          'Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

        self.br = br

    def FindBostonResult(self, first, last):
        """Looks up the given runner and make a record with runner's info."""

        base_url = 'http://registration.baa.org/2010/cf/Public/'
        starting_url = base_url + 'iframe_ResultsSearch.cfm'

        br.open(starting_url)

        root = lxml.html.fromstring(br.response().read())

        # find the page's form
        br.select_form(name='PublicSearch')
        br.form.set_all_readonly(False)

        # set the fields and submit
        br['FirstName'] = first
        br['LastName'] = last
        br.submit()

        # get the response and pull out the right table
        root = lxml.html.fromstring(br.response().read())
        tables = root.cssselect('table')
        table = tables[4]

        record = self.ParseTable(table)

        return record

    def ParseTable(self, table):
        """Pulls the data out of the table and puts it in a record.

        Returns an empty record if the name is not found.
        """
        attrs = ['bib', 'name', 'age', 'gender', 'city', 'state', 'country', 
                 'ctz', 'unk1', 'unk2', '5k', '10k', '15k', '20k', 'half', 
                 '25k', '30k', '35k', '40k', 'pace', 'proj', 'net', 'overall', 
                 'gender', 'div']

        cells = table.cssselect('td')
        record = {}

        for attr, cell in zip(attrs, cells):
            if attr == 'name':
                # unpack the name from the <a> tag
                try:
                    cell = cell.cssselect('a')[0]
                    name = cell.text
                except IndexError:
                    name = cell.text

            record[attr] = cell.text.strip()

        if record['bib']:
            return record
        else:
            return {}



"""
Place Name                 Age Sex No.  City            St Net Time  Pace  Gun T
ime  
===== ==================== === === ==== ================== ========= ===== =====
==== 
    1 DERESE DENIBOB        28 M   2470 BRONX NY           1:05:20.8  5:00 1:05:
21.4 
"""

def CleanLine(line, half=False):
    """Converts a line from coolrunning results to a tuple of values."""
    t = line.split()

    if len(t) < 11:
        return None
    
    net, pace, gun = t[-3:]

    for time in [gun, net, pace]:
        if ':' not in time:
            return None

    name = line[6:26]

    age = line[27:30]
    try:
        age = int(age)
    except ValueError:
        return None

    gender = line[31]

    if gender not in ['M', 'F']:
        return None

    return name, gender, age, net


def ReadResults(filename='NewBedford2010.html'):
    """Reads results from coolrunning and returns a list of tuples."""
    results = []
    for line in open(filename):
        t = CleanLine(line)
        if t:
            results.append(t)
    return results


def ConvertTimeToMinutes(time):
    """Converts time in HH:MM:SS format to minutes."""
    t = time.split('.')
    if len(t) == 2:
        time, fraction = time.split('.')
    
    t = [int(x) for x in time.split(':')]
    if len(t) == 2:
        h = 0
        m, s = t
    else:
        h, m, s = t

    mins = h * 60 + m + s / 60.0
    return mins


def ReadShelf(shelf):
    """Reads data from the shelf and returns a list of pairs.

    Each pair is a half-marathon time and marathon time in minutes.
    """
    unchecked = 0
    no_record = 0
    bad_age = 0
    good = 0

    print 'runners:', len(shelf)

    pairs = []
    for key, t in shelf.iteritems():
        if t is None:
            unchecked += 1
            continue

        first, last, gender, age, net, record = t
        if not record:
            no_record += 1
            continue

        age2 = int(record['age'])
        if age2 < age or age2 > age+1:
            bad_age += 1
            continue

        good += 1

        half = ConvertTimeToMinutes(net)
        full = ConvertTimeToMinutes(record['net'])
        pairs.append((half, full))

        # print key, age, record['age'], net, record['net'], half, full

    print 'unchecked', unchecked
    print 'no record', no_record
    print 'wrong person', bad_age
    print 'good', good

    return pairs


def LookupResults(shelf):
    """Read results from New Bedform and find runners in Boston."""
    # read New Bedford results
    res = ReadResults()

    # try to find the same people in the Boston results
    for name, gender, age, net in res:
        key = name.lower()

        # if we have already looked up this name, skip it
        if key in shelf and shelf[key] != None:
            print 'skipping', key
            continue

        t = key.split()
        first = t[0]
        last = ' '.join(t[1:])

        print first, last, gender, age, net

        record = {}
        record = FindBostonResult(first, last)
        
        print record

        shelf[key] = [first, last, gender, age, net, record]

        # wait a few seconds so we don't annoy them
        delay = random.randint(1, 5)
        print delay
        time.sleep(delay)


def Fit(halfs, fulls):
    """Find the linear least squares fit between halfs and fulls."""
    inter, slope = correlation.LeastSquares(halfs, fulls)
    print '(inter, slope):', inter, slope

    res = correlation.Residuals(halfs, fulls, inter, slope)
    R2 = correlation.CoefDetermination(fulls, res)

    print 'inter', inter
    print 'slope', slope
    print 'R^2', R2
    print

    print 'prediction', inter + slope * ConvertTimeToMinutes('1:34:05')

    return inter, slope, R2


def FitOneParam(halfs, fulls):
    res = [y-x for x, y in zip(halfs, fulls)]
    mu = thinkstats.Mean(res)
    alpha = mu / math.log10(2.0)
    print 'mu, alpha', mu, alpha
    return mu, alpha


def FilterPairs(pairs, low, high, factor):
    """Select only pairs where the half marathon time is in range
    and the marathon time <= factor * half-marathon-time."""
    return [(x, y) for x, y in pairs if low <= x <= high and y <= x*factor]


def LogTransform(pairs):
    return [(math.log10(x), math.log10(y)) for x, y in pairs]


def CountQualifiers(pairs, threshold=200.0):
    """Count people who ran the marathon faster than threshold.

    Returns a tuple of number of qualifiers, number of runners.
    """
    qual = [y for x, y in pairs if y <= threshold]
    return len(qual), len(pairs)


def ScatterPlot(root, xs, ys, fit=None, log_flag=False, alpha=1.0):
    """Makes a scatter plot of ys versus xs."""
    if fit:
        fxs = [min(xs), max(xs)]
        if len(fit) == 3:
            inter, slope, R2 = fit
            fys = [inter + slope * x for x in fxs]
        elif len(fit) == 1:
            fys = [fit[0] for x in fxs]

        pyplot.plot(fxs, fys, color='red')

    pyplot.scatter(xs, ys, alpha=alpha, edgecolors='none')

    if log_flag:
        xlabel = 'Half marathon (log10 min)'
        ylabel = 'Marathon (log10 min)'
    else:
        xlabel = 'Half marathon (min)'
        ylabel = 'Marathon (min)'

    myplot.Plot(root=root,
                xlabel=xlabel,
                ylabel=ylabel,
                title='Boston vs New Bedford',
                legend=False,
                show=True)
        

def MakeScatter(shelf, root='race_predictor1', low=None, high=None,
                fit_flag=False, log_flag=False, line_flag=False):
    pairs = ReadShelf(shelf)

    if low:
        pairs = FilterPairs(pairs, low=low, high=high,
                            factor=100)

    if log_flag:
        pairs = LogTransform(pairs)

    print 'qualifiers', CountQualifiers(pairs)
    halfs, fulls = zip(*pairs)

    fit = Fit(halfs, fulls)
    if log_flag:
        mu, alpha = FitOneParam(halfs, fulls)
        inter = mu
        slope = 1.0
        fit = inter, slope, None

    if line_flag:
        fit = [200.0]

    if fit_flag:
        ScatterPlot(root, halfs, fulls, fit, log_flag=log_flag)
    else:
        ScatterPlot(root, halfs, fulls, log_flag=log_flag)


def MakeFigures(shelf):
    MakeScatter(shelf)
    MakeScatter(shelf, 
                root='race_predictor2',
                fit_flag=True,
                log_flag=True)

    time = ConvertTimeToMinutes('1:34:05')
    spread = 2.0

    MakeScatter(shelf, 
                root='race_predictor3',
                low=time-spread,
                high=time+spread,
                fit_flag=True,
                line_flag=True)


def MakePercentiles(shelf, n=50):
    pairs = ReadShelf(shelf)
    pairs.sort()

    xs = []
    plists = []

    for i in range(0, len(pairs), n):
        subset = pairs[i : i+n]
        print i, len(subset)
        halfs, fulls = zip(*subset)
        cdf = Cdf.MakeCdfFromList(fulls)
        ys = [cdf.Percentile(x) for x in [5, 25, 50, 75, 95]]
        x = thinkstats.Mean(halfs)

        print x, ys

        xs.append(x)
        plists.append(ys)

    # drop the last point
    xs.pop()
    plists.pop()

    ylists = zip(*plists)

    plot_options = [
        dict(color='red', label='5%ile', linestyle='dotted'),
        dict(color='orange', label='25%ile', linestyle='dashed'),
        dict(color='yellow', label='50%ile', linestyle='solid'),
        dict(color='green', label='75%ile', linestyle='dashed'),
        dict(color='cyan', label='95%ile', linestyle='dotted'),
        ]

    pyplot.plot([94, 94], [100, 350])

    for ys, d in zip(ylists, plot_options):
        pyplot.plot(xs, ys, linewidth=3, **d)

    myplot.Plot(root='race_predictor4',
                xlabel='Half marathon (min)',
                ylabel='Marathon (min)',
                show=True)
        
        


def main():
    shelf = shelve.open('race_predictor.db')

    try:
        # LookupResults(shelf)
        MakePercentiles(shelf)
    finally:
        shelf.close()


if __name__ == '__main__':
    main()
