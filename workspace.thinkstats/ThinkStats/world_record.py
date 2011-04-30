"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2011 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import csv
import datetime
import math
import random
import sys

import matplotlib.pyplot as pyplot

import Cdf
import correlation
import myplot
import Pmf
import thinkstats

def ReadData(filename='Marathon_world_record_times.csv', speed=False):
    """Reads a CSV file 

    Args:
      filename: string filename

    Returns:
      list of ...
    """
    fp = open(filename)
    reader = csv.reader(fp)
    reader.next()
    
    distances = {}

    while True:
        distance, gender, data = ReadDistance(reader, speed)
        if distance == None:
            break
        distances[distance, gender] = data

    return distances

def ParseDate(date):
    date = date.split('[')[0]

    formats = ['%m/%d/%Y', '%B %d, %Y', '%b %d, %Y', '%d %B %Y', '%Y-%m-%d']

    for format in formats:
        try:
            return datetime.datetime.strptime(date, format)
        except ValueError:
            continue

    print 'Unparsed date:', date
    return None

def ParseTime(time):
    t = time.split(':')
    seconds = float(t[-1])
    
    try:
        minutes = int(t[-2])
    except IndexError:
        minutes = 0

    try:
        hours = int(t[-3])
    except IndexError:
        hours = 0

    minutes = hours * 60.0 + minutes + seconds / 60.0

    if minutes > 1000:
        print time, seconds, minutes, hours

    return minutes

miles = {
    'marathon' : 26.21875,
    'half marathon' : 13.109375,
    '10000m' : 6.21371192,
    '5000m' : 3.10685596,
    '3000m' : 1.86411358,
    'mile' : 1,
    '1500m' : 0.932056788,
    '800m' : 0.497096954,
    '400m' : 0.248548477,
    '200m' : 0.124274238,
    '100m' :  0.0621371192,
}


def ReadDistance(reader, speed):
    try:
        t = reader.next()
        distance, gender = t
    except StopIteration:
        return None, None, None
    except ValueError:
        print t
        return None, None, None

    data = []
    for t in reader:
        if len(t) == 0:
            break

        time, date = t[0], t[3]
        minutes = ParseTime(time)

        date_obj = ParseDate(date)
        if date_obj is None:
            print t
            print date
        dayofyear = int(date_obj.strftime('%j'))
        years = date_obj.year + dayofyear / 365.24

        if speed:
            speed = miles[distance] / (minutes / 60)
            data.append((years, speed))
        else:
            data.append((years, minutes))

    data.sort()
    return distance, gender, data

def Logistic(z):
    return 1 / (1 + math.exp(-z))

def GeneratePerson(n=10):
    factors = [random.normalvariate(0.0, 1.0) for i in range(n)]
    logs = [Logistic(x) for x in factors]
    return min(logs)

def WorldRecord(m=1000000, n=10, fast=False):
    pmf = Pmf.Pmf()
    data = []
    best = 0.0
    for i in xrange(m):
        person = GeneratePerson(n)
        if not fast:
            pmf.Incr(person)

        if person > best:
            best = person
            data.append((float(i+1)/m, best))
            
    return pmf, data

def PlotTimes(distances, plot_gender='male'):

    for (distance, gender), data in distances.iteritems():
        if gender != plot_gender:
            continue

        pyplot.clf()
        xs, ys = zip(*data)
        pyplot.plot(xs, ys, 'o:')

        root = 'world_record_%s' % distance
        myplot.Plot(root=root,
                    xlabel='year',
                    ylabel='minutes',
                    title='%s world record progression' % distance)

def PlotSpeeds(distances, plot_gender='male'):

    pyplot.rcdefaults()
    pyplot.rc('figure', figsize=(4, 10))
    pyplot.rc('font', size=9.0)
    pyplot.rc('xtick.major', size=0)
    pyplot.rc('ytick.major', size=0)

    pyplot.subplots_adjust(wspace=0.4, hspace=0.4, 
                           right=0.95, left=0.1,
                           top=0.95, bottom=0.05)

    t = miles.items()
    t.sort(key=lambda x: x[1])
    titles = [x[0] for x in t]

    gender = plot_gender
    i=0
    for distance in titles:
        if distance == '200m':
            continue
        i += 1

        data = distances[distance, gender]
        if gender != plot_gender:
            continue

        pyplot.subplot(6, 2, i)
        xs, ys = zip(*data)

        # extend the current record to the present
        first_x = xs[1]
        last_x = xs[-1]
        pyplot.xticks([int(first_x), 1960, 2011])

        last_y = ys[-1]
        pyplot.plot([last_x, 2011.4], [last_y, last_y], 'b:')

        pyplot.plot(xs, ys, 'o:')

        pyplot.title(distance)

    root = 'world_record_speed'
    myplot.Save(root=root)

"""There are six billion ways not to be the fastest marathoner in the world."""

def main(script):
    #distances = ReadData(speed=False)
    #PlotTimes(distances)

    distances = ReadData(speed=True)
    PlotSpeeds(distances)

    #PlotCdfs()
    #PlotSimulations()

    #PlotMarathon(distances['marathon', 'male'])

def PlotCdfs():
    """Plots distribution of ability for different number of factors.
    
    After 100000 people:
    n   max value
    50  0.333842852938
    10  0.6483317765470
    5   0.837633976492
    1   0.983619459771

    """
    cdfs = []
    for n in [50, 10, 5, 1]:
        pmf, data = WorldRecord(m=100000, n=n)
        cdf = Cdf.MakeCdfFromPmf(pmf, name='n=%d' % n)
        print n, max(cdf.Values())
        cdfs.append(cdf)

    options = dict(linewidth=2)
    plot_options = [options] * len(cdfs)

    myplot.Cdfs(cdfs, root='world_record_cdfs', plot_options=plot_options)


def PlotSimulations():
    pyplot.rcdefaults()
    pyplot.rc('figure', figsize=(4, 4.5))
    pyplot.rc('font', size=9.0)
    pyplot.rc('xtick.major', size=0)
    pyplot.rc('ytick.major', size=0)

    pyplot.subplots_adjust(wspace=0.4, hspace=0.4, 
                           right=0.95, left=0.1,
                           top=0.95, bottom=0.05)

    for i in range(1, 5):
        pyplot.subplot(2, 2, i)
        PlotSimulation(100000)

    myplot.Plot(root='world_record_sim')



def PlotSimulation(m=100000):
    pmf, data = WorldRecord(m=m, n=10, fast=True)
    xs, ys = zip(*data)

    # extend to the present
    last_x = xs[-1]
    last_y = ys[-1]
    pyplot.plot([last_x, 1.0], [last_y, last_y], 'b:')

    pyplot.plot(xs, ys, 'o:')
    pyplot.xscale('log')


def PlotMarathon(data):
    xs, ys = zip(*data)

    n = 29
    fxs, fys = xs[:n], ys[:n]
    sxs, sys = xs[n:], ys[n:]

    inter, slope = correlation.LeastSquares(sxs, sys)
    print '(inter, slope):', inter, slope

    target = 26.21875 / 2
    lxs, lys = MakeLine(target, 0, [1965, 2045.5])
    pyplot.plot(lxs, lys, 'r')

    lxs, lys = MakeLine(inter, slope, [1965, 2045.5])
    pyplot.plot(lxs, lys, 'b')

    #pyplot.plot(fxs, fys, 'go')
    pyplot.plot(sxs, sys, 'bo')

    xticks = range(1980, 2060, 20)
    pyplot.xticks(xticks)

    myplot.Plot(root='world_record_predict',
                title='marathon record projection',
                ylabel='mph')


def MakeLine(inter, slope, xs):
    return xs, [inter + slope * x for x in xs]
    

if __name__ == '__main__':
    main(*sys.argv)
