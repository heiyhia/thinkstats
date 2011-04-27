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
    print date

def ParseTime(time):
    try:
        time = datetime.datetime.strptime(time, '%H:%M:%S')
    except ValueError:
        try:
            minutes, seconds = time.split(':')[:2]
            minutes = int(minutes)
            seconds = float(seconds)
        except ValueError:
            minutes = 0
            seconds = float(time)
            
        hours, minutes = divmod(minutes, 60)
        time = datetime.time(hours, minutes, seconds)
    return time

def ReadDistance(reader, speed):
    try:
        distance, gender = reader.next()
    except StopIteration:
        return None, None, None

    data = []
    for t in reader:
        if len(t) == 0:
            break

        time, date = t[0], t[3]
        time = ParseTime(time)
        minutes = time.hour * 60.0 + time.minute + time.second / 60.0

        date_obj = ParseDate(date)
        if date_obj is None:
            print t
            print date
        dayofyear = int(date_obj.strftime('%j'))
        years = date_obj.year + dayofyear / 365.24

        if speed:
            speed = 26.2 / (minutes / 60)
            data.append((years, speed))
        else:
            data.append((years, minutes))

    return distance, gender, data

def Logistic(z):
    return 1 / (1 + math.exp(-z))

def GeneratePerson(n=30):
    factors = [random.normalvariate(0.0, 1.0) for i in range(n)]
    logs = [Logistic(x) for x in factors]
    return min(logs)

def WorldRecord(m=100000):
    pmf = Pmf.Pmf()
    data = []
    best = 0.0
    for i in xrange(m):
        person = GeneratePerson()
        pmf.Incr(person)
        if person > best:
            best = person
            data.append((float(i)/m, best))
            
    return pmf, data

def PlotData(distances, plot_gender='male'):

    for (distance, gender), data in distances.iteritems():
        if gender != plot_gender:
            continue

        pyplot.clf()
        xs, ys = zip(*data)
        pyplot.plot(xs, ys, 'o')

        root = 'world_record_%s' % distance
        myplot.Plot(root=root)

"""There are six billion ways not to be the fastest marathoner in the world."""

def main(script):
    distances = ReadData(speed=False)
    PlotData(distances)


def PlotCdf():
    pmf, data = WorldRecord()
    cdf = Cdf.MakeCdfFromPmf(pmf)

    xs, ys = cdf.Render()
    txs, tys = [], []
    for x, y in zip(xs, ys):
        if y in [0.0, 1.0]:
            continue
        txs.append(x)
        tys.append(math.log(-math.log(y)))

    pyplot.plot(txs, tys)
    myplot.Plot(show=True)
    return

    PlotData(data)
    return


if __name__ == '__main__':
    main(*sys.argv)
