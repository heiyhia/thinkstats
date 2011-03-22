import lxml.html
import mechanize
import random
import shelve
import time

def ParseTable(table):
    """Pulls the data out of the table and puts it in a record.

    Returns an empty record if the name is not found.
    """
    attrs = ['bib', 'name', 'age', 'gender', 'city', 'state', 'country', 'ctz', 
             'unk1', 'unk2', '5k', '10k', '15k', '20k', 'half', '25k', '30k', 
             '35k', '40k', 'pace', 'proj', 'net', 'overall', 'gender', 'div']

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


def FindBostonResult(first, last):
    """Looks up the given runner and make a record with that runner's info."""
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

    record = ParseTable(table)

    return record

br = mechanize.Browser()
br.set_handle_robots(False)

br.addheaders = [('User-agent', 
                  'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1)'
                  'Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

#record = FindBostonResult(first='Ryan', last='Pace')
#print record
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
    """Converts pace in HH:MM:SS format to minutes."""
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
    unchecked = 0
    no_record = 0
    bad_age = 0
    good = 0

    print 'runners:', len(shelf)

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
        print key, age, record['age'], net, record['net']

    print 'unchecked', unchecked
    print 'no record', no_record
    print 'wrong person', bad_age
    print 'good', good


def LookupResults(shelf):
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


def main():
    shelf = shelve.open('race_predictor.db')

    try:
        # LookupResults(shelf)
        ReadShelf(shelf)
    finally:
        shelf.close()


if __name__ == '__main__':
    main()
