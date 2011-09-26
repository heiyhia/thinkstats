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

query = """2012/cf/Public/iframe_EntryLists.cfm?mode=results&criteria=&StoredProcParamsOn=yes&&VarAgeLow=0&VarAgeHigh=0&VarGenderID=1&VarBibNumber=&VarLastName=&VarFirstName=&VarStateList=0&VarCountryOfResList=0&VarCountryOfCtzList=0&VarCityList=&VarZipList=&records=25&headerexists=Yes&bordersize=0&bordercolor=%23ffffff&rowcolorone=%23FFCC33&rowcolortwo=%23FFFFFF&headercolor=%23ffffff&headerfontface=Verdana%2CArial%2CHelvetica%2Csans%2Dserif&headerfontcolor=%23004080&headerfontsize=12px&fontface=Verdana%2CArial%2CHelvetica%2Csans%2Dserif&fontcolor=%23000099&fontsize=10px&linkfield=FormattedSortName&linkurl=OpenDetailsWindow&linkparams=RaceAppID&queryname=SearchResults&tablefields=FullBibNumber%2CWaveAndCorral%2CFormattedSortName%2CAgeOnRaceDay%2CGenderCode%2CCity%2CStateAbbrev%2CCountryOfResAbbrev%2CCountryOfCtzAbbrev%2CDisabilityGroup"""


class Browser(object):

    def __init__(self):
        """Initializes the Browser."""
        br = mechanize.Browser()
        br.set_handle_robots(False)

        br.addheaders = [('User-agent', 
                          'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1)'
                          'Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

        self.br = br

    def StartQuery(self):
        base_url = ('http://registration.baa.org/' )

        url = base_url + query
        self.br.open(url)

        resp = self.br.response().read()
        return resp

    def NextQuery(self, start=25):
        forms = self.br.forms()

        self.br.select_form(nr=2)
        self.br.form.set_all_readonly(False)

        self.br['start'] = str(start)
        self.br.submit()

        resp = self.br.response().read()
        return resp

    def PostForm(self, n=10):

        start = 'http://registration.baa.org/2012/cf/Public/iframe_EntryLists.cfm'

        self.br.open(start)
        resp = self.br.response().read()

        root = lxml.html.fromstring(resp)

        self.br.select_form(name='QuickSearch')
        self.br.form.set_all_readonly(False)

        self.br['GenderID'] = ['1']
        self.br.submit()

        resp = self.br.response().read()
        return resp


def ParseRow(row):
    """Pulls the data out of the row and puts it in a record.

    Returns an empty record if the name is not found.
    """
    attrs = ['bib', 'wave', 'name', 'age', 'gender', 'city', 'state', 'country',
             'ctz', 'unk1']

    cells = row.cssselect('td')

    record = {}

    for attr, cell in zip(attrs, cells):
        record[attr] = cell.text.strip()

    if 'name' in record and record['name'] != '':
        print 'name', record['name']
        return record
    else:
        return None


def ParseResp(resp):
    root = lxml.html.fromstring(resp)

    tables = root.cssselect('table')
    table = tables[2]

    subtables = table.cssselect('table')
    subtable = subtables[1]

    rows = subtable.cssselect('tr')
    print len(rows)

    records = []
    for row in rows:
        record = ParseRow(row)
        if record:
            records.append(record)

    return records


def Download(n):
    br = Browser()
    #resp = br.PostForm(n)
    resp = br.RunQuery(n)

    fp = open('scrape.resp.html', 'w')
    fp.write(resp)
    fp.close()


def Read():
    fp = open('scrape.resp.html')
    resp = fp.read()
    fp.close()


def main():
    br = Browser()
    resp = br.StartQuery()
    records = ParseResp(resp)
    print len(records)

    resp = br.NextQuery(26)
    records = ParseResp(resp)
    print len(records)







if __name__ == '__main__':
    main()
