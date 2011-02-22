
from BeautifulSoup import BeautifulSoup

def ParseResults(filename='Chicago2010.html'):
    soup = BeautifulSoup(open(filename))

    tables = soup.findAll("table", attrs={'class':'ResultsTable'})
    table = tables[0]
             

    for row in table.findAll('tr')[1:]:
        cols = row.findAll('td')[6:9]
        t = [col.contents[0].strip() for col in cols]
        print ' '.join(t) 

ParseResults('Chicago2009.html')
