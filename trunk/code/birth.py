import sys
import gzip

class Respondent: 
    """Represents a respondent."""

class Pregnancy:
    """Represents a pregnancy."""

# these are the attributes to extract from the respondents file
# these are the attributes to extract from the interval file
intattrs = [
    ['id', 1, 8],
    ['pregordr', 9, 10],
    ['nbrnlv', 14, 14],
    ['wks_preg', 22, 23],
    ['kidssex1', 95, 95],
    ['kidssex2', 151, 151],
    ['outcome', 291, 291],
    ['prglngth', 292, 293],
    ['sex1', 325, 325],
    ['sex2', 326, 326],
]

def process_respondent(line):
# take a line from the respondent file and build a Respondent object
    res = Respondent()
    for (attr, start, end) in resattrs:
        setattr(res, attr, line[start-1:end])
    return res

def process_interval(line):
# take a line from the interval file and build an Interval object
    inter = Interval()
    for (attr, start, end) in intattrs:
        setattr(inter, attr, line[start-1:end])
    return inter


class Table(object):
    """Represents a table."""

    def ReadFile(self, filename, attrs, constructor):
        """Reads a compressed data file builds one object per record.
        """
        fp = gzip.open(filename)
        for line in fp:
            record = self.ProcessLine(line, attrs, constructor)
            self.AddRecord(record)
        fp.close()
    
    def ProcessLine(self, line, attrs, constructor):
        """Reads a line and returns an object with the appropriate attrs.

        Uses GetAttrs to get information about the variables we want
        to extract from this record.

        Args:
            line: string line from a data file

            constructor: callable that makes an object for the record.

        Returns:
            object with appropriate attrs.
        """
        obj = constructor()
        for (attr, start, end, cast) in attrs:
            try:
                val = cast(line[start-1:end])
            except ValueError:
                val = 'NA'
            setattr(obj, attr, val)
        return obj

    def AddRecord(self, record):
        self.records.append(record)

class Respondents(Table):
    def __init__(self, filename='2002FemResp.dat.gz'):
        self.records = []
        self.ReadFile(filename, self.GetAttrs(), Respondent)

    def GetAttrs(self):
        return [
            ('caseid', 1, 12, int),
            ]

    def AddRecord(self, record):
        self.records.append(record)
        #print record.caseid

class Pregnancies(Table):
    def __init__(self, filename='2002FemPreg.dat.gz'):
        self.records = []
        self.ReadFile(filename, self.GetAttrs(), Pregnancy)

    def GetAttrs(self):
        return [
            ('caseid', 1, 12, int),
            ('prglength', 275, 276, int),
            ('outcome', 277, 277, int),
            ('birthord', 278, 279, int),
            ('finalwgt', 423, 440, float),
            ]

    def AddRecord(self, record):
        self.records.append(record)
        #print record.outcome, record.prglength, 
        #print record.birthord, record.finalwgt

def main(name):
    resp = Respondents()
    print len(resp.records)
    preg = Pregnancies()
    print len(preg.records)
    
if __name__ == '__main__':
    main(*sys.argv)
