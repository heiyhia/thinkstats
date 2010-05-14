"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import sys
import gzip

class Record(object):
    """Represents a record."""

class Respondent(Record): 
    """Represents a respondent."""

class Pregnancy(Record):
    """Represents a pregnancy."""

class Table(object):
    """Represents a table as a list of objects"""

    def __init__(self):
        self.records = []

    def ReadFile(self, filename, fields, constructor, n=None):
        """Reads a compressed data file builds one object per record.

        Args:
            filename: string name of the file to read

            fields: sequence of (name, start, end, cast) tuples specifying 
            the fields to extract

            constructor: what kind of object to create
        """
        if filename.endswith('gz'):
            fp = gzip.open(filename)
        else:
            fp = open(filename)

        for i, line in enumerate(fp):
            if i == n:
                break
            record = self.MakeRecord(line, fields, constructor)
            self.AddRecord(record)
        fp.close()
    
    def MakeRecord(self, line, fields, constructor):
        """Scans a line and returns an object with the appropriate fields.

        Args:
            line: string line from a data file

            fields: sequence of (name, start, end, cast) tuples specifying 
            the fields to extract

            constructor: callable that makes an object for the record.

        Returns:
            Record with appropriate fields.
        """
        obj = constructor()
        for (field, start, end, cast) in fields:
            try:
                val = cast(line[start-1:end])
            except ValueError:
                val = 'NA'
            setattr(obj, field, val)
        return obj

    def AddRecord(self, record):
        """Adds a record to this table.

        Args:
            record: an object of one of the record types.
        """
        self.records.append(record)

    def ExtendRecords(self, records):
        """Adds records to this table.

        Args:
            records: a sequence of record object
        """
        self.records.extend(records)


class Respondents(Table):
    """Represents the respondent table."""

    def ReadRecords(self, filename='2002FemResp.dat.gz', n=None):
        self.ReadFile(filename, self.GetFields(), Respondent, n)

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.
        """
        return [
            ('caseid', 1, 12, int),
            ]

class Pregnancies(Table):
    """Contains survey data about a Pregnancy."""

    def ReadRecords(self, filename='2002FemPreg.dat.gz', n=None):
        self.ReadFile(filename, self.GetFields(), Pregnancy, n)

    def GetFields(self):
        """Gets information about the fields to extract from the survey data.

        Documentation of the fields is at
        http://nsfg.icpsr.umich.edu/cocoon/WebDocs/NSFG/public/index.htm

        Returns:
            sequence of (name, start, end, type) tuples
        """
        return [
            ('caseid', 1, 12, int),
            ('nbrnaliv', 22, 22, int),
            ('birthwgt_lb', 57, 58, int),
            ('birthwgt_oz', 59, 60, int),
            ('prglength', 275, 276, int),
            ('outcome', 277, 277, int),
            ('birthord', 278, 279, int),
            ('finalwgt', 423, 440, float),
            ]

def main(name):
    resp = Respondents()
    resp.ReadRecords()
    print 'Number of respondents', len(resp.records)

    preg = Pregnancies()
    preg.ReadRecords()
    print 'Number of pregnancies', len(preg.records)
    
if __name__ == '__main__':
    main(*sys.argv)
