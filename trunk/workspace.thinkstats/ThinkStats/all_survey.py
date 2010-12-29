"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import survey

class Respondents1995(Table):
    """Represents the respondent table."""

    def GetFilename(self):
        return '1995FemResp.dat.gz'

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

class Pregnancies1995(Table):
    """Contains survey data about a Pregnancy."""

    def GetFilename(self):
        return '1995FemPreg.dat.gz'

    def GetFields(self):
        """Gets information about the fields to extract from the survey data.

        Documentation of the fields for Cycle 5is at
        ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/
        NSFG/Cycle5Codebook-Pregnancy.pdf

        Returns:
            sequence of (name, start, end, type) tuples
        """
        return [
            ('caseid', 1, 8, int),
            ('nbrnaliv', 14, 14, int),
            ('babysex', 325, 325, int),
            ('outcome', 11, 11, int),
            ('birthord', 9, 10, int),
            ('finalwgt', 403, 412, float),
            ]

    def GetLiveBirthCodes(self):
        return [1]

def main(name, data_dir='.'):
    resp = Respondents()
    resp.ReadRecords(data_dir)
    print 'Number of respondents', len(resp.records)

    preg = Pregnancies()
    preg.ReadRecords(data_dir)
    print 'Number of pregnancies', len(preg.records)
    
if __name__ == '__main__':
    main(*sys.argv)
