"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import sys

import survey

class Respondents1995(survey.Respondents):
    """Represents the respondent table."""

    def GetFilename(self):
        return '1995FemRespData.dat.gz'

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.
        """
        return [
            ('caseid', 1, 8, int),
            ]

class Pregnancies1995(survey.Pregnancies):
    """Contains survey data about a Pregnancy."""

    def GetFilename(self):
        return '1995PregData.dat.gz'

    def GetFields(self):
        """Gets information about the fields to extract from the survey data.

        Documentation of the fields for Cycle 5 is at
        ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/
        NSFG/Cycle5Codebook-Pregnancy.pdf

        Returns:
            sequence of (name, start, end, type) tuples
        """
        return [
            ('caseid', 1, 8, int),
            ('nbrnaliv', 14, 14, int),    # nbrnlv
            ('babysex', 307, 307, int),   # sex1
            ('outcome', 273, 273, int),     # outcome
            ('birthord', 9, 10, int),     # preordr
            ('finalwgt', 385, 394, float), # post_wt
            ]


class Pregnancies1988(survey.Pregnancies):
    """Contains survey data about a Pregnancy."""

    def GetFilename(self):
        return '1988PregData.dat.gz'

    def GetFields(self):
        """Gets information about the fields to extract from the survey data.

        Documentation of the fields for Cycle 4 is at
        ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/
        NSFG/Cycle4Codebook.pdf

        Returns:
            sequence of (name, start, end, type) tuples
        """
        return [
            ('caseid', 1, 5, int),
            ('nbrnaliv', 12, 12, int),     # numbirth
            ('babysex', 75, 75, int),      # b12_1
            ('outcome', 284, 284, int),    # outcome
            ('birthord', 8, 9, int),       # pregnum
            ('finalwgt', 338, 344, float), # w_5
            ]


class Respondents1982(survey.Respondents):
    """Represents the respondent table."""

    def GetFilename(self):
        return '1982FemRespData.dat.gz'

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.
        """
        return [
            ('caseid', 1152, 1156, int),
            ('finalwgt', 756, 761, float),
            ]

class Pregnancies1982(survey.Pregnancies):
    """Contains survey data about a Pregnancy."""

    def GetFilename(self):
        return '1982PregData.dat.gz'

    def GetFields(self):
        """Gets information about the fields to extract from the survey data.

        Documentation of the fields for Cycle 3 is at
        ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Dataset_Documentation/
        NSFG/Cycle3Codebook.pdf

        The SAS file that has the correct column numbers is at
        ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/sas/
        1982PregSetup.sas

        Returns:
            sequence of (name, start, end, type) tuples
        """
        return [
            ('caseid', 276, 280, int),    # id_num
            ('birthord', 4, 5, int),     # pregnancy number
            ('nbrnaliv', 8, 8, int),     # BP1920B
            ('babysex', 26, 26, int),      # B201
            ('outcome', 264, 264, int),    # outcome
            ]


class Respondents1976(survey.Respondents):
    """Represents the respondent table."""

    def GetFilename(self):
        return '1976FemRespData.dat.gz'

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.
        """
        return [
            ('caseid', 1, 5, int),
            ('finalwgt', 721, 726, float),
            ]

class Pregnancies1976(survey.Pregnancies):
    """Contains survey data about a Pregnancy."""

    def GetFilename(self):
        return '1976PregData.dat.gz'

    def GetFields(self):
        """Gets information about the fields to extract from the survey data.

        Documentation of the fields for Cycle 2 is at

        The SAS file that has the correct column numbers is at
        ftp://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/sas/
        1976PregSetup.sas

        Returns:
            sequence of (name, start, end, type) tuples
        """
        return [
            ('caseid', 1, 5, int),    # id
            ('birthord', 7, 8, int),     # BPREC1
            ('babysex', 20, 20, int),      # B5_1
            ('childno_1', 18, 19, int),      # childno_1
            ('pregtype', 14, 14, int),     # BPREC4
            ]

    def Recode(self):
        for record in self.records:
            d = {1:1, 2:4, 3:1, 4:4, 5:1, 6:6}
            record.outcome = d[record.pregtype]

            d = {1:1, 2:0, 3:2, 4:0, 5:1, 6:'NA'}
            record.nbrnaliv = d[record.pregtype]


def main(name, data_dir='.'):
    resp = Respondents1995()
    resp.ReadRecords(data_dir)
    print 'Number of respondents', len(resp.records)

    preg = Pregnancies1995()
    preg.ReadRecords(data_dir)
    print 'Number of pregnancies', len(preg.records)
    
if __name__ == '__main__':
    main(*sys.argv)
