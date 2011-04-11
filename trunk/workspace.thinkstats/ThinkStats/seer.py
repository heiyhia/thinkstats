"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import sys
import gzip
import os

import Pmf
import table

class Records(table.Table):
    """Represents the record table."""

    def ReadRecords(self, data_dir=None, n=None):
        if data_dir == None:
            data_dir = './SEER/incidence'
        filename = self.GetFilename()
        self.ReadFile(data_dir, filename, self.GetFields(), table.Record, n)
        self.Recode()

    def GetFilename(self):
        return 'yr1973_2007.seer9/BREAST.TXT.gz'

    def GetFields(self):
        """Returns a tuple specifying the fields to extract.

        The elements of the tuple are field, start, end, case.

                field is the name of the variable
                start and end are the indices as specified in the NSFG docs
                cast is a callable that converts the result to int, float, etc.

        Codes:
        sex: 1=male, 2=female
        age: years or 999
        seqno: 0 = only one primary in lifetime
        diagmo: 1-12
        diagyr: 1973+
        site:
        follow: 2 or 4 = active followup
        behavior: 3 = malignant
        race:
        hisp:
        stage: 
        survival: YYMM or 9999
        cause:
        status: 1=alive, 4=dead
        deathclass: 0=considered "alive"  1=death due to primary tumor

        The Survival Time Recode is calculated using the date of
        diagnosis and one of the following: date of death, date last known to
        be alive, or follow-up cutoff date used for this file (see title page
        for date for this file). Thus a person diagnosed in May 1976 and who
        died in May 1980 has a Survival Time Recode of 04 years and 00 months.

        """
        return [
            ('caseid', 1, 8, int),
            ('sex', 24, 24, int),
            ('age', 25, 27, int),
            ('seqno', 35, 36, int),
            ('diagmo', 37, 38, int),
            ('diagyr', 39, 42, int),
            ('site', 43, 46, int),
            ('follow', 181, 181, int),
            ('behavior', 215, 215, int),
            ('race', 224, 224, int),
            ('hisp', 225, 225, int),
            ('stage', 229, 230, int),
            ('survival', 241, 244, int),
            ('cause', 245, 249, int),
            ('status', 255, 255, int),
            ('deathclass', 264, 264, int),
            ]

    def MakeHists(self):
        for field in self.GetFields():
            attr = field[0]
            if attr not in ['caseid', 'cause', 'survival']:
                self.MakeHist(attr)

    def MakeHist(self, attr):

        vals = [getattr(record, attr) for record in self.records]
        hist = Pmf.MakeHistFromList(vals)

        print attr
        for val, freq in hist.Items():
            print val, freq
        print


def main(name, data_dir=None):
    table = Records()
    table.ReadRecords(data_dir=data_dir, n=1000000)
    print 'Number of records', len(table.records)

    table.MakeHists()

    
if __name__ == '__main__':
    main(*sys.argv)
