"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import sys
import survey

class Pregnancies(survey.Pregnancies):

    def PartitionRecords(self, constructor):
        """Divides records into two lists: first babies and others.

        Only live births are included

        Args:
            records: sequence of Pregnancy records

            constructor: init method used to make sub-tables
        """
        firsts = constructor()
        others = constructor()

        for p in self.records:
            if p.outcome != 1:
                continue

            if p.birthord == 1:
                firsts.AddRecord(p)
            else:
                others.AddRecord(p)

        return firsts, others

    def Process(self):
        self.lengths = [p.prglength for p in self.records]
        self.mu = Mean(self.lengths)


def PoolRecords(constructor, *tables):
    """Construct a table with records from all tables.
    
    Args:
        constructor: init method used to make the new table
    
        tables: any number of tables

    Returns:
        new table object
    """
    pool = constructor()
    for table in tables:
        pool.ExtendRecords(table.records)
    return pool


def Mean(t):
    """Computes the mean of a sequence of numbers.

    Args:
        t: sequence of numbers

    Returns:
        float
    """
    return float(sum(t)) / len(t)


if __name__ == '__main__':
    main(*sys.argv)
