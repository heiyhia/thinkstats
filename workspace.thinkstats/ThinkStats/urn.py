"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import myplot

from thinkbayes import Suite


class Urn(Suite):
    def Likelihood(self, data, hypo):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer number of blue marbles (out of three)
        data: string 'B' or 'G'
        """
        p = hypo / 3.0
        if data == 'B':
            return p
        else:
            return 1-p
  

def main():
    suite = Urn([0,1,2,3])

    for data in 'BBBBB':
        suite.Update(data)

    suite.Print()
    print suite.Mean() / 3.0

    for data in 'B'*12 + 'G'*2:
        suite.Update(data)

    suite.Print()
    print suite.Mean() / 3.0


if __name__ == '__main__':
    main()
