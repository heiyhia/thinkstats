"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import thinkbayes
import myplot


class Train(thinkbayes.Suite):
    """Suite of hypotheses about the number of trains."""

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer hypothesis about the number of trains
        data: integer serial number of the observed train
        """
        # TODO: write this method!
        return 1


def main():
    hypos = xrange(100, 1001)
    suite = Train(hypos)

    suite.Update(321)
    print suite.Mean()

    myplot.Pmf(suite)
    myplot.Show()


if __name__ == '__main__':
    main()
