"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import matplotlib.pyplot as plt
import Pmf

class Test(unittest.TestCase):

    def testPlot(self):
        t = [1, 2, 2, 3, 5]
        hist = Pmf.MakeHist(t)
        xs, fs = hist.Render()
        list_of_rectangles = plt.bar(xs, fs)

        self.assertEquals(len(list_of_rectangles), 4)

        plt.xlabel('x')
        plt.ylabel('frequency')
        plt.title('Histogram')

        format = 'png'
        filename = 'figure.' + format
        plt.savefig(filename, format=format)


if __name__ == "__main__":
    unittest.main()
