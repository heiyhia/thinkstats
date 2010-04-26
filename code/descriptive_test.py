"""This file contains a solution to an exercise in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import unittest
import matplotlib.pyplot as plt
import Pmf
import descriptive

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

    def testDescriptive(self):
        constructor = descriptive.Pregnancies
        table = constructor()
        table.ReadRecords()

        print
        print 'Number of pregnancies', len(table.records)
        self.assertEquals(len(table.records), 13593)

        firsts, others = table.PartitionRecords(constructor)

        print 'Number of first babies', len(firsts.records)
        print 'Number of others', len(others.records)

        self.assertEquals(len(firsts.records), 4413)
        self.assertEquals(len(others.records), 4735)

        firsts.Process()
        others.Process()

        firsts.MakePmf(name='first babies')
        others.MakePmf(name='others')

        width = 0.4

        xs, fs = firsts.hist.Render()
        xs = Shift(xs, -width)
        plt.bar(xs, fs, width=width, color='0.80', label=firsts.hist.name)

        xs, fs = others.hist.Render()
        plt.bar(xs, fs, width=width, color='0.50', label=others.hist.name)

        plt.axis([20, 50, 0, 2700])
        plt.legend()
        plt.show()

def Shift(xs, shift):
    return [x+shift for x in xs]

if __name__ == "__main__":
    unittest.main()
