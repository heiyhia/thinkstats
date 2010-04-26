"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib.pyplot as plt
import Pmf


def main():
    t = [1, 2, 2, 3, 5]
    hist = Pmf.MakeHist(t)

    xs, fs = hist.Render()
    list_of_rectangles = plt.bar(xs, fs)

    plt.xlabel('x')
    plt.ylabel('frequency')
    plt.title('Histogram')

    filename = 'figure'
    plt.savefig(filename, format='png')

if __name__ == "__main__":
    main()

