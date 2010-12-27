"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import relay_bias
import Pmf
import myplot

def BiasPmf(pmf, speed, name=None):
    """Returns a new PDF representing speeds observed at a given speed.

    The chance of observing a runner is proportional to the difference
    in speed.

    Args:
        pmf: distribution of actual speeds
        speed: speed of the observing runner
        name: string name for the new dist

    Returns:
        Pmf object
    """
    new = pmf.Copy(name=name)
    for val, prob in new.Items():
        diff = abs(val - speed)
        new.Mult(val, diff)
    new.Normalize()
    return new


def main():
    results = relay_bias.ReadResults()
    speeds = relay_bias.GetSpeeds(results)

    pmf = Pmf.MakePmfFromList(speeds, 'actual speeds')
    myplot.Pmf(pmf,
               root='actual_speeds',
               title='PMF of running speed',
               xlabel='speed (mph)',
               ylabel='probability')

    biased = BiasPmf(pmf, 7.5, name='observed speeds')
    myplot.Pmf(biased, 
               root='observed_speeds',
               title='PMF of running speed',
               xlabel='speed (mph)',
               ylabel='probability',
               show=True)


if __name__ == '__main__':
    main()
