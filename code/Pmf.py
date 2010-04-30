# Copyright 2010 Allen B. Downey
#
# License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

"""Functions for building PMFs (probability mass functions)."""

class Hist(object):
    def __init__(self, d=None, name=''):
        if d == None:
            d = {}
        self.d = d
        self.name = name

    def GetDict(self):
        return self.d

    def Count(self, x):
        """Increment the counter associated with the value x.

        Args:
            x: number value
        """
        self.d[x] = self.d.get(x, 0) + 1

    def Freq(self, x):
        """Get the frequency associated with the value x.

        Args:
            x: number value

        Returns:
            int frequency
        """
        return self.d.get(x, 0)

    def Render(self):
        """Generates a sequence of points suitable for plotting.

        Returns:
            sequence of (value, prob) pairs.
        """
        items = self.d.items()
        items.sort()
        return zip(*items)


class Pmf(object):
    """Represents a probability mass function.

    Attributes:
        map: dictionary that maps from values to probabilities
    """
    def __init__(self, d=None, name=''):
        if d == None:
            d = {}
        self.d = d
        self.name = name

    def Normalize(self, denom=None):
        """Normalizes this PMF so the sum of all probs is 1.

        Args:
            denom: float divisor; if None, computes the total of all probs
        """
        denom = denom or self.Total()
        
        for x, p in self.d.iteritems():
            self.d[x] = p / denom
    
    def Total(self):
        """Returns the total of the frequencies in the map."""
        total = sum(self.d.values())
        return float(total)

    def GetDict(self):
        return self.d

    def Prob(self, x):
        """Returns the probability that corresponds to value x.

        Args:
            x: number

        Returns:
            int frequency or float probability
        """
        p = self.d.get(x, 0)
        return p

    def Mean(self):
        """Computes the mean of a PMF.

        Returns:
            float mean
        """

    def Render(self):
        """Generates a sequence of points suitable for plotting.

        Returns:
            sequence of (value, prob) pairs.
        """
        items = self.d.items()
        items.sort()
        return zip(*items)


def MakeHist(t, name=''):
    """Makes a histogram from an unsorted sequence of values.

    Args:
        t: sequence of numbers

        name: string name for this histogram

    Returns:
        Hist object
    """
    hist = Hist(name=name)
    [hist.Count(x) for x in t]
    return hist


def MakePmf(t, name=''):
    """Makes a PMF from an unsorted sequence of values.

    Args:
        t: sequence of numbers

        name: string name for this PMF

    Returns:
        Pmf object
    """
    hist = MakeHist(t, name)
    return MakePmfFromHist(hist)

def MakePmfFromHist(hist):
    """Makes a PMF from a Hist object.

    Args:
        hist: Hist object

    Returns:
        Pmf object
    """
    d = dict(hist.GetDict())
    pmf = Pmf(d, hist.name)
    pmf.Normalize()
    return pmf
