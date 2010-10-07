"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

"""Functions for building PMFs (probability mass functions)."""

class Hist(object):
    def __init__(self, d=None, name=''):
        if d == None:
            d = {}
        self.d = d
        self.name = name

    def GetDict(self):
        return self.d

    def Incr(self, x, count=1):
        """Increments the counter associated with the value x.

        Args:
            x: number value
            count: how much to increment by
        """
        self.d[x] = self.d.get(x, 0) + count

    def Values(self):
        """Gets a sorted iterator of the values."""
        return sorted(self.d.iteritems())

    def Freq(self, x):
        """Gets the frequency associated with the value x.

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
        
    def GetDict(self):
        return self.d

    def Incr(self, x, p=1):
        """Increments the probability associated with the value x.

        Args:
            x: number value
            p: how much to increment by
        """
        self.d[x] = self.d.get(x, 0) + p

    def Mult(self, x, factor):
        """Scales the probability associated with the value x.

        Args:
            x: number value
            factor: how much to multiply by
        """
        self.d[x] = self.d.get(x, 0) * factor

    def Values(self):
        """Gets a sorted iterator of the values.

        Note one source of confusion: the values in this PMF are stored
        as the _keys_ of the dictionary.  The probabilities are stored
        as values in the dictionary.
        """
        return sorted(self.d.iterkeys())

    def Prob(self, x):
        """Gets the probability associated with the value x.

        Args:
            x: number value

        Returns:
            float probability
        """
        return self.d.get(x, 0)

    def Copy(self, name=None):
        """Returns a copy of this Pmf.

        Args:
            name: string name for the new Pmf
        """
        if name is None:
            name = self.name
        return Pmf(dict(self.d), name)

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

    def Mean(self):
        """Computes the mean of a PMF.

        Returns:
            float mean
        """
        mu = 0.0
        for x, p in self.d.iteritems():
            mu += p * x
        return mu

    def Var(self, mu=None):
        if mu is None:
            mu = self.Mean()
            
        var = 0.0
        for x, p in self.d.iteritems():
            var += p * (x - mu)**2
        return var

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
    [hist.Incr(x) for x in t]
    return hist


def MakePmfFromList(t, name=''):
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
