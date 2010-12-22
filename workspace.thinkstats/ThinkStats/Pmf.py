"""This file contains code for use with "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

"""Functions for building PMFs (probability mass functions)."""

class _DictWrapper(object):
    """An object that contains a dictionary."""

    def __init__(self, d=None, name=''):
        if d == None:
            d = {}
        self.d = d
        self.name = name

    def GetDict(self):
        """Gets the dictionary."""
        return self.d

    def Values(self):
        """Gets an iterator of the values.

        Note: one source of confusion is that the keys in this
        dictionaries are the values of the Hist/Pmf, and the
        values are frequencies/probabilities.
        """
        return self.d.iterkeys()

    def Items(self):
        """Gets an iterator of the items (value, freq/prob) pairs."""
        return self.d.iteritems()

    def Render(self):
        """Generates a sequence of points suitable for plotting.

        Returns:
            tuple of (sorted value sequence, freq/prob sequence)
        """
        return zip(*sorted(self.Items()))


class Hist(_DictWrapper):
    def Incr(self, x, count=1):
        """Increments the counter associated with the value x.

        Args:
            x: number value
            count: how much to increment by
        """
        self.d[x] = self.d.get(x, 0) + count

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
    

class Pmf(_DictWrapper):
    """Represents a probability mass function.

    Attributes:
        map: dictionary that maps from values to probabilities
    """
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


def MakeHistFromList(t, name=''):
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


def MakeHistFromDict(d, name=''):
    """Makes a histogram from a map from values to frequencies.

    Args:
        d: dictionary that maps values to frequencies
        name: string name for this histogram

    Returns:
        Hist object
    """
    return Hist(d, name)


def MakePmfFromList(t, name=''):
    """Makes a PMF from an unsorted sequence of values.

    Args:
        t: sequence of numbers
        name: string name for this PMF

    Returns:
        Pmf object
    """
    hist = MakeHistFromList(t, name)
    return MakePmfFromHist(hist)


def MakePmfFromDict(d, name=''):
    """Makes a PMF from a map from values to probabilities.

    Args:
        d: dictionary that maps values to probabilities
        name: string name for this PMF

    Returns:
        Pmf object
    """
    return Pmf(d, name)


def MakePmfFromHist(hist, name=None):
    """Makes a PMF from a Hist object.

    Args:
        hist: Hist object
        name: string name

    Returns:
        Pmf object
    """
    if name is None:
        name = hist.name

    # make a copy of the dictionary
    d = dict(hist.GetDict())
    pmf = Pmf(d, hist.name)
    pmf.Normalize()
    return pmf


def MakePmfFromCdf(cdf, name=None):
    """Makes a PMF from a Cdf object.

    Args:
        cdf: Cdf object

    Returns:
        Pmf object
    """
    if name is None:
        name = cdf.name

    pmf = Pmf()

    prev = 0.0
    for val, prob in cdf.Items():
        pmf.Incr(val, prov-prev)
        prev = prob

    return pmf
