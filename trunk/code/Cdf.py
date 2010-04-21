# Copyright 2010 Allen B. Downey
#
# License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html

"""Functions for building CDFs (cumulative distribution functions)."""

import bisect
import math

Log10 = math.log10
Pow10 = lambda x: math.pow(10.0, x)

class Cdf(object):
    """Represents a cumulative distribution function.

    Attributes:
        vs: sequence of values
        ps: sequence of probabilities
        name: string used as a graph label.
    """
    def __init__(self, vs=None, ps=None, name=''):
        self.vs = vs
        self.ps = ps
        self.name = name

    def Prob(self, v):
        """Returns the probability that corresponds to value v.

        Computes CDF(v)

        Args:
            v: number
        Returns:
            float probability
        """
        if v < self.vs[0]: return 0.0
        index = bisect.bisect(self.vs, v)
        p = self.ps[index-1]
        return p

    def Value(self, p):
        """Returns the value that corresponds to probability p.

        Computes InverseCDF(p)

        Args:
            p: number in the range [0, 1]

        Returns:
            number value
        """
        if p < 0 or p > 1:
            raise ValueError('Probability p must be in range [0, 1]')

        if p == 0: return self.vs[0]
        if p == 1: return self.vs[-1]
        index = bisect.bisect(self.ps, p)
        if p == self.ps[index-1]:
            return self.vs[index-1]
        else:
            return self.vs[index]

    def Percentile(self, p):
        """Returns the value that corresponds to percentile p.

        Args:
            p: number in the range [0, 100]

        Returns:
            number value
        """
        return self.Value(p / 100.0)

    def Mean(self):
        """Computes the mean of a CDF.

        Returns:
            float mean
        """
        old_p = 0
        total = 0.0
        for v, new_p in zip(self.vs, self.ps):
            p = new_p - old_p
            total += p * v
            old_p = new_p
        return total

    def Round(self, multiplier=1000.0):
        """
        An entry is added to the cdf only if the percentile differs
        from the previous value in a significant digit, where the number
        of significant digits is determined by multiplier.  The
        default is 1000, which keeps log10(1000) = 3 significant digits.
        """
        # TODO(write this method)
        pass

    def Render(self):
        """Generates a sequence of points suitable for plotting.

        An empirical CDF is a step function; linear interpolation
        can be misleading.

        Returns:
            sequence of (value, prob) pairs.
        """
        data = [(self.vs[0], 0.0)]
        for i, p in enumerate(self.ps):
            data.append((self.vs[i], p))
            try:
                data.append((self.vs[i+1], p))
            except IndexError:
                pass
        return data

def MakeCdf(items, name=''):
    """Makes a cdf from an unsorted histogram.

    Args:
        items: unsorted sequence of (value, frequency) pairs
        name: string name for this CDF

    Returns:
        cdf: list of (value, fraction) pairs
    """
    runsum = 0
    vs = []
    cs = []

    for value, count in sorted(items):
        runsum += count
        vs.append(value)
        cs.append(runsum)

    total = float(runsum)
    ps = [c/total for c in cs]

    cdf = Cdf(vs, ps, name)
    return cdf

def MakeCdfFromDict(d, name=''):
    return MakeCdf(d.iteritems(), name)

def MakeCdfFromList(seq, name=''):
  """Creates a CDF from an unsorted sequence.
  Args:
    seq: unsorted sequence of sortable values
    name: string name for the cdf
  Returns:
    CDF object
  """
  hist = {}
  for x in seq:
    hist[x] = hist.get(x, 0) + 1
  return MakeCdfFromDict(hist, name)

def ProbLess(cdf1, cdf2):
  """Probability that a value from cdf1 is less than one from cdf2.
  For continuous distributions F and G, the chance that a sample
  from F is less than a sample from G is
  Integral(x): pdf_F(x) * (1 - cdf_G(x))
  This function computes an approximation of this Integral using
  discrete CDFs.
  Args:
    cdf1: CDF object
    cdf2: CDF object
  Returns:
    float probability
  """
  total = 0.0
  i = 0
  j = 0
  x = float('-Inf')
  while True:
    # sweep through cdf1 and compute p, the marginal prob of v2
    unused_v1, p1 = cdf1.data[i]
    v2, p2 = cdf1.data[i+1]
    p = p2 - p1
    # incr through cdf2 to find Prob{x < v2}
    while x <= v2:
      x, y = cdf2.data[j]
      if j == len(cdf2.data)-1:
        break
      else:
        j += 1
    # add up the integral
    total += p * (1 - y)
    i += 1
    if i == len(cdf1.data)-1:
      break
  return total

