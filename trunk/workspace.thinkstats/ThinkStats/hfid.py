import scipy.stats

from math import sqrt

class Gaussian:
    """Represents a Gaussian distribution."""

    def __init__(self, mu, sigma2):
        """Constructs a Gaussian distribution with given mu and sigma.

        mu: mean
        sigma2: variance (square of stdev)
        """
        self.mu = mu
        self.sigma2 = sigma2

    def __str__(self):
        return 'N(%f, %f)' % (self.mu, self.sigma2) 


def SubtractGaussian(g1, g2):
    """Difference between two Gaussians.

    returns: Gaussian object
    """
    return Gaussian(g1.mu - g2.mu, g1.sigma2 + g2.sigma2)


def Pooled(g1, n1, g2, n2):
    """Pooled standard deviation.

    returns: standard deviation
    """
    Sp2 = (n1 * g1.sigma2 + n2 * g2.sigma2) / (n1 + n2)
    return sqrt(Sp2)


def main():
    # fill in these values
    m1 = Gaussian(4.0, 0.5)
    f1 = Gaussian(3.0, 0.5)
    n1 = 30.0

    m2 = Gaussian(3.75, 0.5)
    f2 = Gaussian(3.25, 0.5)
    n2 = 90.0

    # don't need to change anything else below here...
    d1 = SubtractGaussian(m1, f1)
    d2 = SubtractGaussian(m2, f2)

    print 'd1 =', d1
    print 'd2 =', d2

    delta = SubtractGaussian(d2, d1).mu
    print 'delta', delta

    Sp = Pooled(d1, n1, d2, n2)

    t = delta / Sp / sqrt(1/n1 + 1/n2)
    print 't-stat =', t

    df = n1 + n2
    p = scipy.stats.t.cdf(t, df)
    print 'p-value =', p


if __name__ == '__main__':
    main()

