"""This file contains code used in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from math import pow


class Bayes(object):

    def Update(self, prior, evidence, verbose=False):
        """Updates a prior based on evidence.

        Args:
            prior: Pmf that maps hypotheses to probabilities
            evidence: whatever kind of object Likelihood expects

        Returns:
            posterior: Pmf that maps hypotheses to probabilities
        """
        posterior = prior.Copy()

        for hypo in prior.Values():
            likelihood = self.Likelihood(evidence, hypo)
            if verbose:
                print hypo, likelihood
            posterior.Mult(hypo, likelihood)
        posterior.Normalize()
        return posterior

    def Likelihood(self, evidence, hypo):
        """Computes the likelihood of the evidence assuming the hypothesis.

        Args:
            evidence: some representation of the evidence
            hypo: some representation of the hypothesis

        Returns:
            probability (or likelihood) of the evidence given the hypothesis
        """
        raise Error('');


class BinomialBayes(Bayes):

    def Likelihood(self, evidence, hypo):
        """Computes the likelihood of the evidence assuming the hypothesis.

        Args:
            evidence: a tuple of (number of heads, number of tails)
            hypo: float probability of heads

        Returns:
            probability of tossing the given number of heads and tails with a
            coin that has p probability of heads
        """
        heads, tails = evidence
        p = hypo
        print p, heads, tails
        return pow(p, heads) * pow(1-p, tails)


def main():
    return

if __name__ == '__main__':
    main()
