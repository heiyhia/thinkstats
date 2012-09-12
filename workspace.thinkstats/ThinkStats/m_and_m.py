"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

from thinkbayes import Pmf


class M_and_M(Pmf):
    def __init__(self, hypos):
        Pmf.__init__(self)
        for hypo in hypos:
            self.Set(hypo, 1)
        self.Normalize()

    def Update(self, data):
        for hypo in self.Values():
            like = self.Likelihood(hypo, data)
            self.Mult(hypo, like)
        self.Normalize()

    mix94 = dict(brown=30,
                 yellow=20,
                 red=20,
                 green=10,
                 orange=10,
                 tan=10)

    mix96 = dict(blue=24,
                 green=20,
                 orange=16,
                 yellow=14,
                 red=13,
                 brown=13)

    hypo1 = dict(bag1=mix94, bag2=mix96)
    hypo2 = dict(bag1=mix96, bag2=mix94)

    hypotheses = dict(A=hypo1, B=hypo2)

    def Likelihood(self, hypo, data):
        bag, color = data
        bags = self.hypotheses[hypo]
        mix = bags[bag]
        like = mix[color]
        return like


def main():
    hypos = ['A', 'B']
    pmf = M_and_M(hypos)

    dataset = [('bag1', 'yellow'),
              ('bag2', 'green')]

    for data in dataset:
        pmf.Update(data)

    for hypo, prob in pmf.Items():
        print hypo, prob


if __name__ == '__main__':
    main()
