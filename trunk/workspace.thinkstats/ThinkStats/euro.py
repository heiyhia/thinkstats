"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

"""This file contains a partial solution to a problem from
MacKay, "Information Theory, Inference, and Learning Algorithms."

    Exercise 3.15 (page 50): A statistical statement appeared in
    "The Guardian" on Friday January 4, 2002:

        When spun on edge 250 times, a Belgian one-euro coin came
        up heads 140 times and tails 110.  'It looks very suspicious
        to me,' said Barry Blight, a statistics lecturer at the London
        School of Economics.  'If the coin were unbiased, the chance of
        getting a result as extreme as that would be less than 7%.'

MacKay asks, "But do these data give evidence that the coin is biased
rather than fair?"

"""

from thinkbayes import Suite
import myplot


class Euro(Suite):

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer value of x, the probability of heads (0-100)
        data: string 'H' or 'T'
        """
        x = hypo / 100.0
        if data == 'H':
            return x
        else:
            return 1-x


def RunUniformPrior():
    suite = Euro(xrange(0, 101))
    dataset = 'H' * 140 + 'T' * 110

    for data in dataset:
        suite.Update(data)

    myplot.Clf()
    myplot.Pmf(suite)
    myplot.Save(root='euro1',
                xlabel='x',
                ylabel='Probability')


def Main():
    RunUniformPrior()


def MLE(pmf):
    prob, val = max((prob, val) for val, prob in pmf.Items())
    return val


def TrianglePrior():
    pmf = Pmf.Pmf()
    for x in range(0, 51):
        pmf.Set(x, x)
    for x in range(51, 101):
        pmf.Set(x, 100-x) 
    pmf.Normalize()
    return pmf


def RunUpdate(pmf):
    evidence = 'H' * 140 + 'T' * 110

    for outcome in evidence:
        Update(pmf, outcome)

    print pmf.Prob(50)
    # 0.0209765261295

    print 'MLE', MLE(pmf)

    print 'Mean', Mean(pmf)

    print '5th %ile', Percentile(pmf, 0.05) 
    print '95th %ile', Percentile(pmf, 0.95) 


def Ratio():
    biased = Pmf.Pmf()
    for x in range(0, 101):
        biased.Set(x, 1)
    biased.Set(50, 0)
    biased.Normalize()

    print Likelihood(biased, 140, 110)

    fair = Pmf.Pmf()
    fair.Set(50, 1)
    print Likelihood(fair, 140, 110)


def Likelihood(pmf, heads, tails):
    total = 0
    for x in pmf.Values():
        p = x / 100.0
        likelihood = p**heads * (1-p)**tails
        total += pmf.Prob(x) * likelihood
    return total


def Leftovers():
    Ratio()
    return

    pmf1 = Pmf.Pmf()
    for x in range(0, 101):
        pmf1.Set(x, 1)
    pmf1.Normalize()

    pmf2 = TrianglePrior()

    # plot the priors
    myplot.Clf()
    myplot.Pmfs([pmf1, pmf2])
    myplot.Save(root='simple_coin_both_prior',
                title='Biased coin',
                xlabel='x',
                ylabel='Probability')

    RunUpdate(pmf1)
    RunUpdate(pmf2)

    # plot the posterior distributions
    myplot.Clf()
    myplot.Pmfs([pmf1, pmf2])
    myplot.Save(root='simple_coin_both_post',
               title='Biased coin',
               xlabel='x',
               ylabel='Probability')

if __name__ == '__main__':
    Main()
