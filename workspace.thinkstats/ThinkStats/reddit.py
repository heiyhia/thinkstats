"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

"""
This is a modified version of liar.py

This file contains a solution to a problem posed on reddit:

http://www.reddit.com/r/statistics/comments/15rurz/
question_about_continuous_bayesian_inference/

I came across this while working on a persona project. Though the
actual problem has to do with estimating article quality on a
reddit-like site I think it's best way to explain it like this:

I have a biased coin that comes up heads X% of the time. However I
don't know the value X. I know that if I flip it myself I can estimate
the probability of any given bias X by using posterior probability
density function.

However, let's say that I can't look at the coin when it lands. I have
to have an untrustworthy friend check for me. I know my friend will
lie to me Y% of the time. Now obviously if he lies exactly 50% of the
time I'll never be able to determine anything about the coin. But if 
Y!=50% then I can.

Given a prior probably density for the coin Pr(X) and a probability
density for the trustworthiness of my friend Pr(Y) how can I determine
the new probability density for the coin bias after my friend has
reported that he sees it land heads?


And a followup message:

The motivation for this is a voting system I'm trying to develop for
sites like reddit. The perennial complaint is that when a subreddit
first starts it has great content but then newbs move in and upvote
stupid memes and stuff. This new voting system would seek to
ameliorate this effect by weighting different voters
differently. Newbs would have votes with zero weight until they've
proven themselves trustworthy. We judge a newbs trustworthiness by how
well their votes predict the votes of founders.

The model is this:

* Each article has a quality Q
* Each user has a trustworthiness T
* A group of users are "founders" and have trustworthiness T=1
* When a founder votes on an article the prob of an upvote is Q
* When a non-founder votes on something they report the wrong vote 
  T percent of the time

Then we have to come up with:

* A system for rating a user's trustworthiness
* A system for judging article quality given the upvote of a user with
  trustworthiness T.

"""

import thinkbayes
import myplot


class Redditor(thinkbayes.Suite):
    """Represents hypotheses about the trustworthiness of a redditor."""

    def Likelihood(self, data, hypo):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer value of x, the prob of a correct vote (0-100)
        data: (vote, q) pair, where vote is 'up' or 'down' and
              q is the mean quality of the link
        """
        x = hypo / 100.0
        vote, q = data

        if vote == 'up':
            return x * q + (1-x) * (1-q)
        elif vote == 'down':
            return x * (1-q) + (1-x) * q
        else:
            return 0

class Item(thinkbayes.Suite):
    """Represents hypotheses about the trustworthiness of a redditor."""

    def Likelihood(self, data, hypo):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer value of x, the prob of garnering an upvote
        data: (vote, t) pair, where vote is 'up' or 'down' and
              t is the mean trustworthiness of the redditor
        """
        x = hypo / 100.0
        vote, t = data

        if vote == 'up':
            return x * t + (1-x) * (1-t)
        elif vote == 'down':
            return x * (1-t) + (1-x) * t
        else:
            return 0


def Summarize(suite):
    """Prints summary statistics for the suite."""
    print 'MLE', suite.MaximumLikelihood()

    print 'Mean', suite.Mean()
    print 'Median', thinkbayes.Percentile(suite, 50) 

    print '5th %ile', thinkbayes.Percentile(suite, 5) 
    print '95th %ile', thinkbayes.Percentile(suite, 95) 

    print 'CI', thinkbayes.CredibleInterval(suite, 90)


def Main():
    # make a redditor with some trustworthiness (mean_t = 0.67)
    founder = Redditor(name='redditor')
    beta = thinkbayes.Beta(2, 1)
    for val, prob in beta.MakePmf().Items():
        founder.Set(val*100, prob)

    # make a new item with unknown quality (mean_q = 0.5)
    item = Item(range(0, 101), name='item')

    # compute the means
    mean_t = founder.Mean() / 100.0
    mean_q = item.Mean() / 100.0

    print mean_t
    print mean_q

    # perform simultaneous updates
    founder.Update(('up', mean_q))
    item.Update(('up', mean_t))

    Summarize(item)

    # display the posterior distributions
    myplot.Pmf(founder)
    myplot.Pmf(item)
    myplot.Show()


if __name__ == '__main__':
    Main()
