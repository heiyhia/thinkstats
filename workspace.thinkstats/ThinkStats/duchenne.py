"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2013 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

"""
Duchenne experiment data accompanying the paper by Whitehill, et al., 
Whose Vote Should Count More: Optimal Integration of Labels from 
Labelers of Unknown Expertise", NIPS 2009:

Data downloaded from 
http://mplab.ucsd.edu/~jake/DuchenneExperiment/DuchenneExperiment.html
on March 7, 2013
"""

import thinkbayes
import myplot


class Labeler(thinkbayes.Suite):
    """Represents hypotheses about the trustworthiness of a labeler."""

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer value of x, the prob of a correct vote (0-100)
        data: (vote, q) pair, where vote is 'yes' or 'no' and
              q is the mean quality of the link
        """
        x = hypo / 100.0
        vote, q = data

        if vote == 'yes':
            return x * q + (1-x) * (1-q)
        elif vote == 'no':
            return x * (1-q) + (1-x) * q
        else:
            return 0

class Photo(thinkbayes.Suite):
    """Represents hypotheses about the trustworthiness of a labeler."""

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer value of x, the prob of garnering an upvote
        data: (vote, t) pair, where vote is 'yes' or 'no' and
              t is the mean trustworthiness of the labeler
        """
        x = hypo / 100.0
        vote, t = data

        if vote == 'yes':
            return x * t + (1-x) * (1-t)
        elif vote == 'no':
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
    # make a labeler with some trustworthiness (mean_t = 0.67)
    labeler = Labeler(name='labeler')
    beta = thinkbayes.Beta(2, 1)
    for val, prob in beta.MakePmf().Items():
        labeler.Set(val*100, prob)

    # make a new photo with unknown quality (mean_q = 0.5)
    photo = Photo(range(0, 101), name='photo')

    # compute the means
    mean_t = labeler.Mean() / 100.0
    mean_q = photo.Mean() / 100.0

    print mean_t
    print mean_q

    # perform simultaneous updates
    labeler.Update(('yes', mean_q))
    photo.Update(('yes', mean_t))

    Summarize(photo)

    # display the posterior distributions
    myplot.Pmf(labeler)
    myplot.Pmf(photo)
    myplot.Show()


if __name__ == '__main__':
    Main()
