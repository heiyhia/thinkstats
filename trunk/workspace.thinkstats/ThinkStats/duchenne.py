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


def ReadLabels(filename="duchenne/mturklabels.txt"):
    """Returns a list of (photo, labeler, label) tuples.
    """
    labels = []
    for line in open(filename):
        #photo, labeler, label = line.split()
        labels.append(line.split())
    return labels


def MakeObjects(labels):
    """Make Photo and Labeler objects.

    Return: (photos, labelers), a map from photo code to Photo
            and a map from labeler code to Labeler
    """
    photos = {}
    labelers = {}
    for pcode, lcode, label in labels:
        photos[pcode] = Photo()
        labelers[lcode] = Labeler()

    return photos, labelers


class Labeler(thinkbayes.Suite):
    """Represents hypotheses about the trustworthiness of a labeler."""

    def __init__(self):
        thinkbayes.Suite.__init__(self)
        beta = thinkbayes.Beta(2, 1)
        for val, prob in beta.MakePmf().Items():
            self.Set(val, prob)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer value of x, the prob of a correct vote (0-100)
        data: (vote, q) pair, where vote is '1' or '0' and
              q is the mean quality of the link
        """
        x = hypo
        vote, q = data

        if vote == '1':
            return x * q + (1-x) * (1-q)
        elif vote == '0':
            return x * (1-q) + (1-x) * q
        else:
            raise ValueError


class Photo(thinkbayes.Suite):
    """Represents hypotheses about the trustworthiness of a labeler."""

    def __init__(self):
        thinkbayes.Suite.__init__(self)
        beta = thinkbayes.Beta(1, 1)
        for val, prob in beta.MakePmf().Items():
            self.Set(val, prob)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under the hypothesis.

        hypo: integer value of x, the prob of garnering an upvote
        data: (vote, t) pair, where vote is '1' or '0' and
              t is the mean trustworthiness of the labeler
        """
        x = hypo
        vote, t = data

        if vote == '1':
            return x * t + (1-x) * (1-t)
        elif vote == '0':
            return x * (1-t) + (1-x) * t
        else:
            raise ValueError


def Summarize(suite):
    """Prints summary statistics for the suite."""
    print 'MLE', suite.MaximumLikelihood()

    print 'Mean', suite.Mean()
    print 'Median', thinkbayes.Percentile(suite, 50) 

    print '5th %ile', thinkbayes.Percentile(suite, 5) 
    print '95th %ile', thinkbayes.Percentile(suite, 95) 

    print 'CI', thinkbayes.CredibleInterval(suite, 90)


def RunUpdates(photos, labelers, labels):
    """Runs the mutual update process.

    photos: map from id to Photo
    labelers: map from id to Labeler
    labels: label applied to the Photo by the Labeler
    """
    for pcode, lcode, label in labels:
        photo = photos[pcode]
        labeler = labelers[lcode]
        print label, pcode, lcode
        Update(photo, labeler, label)


def Update(photo, labeler, label):
    """Updates the photo and labeler.

    photo: Photo
    labeler: Labeler
    label: string '1' or '0'
    """
    mean_t = labeler.Mean()
    mean_q = photo.Mean()

    print 'trustworthiness', mean_t
    print 'quality', mean_q

    # perform simultaneous updates
    labeler.Update((label, mean_q))
    photo.Update((label, mean_t))


def PlotPosteriorMeans(d, name):
    """Plots the CDF of the means of the posteriors.

    d: map from code to posterior Suite
    name: label for the cdf
    """
    means = [item.Mean() for item in d.itervalues()]
    cdf = thinkbayes.MakeCdfFromList(means, name=name)
    myplot.Cdf(cdf)


def Main():
    labels = ReadLabels()
    photos, labelers = MakeObjects(labels)
    for code, photo in photos.iteritems():
        print code
    for code, labeler in labelers.iteritems():
        print code

    RunUpdates(photos, labelers, labels)

    myplot.Clf()
    PlotPosteriorMeans(photos, 'photos')
    PlotPosteriorMeans(labelers, 'labelers')
    myplot.Show()


if __name__ == '__main__':
    Main()
