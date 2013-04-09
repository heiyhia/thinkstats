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
    """Returns a list of (photo, labeler, label) tuples."""
    return ReadTextFile(filename)


def ReadTruth(filename="duchenne/groundtruth.txt"):
    """Returns a list of (photo, label) tuples."""
    return ReadTextFile(filename)


def ReadTextFile(filename):
    labels = []
    for line in open(filename):
        labels.append(line.split())
    return labels


def MakeObjects(labels):
    """Make Photo and Labeler objects.

    Return: (photo_map, labeler_map), a map from photo code to Photo
            and a map from labeler code to Labeler
    """
    photo = Photo()
    labeler = Labeler()

    photo_map = {}
    labeler_map = {}
    for pcode, lcode, label in labels:
        photo_map[pcode] = photo.Copy()
        labeler_map[lcode] = labeler.Copy()

    return photo_map, labeler_map


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


def RunUpdates(photo_map, labeler_map, labels):
    """Runs the mutual update process.

    photo_map: map from id to Photo
    labeler_map: map from id to Labeler
    labels: label applied to the Photo by the Labeler
    """
    for pcode, lcode, label in labels:
        photo = photo_map[pcode]
        labeler = labeler_map[lcode]
        # print label, pcode, lcode
        Update(photo, labeler, label)


def Update(photo, labeler, label):
    """Updates the photo and labeler.

    photo: Photo
    labeler: Labeler
    label: string '1' or '0'
    """
    mean_t = labeler.Mean()
    mean_q = photo.Mean()

    #print 'trustworthiness', mean_t
    #print 'quality', mean_q

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
    truth = ReadTruth()
    truth_map = {}
    for pcode, label in truth:
        truth_map[pcode] = label

    labels = ReadLabels()
    photo_map, labeler_map = MakeObjects(labels)

    RunUpdates(photo_map, labeler_map, labels)

    yes = []
    no = []
    for pcode, photo in photo_map.iteritems():
        if pcode in truth_map:
            mean = photo.Mean()

            if truth_map[pcode] == '1':
                yes.append(mean)
            else:
                no.append(mean)

    myplot.Clf()
    cdf_yes = thinkbayes.MakeCdfFromList(yes, name='yes')
    cdf_no = thinkbayes.MakeCdfFromList(no, name='no')
    myplot.Cdfs([cdf_yes, cdf_no])
    myplot.Show()

    return

    myplot.Clf()
    PlotPosteriorMeans(photo_map, 'photos')
    PlotPosteriorMeans(labeler_map, 'labelers')
    myplot.Show()


if __name__ == '__main__':
    Main()
