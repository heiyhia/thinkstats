"""This file contains code used in "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""

import matplotlib.pyplot as pyplot
import myplot
import numpy

import csv
import HTML
import math
import random
import shelve
import sys
import time

import thinkbayes


formats = ['pdf', 'eps', 'png']


class Locker(object):
    def __init__(self, shelf_file):
        self.shelf = shelve.open(shelf_file)

    def Close(self):
        self.shelf.close()

    def Add(self, key, value):
        self.shelf[key] = value

    def Lookup(self, key):
        return self.shelf.get(key)

    def Keys(self):
        return self.shelf.iterkeys()

    def Read(self):
        return dict(self.shelf)


class Subject(object):
    """Represents a subject from the belly button study."""

    def __init__(self, code):
        self.code = code
        self.species = []

    def Add(self, species, count):
        """Add a species-count pair.

        species: string species/genus name
        count: int number of individuals
        """
        self.species.append((count, species))

    def GetM(self):
        """Gets number of observed species."""
        return len(self.species)
        
    def Sort(self):
        """Sorts the count-species pairs in increasing order."""
        self.species.sort()        

    def GetCounts(self):
        """Gets the list of species counts

        Should be in increasing order, if Sort() has been invoked.
        """
        return [count for count, _ in self.species]

    def GetNames(self):
        """Gets the names of the seen species."""
        return [name for _, name in self.species]

    def PrintCounts(self):
        """Prints the counts and species names."""
        for count, name in reversed(self.species):
            print count, name

    def GetSpecies(self, index):
        """Gets the count and name of the indicated species.

        Returns: count-species pair
        """
        return self.species[index]

    def GetCdf(self):
        counts = self.GetCounts()
        items = enumerate(counts)
        cdf = thinkbayes.MakeCdfFromItems(items)
        return cdf

    def Process(self, iterations=300):
        """Computes the posterior distribution of n and the prevalences.

        Sets: self.suite
        """
        #self.PrintCounts()

        counts = self.GetCounts()
        m = len(counts)

        n = 400
        ns = range(m, n)
        self.suite = Species5(ns, iterations=iterations)

        start = time.time()    
        self.suite.Update(counts)
        end = time.time()
        print 'time', end-start

    def MakeFigures(self):
        """Makes figures showing distribution of n and the prevalences."""
        self.PlotDistOfN()
        self.PlotPrevalences()

    def PlotDistOfN(self):
        """Plots distribution of N."""
        pmf = self.suite.DistOfN()
        print '90% CI for N:', pmf.CredibleInterval(90)
        pmf.name = self.code

        myplot.Clf()
        myplot.PrePlot(num=1)

        myplot.Pmf(pmf)

        root = 'species-ndist-%s' % self.code
        myplot.Save(root=root,
                    xlabel='Number of species',
                    ylabel='Prob',
                    formats=formats,
                    )

    def PlotPrevalences(self, num=5):
        """Plots dist of prevalence for several species.

        num: how many species (starting with the highest prevalence)
        """
        myplot.Clf()
        myplot.PrePlot(num=5)

        for rank in range(1, num+1):
            self.PlotPrevalence(rank)

        root = 'species-prev-%s' % self.code
        myplot.Save(root=root,
                    xlabel='Prevalence',
                    ylabel='Prob',
                    formats=formats,
                    axis=[0, 0.3, 0, 1],
                    )

    def PlotPrevalence(self, rank=1, cdf_flag=True):
        """Plots dist of prevalence one species.

        rank: rank order of the species to plot.
        cdf_flag: whether to plot the CDF
        """
        # convert rank to index
        index = self.GetM() - rank

        pmfs, mix = self.suite.DistOfPrevalence(index)
        count, species = self.GetSpecies(index)
        mix.name = '%d (%d)' % (rank, count)

        print '90%% CI for prevalence of species %d:' % rank, 
        print mix.CredibleInterval(90)

        if cdf_flag:
            cdf = thinkbayes.MakeCdfFromPmf(mix)
            myplot.Cdf(cdf)
        else:
            myplot.Pmf(pmf)

    def PlotMixture(self, rank=1):
        """Plots dist of prevalence for all n, and the mix.

        rank: rank order of the species to plot
        """
        # convert rank to index
        index = self.GetM() - rank

        print self.GetSpecies(index)
        print self.GetCounts()[index]

        pmfs, mix = self.suite.DistOfPrevalence(index)

        myplot.Clf()
        for pmf in pmfs.Values():
            myplot.Pmf(pmf, color='blue', alpha=0.2, linewidth=0.5)

        myplot.Pmf(mix, color='blue', alpha=0.9, linewidth=2)

        root = 'species-mix-%s' % self.code
        myplot.Save(root=root,
                    xlabel='Prevalence',
                    ylabel='Prob',
                    formats=formats,
                    axis=[0, 0.3, 0, 0.3],
                    legend=False)

    def GetSeenSpecies(self):
        """Makes a set of the names of seen species.

        Returns: number of species, set of string species names
        """
        names = self.GetNames()
        m = len(names)
        seen = set(SpeciesGenerator(names, m))
        return m, seen

    def GenerateObservations(self, num_samples):
        """Generates a series of random observations.

        Returns: number of species, sequence of string species names
        """
        n, prevalences = self.suite.Sample()

        names = self.GetNames()
        name_iter = SpeciesGenerator(names, n)

        d = dict(zip(name_iter, prevalences))
        cdf = thinkbayes.MakeCdfFromDict(d)
        observations = cdf.Sample(num_samples)

        return n, observations

    def RunSimulation(self, num_samples, frac_flag=False, delta=0.01):
        """Simulates additional observations and returns a rarefaction curve.

        k is the number of additional observations
        num_new is the number of new species seen

        num_samples: how many new samples to simulate
        frac_flag: whether to convert to fraction of species seen

        Returns: list of (k, num_new) pairs
        """
        m, seen = self.GetSeenSpecies()
        n, observations = self.GenerateObservations(num_samples)

        curve = []
        for k, obs in enumerate(observations):
            seen.add(obs)

            if frac_flag:
                frac_seen = len(seen) / float(n)
                frac_seen += random.uniform(-delta, delta)
                curve.append((k+1, frac_seen))
            else:
                num_new = len(seen) - m
                curve.append((k+1, num_new))

        return curve

    def RunSimulations(self, num_sims, num_samples, frac_flag=False):
        """Runs simulations and returns a list of curves.

        Each curve is a sequence of (k, num_new) pairs.

        num_sims: how many simulations to run
        num_samples: how many samples to generate in each simulation
        frac_flag: whether to convert num_new to fraction of total
        """
        curves = [self.RunSimulation(num_samples, frac_flag) 
                  for i in range(num_sims)]
        return curves

    def MakeJointPredictive(self, curves):
        """Makes a joint distribution of k and num_new.

        curves: list of (k, num_new) curves 

        Returns: joint Pmf of (k, num_new)
        """
        joint = thinkbayes.Joint()
        for curve in curves:
            for k, num_new in curve:
                joint.Incr((k, num_new))
        joint.Normalize()
        return joint

    def MakeFracCdfs(self, curves):
        """Makes Cdfs of the fraction of species seen.

        curves: list of (k, num_new) curves 

        Returns: list of Cdfs
        """
        d = {}
        for curve in curves:
            for k, frac in curve:
                d.setdefault(k, []).append(frac)

        cdfs = {}
        for k, fracs in d.iteritems():
            cdf = thinkbayes.MakeCdfFromList(fracs)
            cdfs[k] = cdf

        return cdfs

    def MakeConditionals(self, curves, ks):
        """Makes Cdfs of the distribution of num_new conditioned on k.

        curves: list of (k, num_new) curves 
        ks: list of values of k

        Returns: list of Cdfs
        """
        joint = self.MakeJointPredictive(curves)

        cdfs = []
        for k in ks:
            pmf = joint.Conditional(1, 0, k)
            pmf.name = 'k=%d' % k
            cdf = thinkbayes.MakeCdfFromPmf(pmf)
            cdfs.append(cdf)
            print '90%% credible interval for %d' % k,
            print cdf.CredibleInterval(90)
        return cdfs


def SpeciesGenerator(names, num):
    """Generates a series of names, starting with the given names.

    Additional names are 'unseen' plus a serial number.

    names: list of strings

    Returns: string iterator
    """
    i = 0
    for name in names:
        yield '%s-%d' % (name, i)
        i += 1

    while i < num:
        yield 'unseen-%d' % i
        i += 1
            

def ReadRarefactedData(filename='journal.pone.0047712.s001.csv'):
    """Reads a data file and returns a list of Subjects.

    Data from http://www.plosone.org/article/
    info%3Adoi%2F10.1371%2Fjournal.pone.0047712#s4

    Returns: list of Subjects
    """
    fp = open(filename)
    reader = csv.reader(fp)
    header = reader.next()
    
    subject = Subject('')
    subject_map = {}

    for t in reader:
        code = t[0]
        if code != subject.code:
            # sort the old subject
            subject.Sort()

            # start a new subject
            subject = Subject(code)
            subject_map[code] = subject

        species = t[1]
        count = int(t[2])
        subject.Add(species, count)

    ComputeNumReads(subject_map)
    return subject_map


def ReadCompleteDataset(filename='BBB_data_from_Rob.csv'):
    """Reads a data file and returns a list of Subjects.

    Data from personal correspondence with Rob Dunn, received 2-7-13.

    Converted from xlsx to csv.

    Returns: map from code to Subject
    """
    fp = open(filename)
    reader = csv.reader(fp)
    header = reader.next()
    header = reader.next()

    subject_codes = header[1:-1]
    subject_codes = ['B'+code for code in subject_codes]

    # print subject_codes

    subject_map = {}
    for code in subject_codes:
        subject_map[code] = Subject(code)
    
    for t in reader:
        otu_code = t[0]
        if otu_code == '':
            continue

        otu_names = t[-1]
        taxons = otu_names.split(';')
        species = taxons[-1]
        counts = [int(x) for x in t[1:-1]]

        # print otu_code, species

        for code, count in zip(subject_codes, counts):
            if count > 0:
                subject_map[code].Add(species, count)

    for code, subject in subject_map.iteritems():
        subject.Sort()

    ComputeNumReads(subject_map)
    return subject_map
        

def ComputeNumReads(subject_map):
    """Computes the number of reads for each subject.

    Adds attributes named num_species and num_reads.

    subject_map: map from subject code to Subject.
    """
    for code, subject in subject_map.iteritems():
        counts = subject.GetCounts()
        subject.num_species = len(counts)
        subject.num_reads = sum(counts)


def JoinSubjects():
    """Reads both datasets and computers their inner join.

    Finds all subjects that appear in both datasets.

    For subjects in the rarefacted dataset, looks up the total
    number of reads and stores it as total_reads.  num_reads
    is normally 400.
    
    """

    # read the rarefacted dataset
    sampled_subjects = ReadRarefactedData()

    # read the complete dataset
    all_subjects = ReadCompleteDataset()

    count = 0
    for code, subject in sampled_subjects.iteritems():
        if code in all_subjects:
            match = all_subjects[code]

            subject.total_reads = match.num_reads
            subject.total_species = match.num_species
            count += 1

    return sampled_subjects


def JitterCurve(curve, dx=0.2, dy=0.3):
    """Adds random noise to the pairs in a curve.

    dx and dy control the amplitude of the noise in each dimension.
    """
    curve = [(x+random.uniform(-dx, dx), 
              y+random.uniform(-dy, dy)) for x, y in curve]
    return curve


def OffsetCurve(curve, i, n, dx=0.3, dy=0.3):
    """Adds random noise to the pairs in a curve.

    i is the index of the curve
    n is the number of curves

    dx and dy control the amplitude of the noise in each dimension.
    """
    xoff = -dx + 2 * dx * i / (n-1)
    yoff = -dy + 2 * dy * i / (n-1)
    curve = [(x+xoff, y+yoff) for x, y in curve]
    return curve


def PlotCurves(curves, root='species-rare'):
    """Plots a set of curves.

    curves is a list of curves; each curve is a list of (x, y) pairs.
    """
    myplot.Clf()
    color = '#225EA8'

    n = len(curves)
    for i, curve in enumerate(curves):
        curve = OffsetCurve(curve, i, n)
        xs, ys = zip(*curve)
        myplot.Plot(xs, ys, color=color, alpha=0.3, linewidth=0.5)

    myplot.Save(root=root,
                xlabel='# samples',
                ylabel='# species',
                formats=formats,
                legend=False)


def PlotConditionals(cdfs, root='species.cond'):
    """Plots cdfs of num_new conditioned on k.

    cdfs: list of Cdf
    root: string filename root
    """
    myplot.Clf()
    myplot.PrePlot(num=len(cdfs))

    myplot.Cdfs(cdfs)

    myplot.Save(root=root,
                xlabel='# new species',
                ylabel='Prob',
                formats=formats)


def PlotFracCdfs(cdfs, root='species.frac'):
    """Plots CDFs of the fraction of species seen.

    cdfs: map from k to CDF of fraction of species seen after k samples
    """
    myplot.Clf()
    color = '#225EA8'

    for k, cdf in cdfs.iteritems():
        if k not in [10, 50] and k % 100:
            continue

        print k
        xs, ys = cdf.Render()
        ys = [1-y for y in ys]
        myplot.Plot(xs, ys, color=color, linewidth=1)

        x = 0.9
        y = 1 - cdf.Prob(x)
        pyplot.text(x, y, str(k), fontsize=9, color=color,
                    horizontalalignment='center',
                    verticalalignment='center',
                    bbox=dict(facecolor='white', edgecolor='none'))

    myplot.Save(root=root,
                xlabel='Fraction of species seen',
                ylabel='Probability',
                formats=formats,
                legend=False)


class Species(thinkbayes.Suite):
    """Represents hypotheses about the number of species."""
    
    def __init__(self, ns, iterations=1000):
        hypos = [thinkbayes.Dirichlet(n) for n in ns]
        thinkbayes.Suite.__init__(self, hypos)
        self.iterations = iterations

    def Update(self, data):
        """Updates the suite based on the data.

        data: list of observed frequencies
        """
        # call Update in the parent class, which calls Likelihood
        thinkbayes.Suite.Update(self, data)

        # update the next level of the hierarchy
        for hypo in self.Values():
            hypo.Update(data)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under this hypothesis.

        hypo: Dirichlet object
        data: list of observed frequencies
        """
        dirichlet = hypo

        # draw sample Likelihoods from the hypothetical Dirichlet dist
        # and add them up
        like = 0
        for i in range(self.iterations):
            like += dirichlet.Likelihood(data)

        # correct for the number of ways the observed species
        # might as been chosen from all species
        m = len(data)
        like *= thinkbayes.BinomialCoef(dirichlet.n, m)

        return like

    def DistOfN(self):
        """Computes the distribution of n."""
        pmf = thinkbayes.Pmf()
        for hypo, prob in self.Items():
            pmf.Set(hypo.n, prob)
        return pmf
        

class Species2(object):
    """Represents hypotheses about the number of species.

    Combines two layers of the hierarchy into one object.

    ns and probs represent the distribution of N

    params represents the parameters of the Dirichlet distributions
    """
    
    def __init__(self, ns, iterations=1000):
        self.ns = ns
        self.probs = numpy.ones(len(ns), dtype=numpy.double)
        self.params = numpy.ones(self.ns[-1], dtype=numpy.int)
        self.iterations = iterations

    def GetPmf(self):
        """Returns a Pmf of possible number of species."""
        items = zip(self.ns, self.probs)
        pmf = thinkbayes.MakePmfFromItems(items)
        return pmf

    def Update(self, data):
        like = numpy.zeros(len(self.ns), dtype=numpy.double)
        for i in range(1000):
            like += self.SampleLikelihood(data)

        self.probs *= like
        self.probs /= self.probs.sum()

        m = len(data)
        self.params[:m] += data

    def SampleLikelihood(self, data):
        """Computes the likelihood of the data for all values of n.

        Draws one sample from the distribution of prevalences.

        data: sequence of observed counts

        Returns: numpy array of m likelihoods
        """
        gammas = numpy.random.gamma(self.params)

        m = len(data)
        row = gammas[:m]
        col = numpy.cumsum(gammas)

        log_likes = []
        for n in self.ns:
            ps = row / col[n-1]
            terms = numpy.log(ps) * data
            log_like = terms.sum()
            log_likes.append(log_like)

        log_likes -= numpy.max(log_likes)
        likes = numpy.exp(log_likes)

        coefs = [thinkbayes.BinomialCoef(n, m) for n in self.ns]
        likes *= coefs

        return likes

    def DistOfN(self):
        """Computes the distribution of n.

        Returns: new Pmf object
        """
        pmf = thinkbayes.MakePmfFromDict(dict(zip(self.ns, self.probs)))
        return pmf

    def MarginalBeta(self, n, index):
        """Computes the conditional distribution of the indicated species.
        
        n: conditional number of species
        index: which species

        Returns: Beta object representing a distribution of prevalence.
        """
        alpha0 = self.params[:n].sum()
        alpha = self.params[index]
        return thinkbayes.Beta(alpha, alpha0-alpha)

    def DistOfPrevalence(self, index):
        """Computes the distribution of prevalence for the indicated species.

        index: which species

        Returns: (pmfs, mix) where pmfs is a MetaPmf and mix is a Pmf
        """
        pmfs = thinkbayes.Pmf()

        for n, prob in zip(self.ns, self.probs):
            beta = self.MarginalBeta(n, index)
            pmf = beta.MakePmf()
            pmfs.Set(pmf, prob)

        mix = thinkbayes.MakeMixture(pmfs)
        return pmfs, mix
        
    def Sample(self):
        """Draws random n and prevalences.

        Returns: (n, prevalences)
        """
        pmf = self.DistOfN()
        n = pmf.Random()
        prevalences = self.SampleConditional(n)
        return n, prevalences

    def SampleConditional(self, n):
        """Draws a sample of prevalences given n.

        n: the number of species assumed in the conditional

        Returns: numpy array of n prevalences
        """
        params = self.params[:n]
        gammas = numpy.random.gamma(params)
        gammas /= gammas.sum()
        return gammas
        

class Species3(Species2):
    """Represents hypotheses about the number of species."""
    
    def Update(self, data):
        """Updates the suite based on the data.

        data: list of observations
        """
        # sample the likelihoods and add them up
        like = numpy.zeros(len(self.ns), dtype=numpy.double)
        for i in range(self.iterations):
            like += self.SampleLikelihood(data)

        self.probs *= like
        self.probs /= self.probs.sum()

        m = len(data)
        self.params[:m] += data

    def SampleLikelihood(self, data):
        """Computes the likelihood of the data under all hypotheses.

        data: list of observations
        """
        # get a random sample
        gammas = numpy.random.gamma(self.params)

        # row is just the first m elements of gammas
        m = len(data)
        row = gammas[:m]

        # col is the cumulative sum of gammas
        col = numpy.cumsum(gammas)[self.ns[0]-1:]

        # each row of the array is a set of ps, normalized
        # for each hypothetical value of n
        array = row / col[:, numpy.newaxis]

        # computing the multinomial PDF under a log transform
        # take the log of the ps and multiply by the data
        terms = numpy.log(array) * data

        # add up the rows
        log_likes = terms.sum(axis=1)

        # before exponentiating, scale into a reasonable range
        log_likes -= numpy.max(log_likes)
        likes = numpy.exp(log_likes)

        # correct for the number of ways we could see m species
        # out of a possible n
        coefs = [thinkbayes.BinomialCoef(n, m) for n in self.ns]
        likes *= coefs

        return likes


class Species4(Species):
    """Represents hypotheses about the number of species."""
    
    def Update(self, data):
        """Updates the suite based on the data.

        data: list of observed frequencies
        """
        m = len(data)

        # loop through the species and update one at a time
        for i in range(m):
            one = numpy.zeros(i+1)
            one[i] = data[i]
            
            # call the parent class
            Species.Update(self, one)

    def Likelihood(self, hypo, data):
        """Computes the likelihood of the data under this hypothesis.

        Note: this only works correctly if we update one species at a time.

        hypo: Dirichlet object
        data: list of observed frequencies
        """
        dirichlet = hypo
        like = 0
        for i in range(self.iterations):
            like += dirichlet.Likelihood(data)

        # correct for the number of unseen species the new one
        # could have been
        m = len(data)
        num_unseen = dirichlet.n - m + 1
        like *= num_unseen

        return like


class Species5(Species2):
    """Represents hypotheses about the number of species.

    Combines two laters of the hierarchy into one object.

    ns and probs represent the distribution of N

    params represents the parameters of the Dirichlet distributions
    """
    
    def Update(self, data):
        """Updates the suite based on the data.

        data: list of observed frequencies
        """
        # loop through the species and update one at a time
        m = len(data)
        for i in range(m):
            self.UpdateOne(i+1, data[i])
            self.params[i] += data[i]

    def UpdateOne(self, m, count):
        """Updates the suite based on the data.

        Evaluates the likelihood for all values of n.

        m: which species was observed
        count: how many were observed
        """
        # sample the likelihoods and add them up
        likes = numpy.zeros(len(self.ns), dtype=numpy.double)
        for i in range(self.iterations):
            likes += self.SampleLikelihood(m, count)

        # correct for the number of unseen species the new one
        # could have been
        unseen_species = [n-m+1 for n in self.ns]
        likes *= unseen_species

        # multiply the priors by the likelihoods and renormalize
        self.probs *= likes
        self.probs /= self.probs.sum()

    def SampleLikelihood(self, m, count):
        """Computes the likelihood of the data under all hypotheses.

        m: which species was observed
        count: how many were observed
        """
        # get a random sample of p
        gammas = numpy.random.gamma(self.params)

        # sums is the cumulative sum of p, for each value of n
        sums = numpy.cumsum(gammas)[self.ns[0]-1:]

        # get p for the mth species, for each value of n
        ps = gammas[m-1] / sums
        log_likes = numpy.log(ps) * count

        # before exponentiating, scale into a reasonable range
        log_likes -= numpy.max(log_likes)
        likes = numpy.exp(log_likes)

        return likes


class Species6(Species):
    """Represents hypotheses about the number of species."""
    
    def __init__(self, ns):
        hypos = [OversampledDirichlet(n) for n in ns]
        thinkbayes.Suite.__init__(self, hypos)

    def Update(self, data):
        """Updates the suite based on the data.

        data: list of observed frequencies
        """
        for hypo in self.Values():
            hypo.Peek(data)

        # call Update in the parent class, which calls Likelihood
        thinkbayes.Suite.Update(self, data)

        for hypo in self.Values():
            hypo.Update(data)


class OversampledDirichlet(thinkbayes.Dirichlet):
    """Provides oversampled random selection from a Dirichlet distribution.

    By peeking at the posterior distribution, we can see where the
    non-negligible mass will be, then restrict the uniform prior to
    just this range.

    This version doesn't provide Update, because it only handles the
    special case when the prior distribution is still uniform.

    """
    
    def Likelihood(self, data):
        like = thinkbayes.Dirichlet.Likelihood(self, data)
        return like * self.weight

    def Peek(self, data):
        """Figures out the bounds in the posterior distribution that
        contain non-negligible mass.
        """
        m = len(data)
        priors = self.MakeConditionals(self.params)
        pmfs = [prior.MakePmf() for prior in priors]
        self.conditionals = pmfs

        self.TrimMarginals(data)

    def TrimMarginals(self, data):
        m = len(data)
        params = numpy.copy(self.params)
        params[:m] += data

        posteriors = self.MakeMarginals(params)
        cdfs = [beta.MakeCdf() for beta in posteriors]

        lows = [cdf.Percentile(2) for cdf in cdfs]
        highs = [cdf.Percentile(98) for cdf in cdfs]
        
        p_hit = 1
        pmfs = self.conditionals
        for pmf, low, high in zip(pmfs, lows, highs):
            mass = self.TrimMarginal(pmf, low, high)
            print self.n, low, high, mass
            p_hit *= mass

        self.weight = p_hit
        print self.n, self.weight
        conditionals = [thinkbayes.MakeCdfFromPmf(pmf) for pmf in pmfs]
        self.conditionals = conditionals

    def TrimMarginal(self, pmf, low, high):
        for val in pmf.Values():
            if val < low or val > high:
                pmf.Remove(val)
        return pmf.Normalize()

    def MakeMarginals(self, params):
        total = sum(params)
        betas = [thinkbayes.Beta(x, total-x) for x in params]
        return betas

    def MakeConditionals(self, params):
        total = sum(params)
        betas = []
        for x in params[:-1]:
            total -= x
            beta = thinkbayes.Beta(x, total)
            betas.append(beta)
        return betas

    def Random(self):
        """Generates a random variate from this distribution.

        Returns: normalized array of probabilities
        """
        fraction = 1.0
        ps = numpy.zeros(self.n)

        for i, cond in enumerate(self.conditionals):
            p = cond.Random()
            ps[i] = p * fraction
            fraction *= 1-p

        ps[-1] = fraction
        return ps


def MakePosterior(constructor, data, ns, iterations=1000):
    """Makes a suite, updates it and returns the PMF of N.

    Prints the elapsed time.

    data: observed species and their counts
    ns: sequence of hypothetical ns
    iterations: how many samples to draw
    """
    suite = constructor(ns, iterations)

    print constructor.__name__
    start = time.time()
    suite.Update(data)
    end = time.time()

    print end-start

    pmf = suite.DistOfN()
    return pmf


def PlotAllVersions():
    """Makes a graph of posterior distributions of N."""
    data = [1, 2, 3]
    m = len(data)
    n = 20
    ns = range(m, n)

    for constructor in [Species, Species2, Species3, Species4, Species5]:
        pmf = MakePosterior(constructor, data, ns)
        pmf.name = '%s' % (constructor.__name__)
        myplot.Pmf(pmf)

    myplot.Show()
    return

    myplot.Save(root='species3',
                xlabel='Number of species',
                ylabel='Prob')


def PlotMedium():
    """Makes a graph of posterior distributions of N."""
    data = [1, 1, 1, 1, 2, 3, 5, 9]
    m = len(data)
    n = 20
    ns = range(m, n)

    for constructor in [Species, Species2, Species3, Species4, Species5]:
        pmf = MakePosterior(constructor, data, ns)
        pmf.name = '%s' % (constructor.__name__)
        myplot.Pmf(pmf)

    myplot.Show()


def CompareHierarchicalExample():
    """Makes a graph of posterior distributions of N."""
    data = [3, 2, 1]
    m = len(data)
    n = 30
    ns = range(m, n)

    constructors = [Species, Species5]
    iterations = [1000, 100]

    for constructor, iterations in zip(constructors, iterations):
        pmf = MakePosterior(constructor, data, ns, iterations)
        pmf.name = '%s' % (constructor.__name__)
        myplot.Pmf(pmf)

    myplot.Show()


def SimpleDirichletExample():
    """Makes a plot showing posterior distributions for three species.

    This is the case where we know there are exactly three species.
    """
    myplot.Clf()
    myplot.PrePlot(3)

    names = ['lions',  'tigers', 'bears']
    data = [3, 2, 1]

    dirichlet = thinkbayes.Dirichlet(3)
    for i in range(3):
        beta = dirichlet.MarginalBeta(i)
        print 'mean', names[i], beta.Mean()

    dirichlet.Update(data)
    for i in range(3):
        beta = dirichlet.MarginalBeta(i)
        print 'mean', names[i], beta.Mean()

        pmf = beta.MakePmf(name=names[i])
        myplot.Pmf(pmf)

    myplot.Save(root='species1',
                xlabel='Prevalence',
                ylabel='Prob',
                formats=formats,
                )


def HierarchicalExample():
    """Shows the posterior distribution of n for lions, tigers and bears.
    """
    ns = range(3, 30)
    suite = Species(ns, iterations=8000)

    data = [3, 2, 1]
    suite.Update(data)

    myplot.Clf()
    myplot.PrePlot(num=1)

    pmf = suite.DistOfN()
    myplot.Pmf(pmf)
    myplot.Save(root='species2',
                xlabel='Number of species',
                ylabel='Prob',
                formats=formats,
                )


def ProcessSubjects(codes):
    """Process subjects with the given codes and plot their posteriors.

    code: sequence of string codes
    """
    myplot.Clf()
    myplot.PrePlot(len(codes))

    subjects = ReadRarefactedData()
    pmfs = []
    for code in codes:
        subject = subjects[code]

        subject.Process()
        pmf = subject.suite.DistOfN()
        pmf.name = subject.code
        myplot.Pmf(pmf)

        pmfs.append(pmf)

    print 'ProbGreater', thinkbayes.PmfProbGreater(pmfs[0], pmfs[1])
    print 'ProbLess', thinkbayes.PmfProbLess(pmfs[0], pmfs[1])

    myplot.Save(root='species4',
                xlabel='Number of species',
                ylabel='Prob',
                formats=formats,
                )


def RunSubject(code):
    """Run the analysis for the subject with the given code.

    code: string code
    """
    subjects = ReadRarefactedData()
    subject = subjects[code]
    subject.Process()
    subject.MakeFigures()

    num_samples = 400
    curves = subject.RunSimulations(100, num_samples)
    root = 'species-rare-%s' % subject.code
    PlotCurves(curves, root=root)

    num_samples = 800
    curves = subject.RunSimulations(500, num_samples)
    ks = [100, 200, 400, 800]
    cdfs = subject.MakeConditionals(curves, ks)
    root = 'species-cond-%s' % subject.code
    PlotConditionals(cdfs, root=root)

    curves = subject.RunSimulations(500, num_samples, frac_flag=True)
    cdfs = subject.MakeFracCdfs(curves)
    root = 'species-frac-%s' % subject.code
    PlotFracCdfs(cdfs, root=root)


def MakePredictions(lockername='species_locker.db', num=3, replace=False):
    """Make predictions for each subject in the list and store in a locker.

    locker: string locker file name
    num: how many subjects to process
    """
    subject_map = JoinSubjects()
    locker = Locker(lockername)

    i = 0
    for code, subject in subject_map.iteritems():
        print code, '...',

        if not replace:
            processed = locker.Lookup(code)
            if processed is not None:
                print 'In cache'
                continue

        print 'Processing'
        MakePrediction(subject)
        locker.Add(subject.code, subject)

        i += 1
        if i == num:
            break

    locker.Close()


def MakePrediction(subject, num_sims=300):
    """Make predictions for the given subject.

    subject: Subject object
    num_sims: how many simulations to run for predictions

    Adds attributes
    ps: list of probabilities
    cis: list of credible intervals
    """
    subject.Process()

    print subject.code, subject.num_reads, subject.total_reads

    add_reads = subject.total_reads - subject.num_reads
    curves = subject.RunSimulations(num_sims, add_reads)
    #PlotCurves(curves, root='species.subject')
    cdfs = subject.MakeConditionals(curves, [add_reads])
    cdf = cdfs[0]

    ps = range(10, 100, 10)
    cis = [cdf.CredibleInterval(p) for p in ps]

    for p, ci in zip(ps, cis):
        print p, ci

    subject.ps = ps
    subject.cis = cis


def MakePredictionTable(subject_map):
    """Makes a table of predictions in LaTeX format.

    subject_map: map from code to Subject

    Precondition: subject have attributes ps and cis
    """
    def PrintRow(t):
        print '&',
        t = [str(x) for x in t]
        print ' & '.join(t),
        print r'\\'
        
    i = 0
    for subject in subject_map.itervalues():
        if i == 0:
            PrintRow(subject.ps)

        try:
            PrintRow(subject.cis)
        except AttributeError:
            break

        i += 1


def WritePredictionHTML(subject_map, filename='species_table.html'):
    """Write the prediction table to a file.

    subject_map: map from code to Subject
    """
    t = MakePredictionHTML(subject_map)
    fp = open(filename, 'w')
    fp.write(str(t))
    fp.close()


def MakePredictionHTML(subject_map):
    """Makes a table of predictions in LaTeX format.

    subject_map: map from code to Subject

    Precondition: subject have attributes ps and cis
    """
    subject = subject_map.values()[0]
    header = ['Code', '# reads', '# species'] + subject.ps
    t = HTML.Table(header_row=header)

    for code, subject in sorted(subject_map.iteritems()):
        names = subject.GetNames()
        m = len(names)

        row = [code, subject.total_reads, m] + subject.cis
        t.rows.append(row)

    return t


def SummarizeData():
    """Read data and print subject codes and number of species."""
    subject_map = JoinSubjects()

    for code, subject in subject_map.iteritems():
        print subject.code, subject.total_species
        #print subject.num_reads, subject.num_species
        #print subject.total_reads, subject.total_species


def ValidatePredictions(lockername='species_locker.db'):
    locker = Locker(lockername)
    subject_map = locker.Read()
    locker.Close()

    for code, subject in subject_map.iteritems():
        print subject.code
        print subject.num_reads, subject.num_species
        print subject.total_reads, subject.total_species

        pmf = subject.suite.GetPmf()
        print pmf.Mean(), pmf.CredibleInterval(90)

        myplot.Pmf(pmf)
        myplot.Show()


def main(script, *args):
    random.seed(17)

    # MakePredictions(num=60, replace=True)
    ValidatePredictions()
    return

    SummarizeData()
    return

    SimpleDirichletExample()
    HierarchicalExample()
    return

    RunSubject(4)
    return

    SummarizeData()
    return

    CompareHierarchicalExample()
    return

    PlotMedium()
    return

    random.seed(17)
    PlotAllVersions()
    return

    pmf = MakePosterior(Species)
    myplot.Pmf(pmf)
    myplot.Show()
    return


    random.seed(17)
    p = dirichlet.Likelihood(data)
    print p

    random.seed(17)
    lp = dirichlet.LogLikelihood(data)
    print lp, math.exp(lp)

    return


    xs = []
    for i in range(ITERATIONS):
        x = dirichlet.Random()
        xs.append(x)

    cols = zip(*xs)
    for col in cols:
        cdf = thinkbayes.MakeCdfFromList(col)
        myplot.Cdf(cdf)

    myplot.Show()
    return


    beta = thinkbayes.Beta()
    beta.Update((140, 110))
    print beta.Mean()

    pmf = beta.MakePmf(101)
    myplot.Pmf(pmf)
    myplot.Show()
    return


    flag = int(flag)

    d = {
        0: [1],
        1: [9, 6],
        2: [8, 3],
        3: [11, 4, 2],
        4: [3, 2, 2, 1],
        5: [7, 3, 1, 1, 1, 1],
        6: [7, 2, 1],
        7: [4, 3, 9, 1, 1],
        8: [14, 5, 2, 2, 1,1,1,1,1,1],
        9: [11, 2, 1,1],
        }
    freqs = d[flag]

    sample = MakeSample(freqs)
    meta = MakeMeta(sample, n=31)

    meta.Print()

    root='rare%d' % flag
    MakeFigures(root, sample, meta)


if __name__ == '__main__':
    profile = False
    if profile:
        import cProfile
        cProfile.run('main(*sys.argv)')
    else:
        main(*sys.argv)
