# the following is needed to import Gui
import sys
sys.path.insert(1, '/home/downey/python')
from Gui import *

from math import *
from random import *
import time
import getopt

"""
p(i+) is the probability that the last changepoint is at i
      (meaning that i is the first observation after the changepoint)

p(j++) is the probability that the second-to-last changepoint is
       at j

"""

class Slider:
    """each slider keeps track of the cumulative sums starting
    from the given index and computes the probabilities p(i+|index++)
    for each i between index+2 and t-2, where t is the index of the
    most recent observation
    """
    def __init__(self, index):
        self.index = index
        self.n = 0
        self.S1 = []
        self.S2 = []

    def update(self, x):
        self.n += 1
        try:
            self.S1.append(self.S1[-1] + x)
            self.S2.append(self.S2[-1] + x*x)
        except IndexError:
            self.S1.append(x)
            self.S2.append(x*x)

    def moments(self, i):
        """return the moments of the set of points to the left of i
        (not including i)"""
        return [ i, self.S1[i-1], self.S2[i-1] ]

    def normalize(self, factor=1.0):
        nc = sum(self.like) / factor
        if nc == 0: return
        for i in xrange(self.n):
            self.like[i] /= nc

    def dump(self):
        print '\n#', self.index
        total = 0
        for i, x in enumerate(self.like):
            total += x
            print i+self.index, total

    def pjpp(self, i):
        """return p(i+|j++)

        Precondition: likelihood has run since the last update.
        """
        return self.like[i-self.index]
    
    
    def pcpp(self, i):
        """return p(i+|jc), where jc means there have been 0-1
        changepoints

        Precondition: likelihood0 has run since the last update.
        """
        return self.cpp[i-self.index]
    
    def likelihood0(self, f):
        """compute likelihood(i+|jc), where jc is the hypothesis
        that there are 0-1 changepoints in the region.
        """
        self.likelihood(f)
        self.cpp = self.like

        #n = len(self.cpp)-1
        #if 26 <= n <= 31:
        #    dump_list(self.cpp, 'cpp.', n)

        #if n == 26 and False:
        #    print_flag = True
        #else:
        #    print_flag = False

            
    def likelihood(self, f=None):
        """compute likelihood(i+|j++), where j = self.index.
        f is the prior probability of a changepoint at a given
        index; if f is provided, it means that we should also
        consider the hypothesis that there is no changepoint
        in the region.
        """
        n = self.n
        self.like = [float('-Inf')] * n
        if n < 4: return

        if simple:
            ll = loglike
        else:
            ll = loglike2
                
        m0 = self.moments(n)
        if f:
            if sigma_known:
                s2 = known_sigma**2

                if mu_known:
                    like0 = ll(m0[0], m0[1], m0[2],
                                           s2=s2, mu0=known_mu)
                else:
                    like0 = ll(m0[0], m0[1], m0[2], s2=s2)

            elif mu_known:
                like0 = ll(m0[0], m0[1], m0[2], mu0=known_mu)
                
            elif pooled:
                # use the pooled variance
                sp2 = spooled(*m1 + m2)
                like0 = ll(m0[0], m0[1], m0[2], s2=sp2)

            else:
                # both parameters unknown -- drawn them from
                # posterior distributions based on estimates
                like0 = ll(*m0)

        for i in xrange(2,n-1):
            m1 = self.moments(i)
            m2 = sub_moments(m0, m1)

            if f and self.n==68 and i==66 and False:
                flag = True
            else:
                flag = False

            if sigma_known:
                s2 = known_sigma**2

                if mu_known:
                    like1 = ll(m1[0], m1[1], m1[2],
                                        s2=s2, mu0=known_mu)
                    
                else:
                    like1 = ll(m1[0], m1[1], m1[2], s2=s2)
                    
                like2 = ll(m2[0], m2[1], m2[2], s2=s2)

            elif mu_known:
                like1 = ll(m1[0], m1[1], m1[2], mu0=known_mu)
                like2 = ll(m2[0], m2[1], m2[2])
                
            elif pooled:
                # use the pooled variance
                sp2 = spooled(*m1 + m2)

                # uncomment these lines to choose from the posterior
                # distribution of sigma^2
                #n = m0[0]
                #sp2 = scaled_inverse_chi2variate(n-1.0, sp2)

                #like1 = ll(m1[0], m1[1], m1[2], s2=sp2)
                #like2 = ll(m2[0], m2[1], m2[2], s2=sp2)
                like1 = ll(m1[0], m1[1], m1[2], s2=sp2)
                like2 = ll(m1[0], m1[1], m1[2], s2=sp2)

            else:
                # both parameters unknown -- drawn them from
                # posterior distributions based on estimates
                
                like1 = ll(*m1)
                like2 = ll(*m2)
                if f and print_flag:
                    print 'like1', i, like1
                    print 'like2', i, like2
                    print 'sum', i, like1 + like2

            self.like[i] = like1+like2

        # the first two and last two entries are 0 because we
        # can't compute the likelihood of 0 or 1 data points.

        if f:
            # the likelihood for the zero-changepoint hypothesis
            # gets tucked into like[0]
            self.like[0] = like0

        #self.like[1] = self.like[2] 
        #self.like[n-1] = self.like[n-2] 
        
        # convert likelihoods from log to linear
        # (scale the likelihoods to avoid overflow when
        # we exponentiate)
        m = max(self.like)
        self.like = [exp(x-m) for x in self.like]

        if f:
            prior = (1-f) / f

            # need to get back and see what's wrong with the following
            #fc = pow(1-f, n)
            #prior = n * fc / (1-fc)
            #self.like[0] *= prior

        # the following factor compensates for the lumping
        # of probability when n is small, but I'm not sure
        # it's not a kludge
        # factor = (n-3.0) / (n-1.0)
        factor = 1.0
        
        self.normalize(factor)
        self.like[0] = 0


def spooled(m0, m1, m2, n0, n1, n2):
    # compute the spooled variance of the two data sets described
    # by moments (m0, m1, m2) and (n0, n1, n2)
    mum = m1 / m0
    mun = n1 / n0
    tm = m2 - mum**2 / m0
    tn = n2 - mun**2 / n0
    sp2 = (tm + tn) / (n0 + m0 - 1)
    return sp2


def var(m0, m1, m2):
    mu = m1 / m0
    return m2 / m0 - mu**2

def spooled2(m0, m1, m2, n0, n1, n2):
    # compute the spooled variance of the two data sets described
    # by moments (m0, m1, m2) and (n0, n1, n2)
    v1 = var(m0, m1, m2)
    v2 = var(n0, n1, n2)
    sp2 = (m0-1) * v1 + (n0-1) * v2 / (n0 + m0 - 2)
    return sp2


def loglike(n, m1, m2, mu0=None, s2=None):
    mu = m1 / n
    var = m2 / n - mu**2
    #like = -n*log(var)/2 - n/2.0
    like = -n*log(var)/2
    return like


def loglike2(n, m1, m2, flag=False, s2=None, mu0=None):
    """compute the log likelihood of the given moments about the origin
    """
    
    # compute the sample mean
    ybar = m1 / n

    # if we don't know the variance, estimate it
    if s2 == None:
        s2 = m2 / n - ybar**2
        s2 = s2 * n / (n-1)
        s2 = max(s2, 1e-16) 
        
        # draw sig2 from the its posterior distribution
        sig2 = scaled_inverse_chi2variate(n-1.0, s2)
    else:
        # otherwise use the known value
        sig2 = s2

    # if we don't know mu0, draw it from the posterior distribution
    if mu0 == None:
        mu = normalvariate(ybar, sqrt(sig2/n))
    else:
        mu = mu0

    # when
    # mu = ybar and
    # sig2 = m2 / n - ybar**2
    # then cm/sig2/2 = n/2

    # compute log likelihood, ignoring the -n/2 log(2 pi) term,
    # which would get divided out during normalization anyway.
    cm = m2 - 2*mu*m1 + mu*mu*n
    like = -n*log(sig2)/2 - cm/sig2/2

    if print_flag:
        print -n*log(sig2)/2, - cm/sig2/2,

    #if flag:
    #    print '\nlike2'
    #    print n, m1, m2
    #    print ybar, s2
    #    print mu, sig2, sqrt(sig2/n)
    #    print m2, - 2*mu*m1, mu*mu*n, cm
    #    print -n*log(sig2)/2, -cm/sig2/2
    #    print like

    return like


print_flag = False

def print_stats(n, s1, s2):
    mu = s1 / n
    var = s2 / n - mu**2
    like = -n * log(var) / 2
    print mu, var, like


def calc_moments(t):
    n = len(t)
    s1 = sum(t)
    s2 = sum([z*z for z in t])
    return [n, s1, s2]

def add_moments(c1, c2):
    """return the sum of two moment vectors"""
    return [m1+m2 for m1, m2 in zip(c1, c2)]

def sub_moments(c1, c2):
    """return the difference of two moment vectors"""
    a, b, d = c1
    e, f, g = c2
    return [a-e, b-f, d-g]
    #return [c1[i]-c2[i] for i in range(3)]
    #return [m1-m2 for m1, m2 in zip(c1, c2)]


class Estimate:
    def __init__(self, n, p, pp):
        self.n = n
        self.p = p
        self.pp = pp

    def ppkp(self, i):
        return self.pp[i]

    def dump(self):
        print '\n#', self.n
        for i, x in enumerate(self.p):
            print i, x

        print '\n#', self.n
        for i, x in enumerate(self.pp):
            print i, x


def copy_est(previous, f):
    n = previous.n + 1
    p = previous.p[:] + [0]
    pp = previous.pp[:] + [0]
    return Estimate(n, p, pp)


class Series:
    def __init__(self, f=0.02):
        self.f = f * fudge
        self.data = []
        self.sliders = {}
        self.estimates = {}
        self.init_est()

    def init_est(self):
        """initialize the first three elements of self.estimates
        """
        for i in range(0, 4):
            self.estimates[i] = self.make_est(i)
        
    def make_est(self, n):
        """make estimates of p and pp based on no data, only the
        prior probability, f, of a changepoint at any time.
        """
        f = self.f
        p = [f * pow(1-f, n-i-1) for i in range(n)]
        pp = [0] * n
        
        for i in range(n-1):
            t = [p[k] * f * pow(1-f, k-i-1)
                 for k in range(i+1, n)]
            pp[i] = sum(t)

        return Estimate(n, p, pp)

    
    def add(self, x):
        self.data.append(x)
        i = len(self.data)-1
        self.sliders[i] = Slider(i)

        for sl in self.sliders.itervalues():
            sl.update(x)


    def process(self, p=False, pp=False):
        # the earliest slider needs to compute the
        # zero-changepoint probabilities
        i0 = min(self.sliders)
        s0 = self.sliders[i0]
        s0.likelihood0(self.f)

        for sl in self.sliders.itervalues():
            sl.likelihood()

        self.calc_prob()
    

    def calc_prob(self):
        """each slider corresponds to the hypothesis that the second
        breakpoint is at j.

        The jth slider can compute p(i+|j++) for any i
        
        p(i+) = sum_j p(i+|j++) * p(j++)  for all j less than i
        """
        i = len(self.data)-1
        
        if i<3: return
        
        previous = self.estimates[i-1]
        
        self.estimates[i] = copy_est(previous, self.f)

        #print i, self.estimates[i].p

        for j in range(iters):
            self.iterate_est(i)


    def iterate_est(self, n):
        est = self.estimates
        p = est[n].p
        pp = est[n].pp
        f = self.f
        sliders = self.sliders

        # update p
        cpp = 1 - sum(pp)
        #print cpp

        #if n==25 or n==26 or n==27:
        #    t = [sliders[0].pcpp(i+1) for i in range(0, n-1)]
        #    dump_list(t, 'pcpp.', n)


        for i in range(0, n-1):
            t = [pp[j] * sliders[j].pjpp(i) for j in range(0, i-1)]
            p[i] = sum(t)
            p[i] += cpp * sliders[0].pcpp(i+1)

                    
        # update pp
        for i in range(0, n-2):
            t = [p[k] * est[k].p[i] for k in range(i+1, n)]
            # t2 = [est[k].p[i] for k in range(i+1, n)]
            pp[i] = sum(t)


    def make_cumulants(self, win=10):
        """compute cump, which is a cumulative sum of the
        estimated p(i+), added backwards starting with the
        most recent observation, and cumpp, which is the
        cumulative sum of p(i++)
        """
        n = len(self.data)-1

        p = self.estimates[n].p
        pp = self.estimates[n].pp
        winp = [0] * n
        cump = [0] * n
        cumpp = [0] * n

        t = 0
        for i in range(n):
            t += p[i]
            winp[i] = t
            if i >= win:
                t -= p[i-win]
            
        tp = 0
        tpp = 0
        for i in range(n-1, -1, -1):
            tp += p[i]
            cump[i] = tp
            tpp += pp[i]
            cumpp[i] = tpp

        self.winp, self.cump, self.cumpp = winp, cump, cumpp
            

    def print_like(self, step=10):
        n = self.n
        for i in xrange(0, n, step):
            self.sliders[i].dump()


    def print_p(self):
        n = self.n
        p = self.estimates[n].p

        print '\n#prob'
        total = 0
        for i in range(n-1, -1, -1):
            total += p[i]
            #print i, p[i]
            print i, total


    def print_pp(self):
        n = self.n
        pp = self.pp

        print '\n#prob+'
        total = 0
        for i in range(n-1, -1, -1):
            total += pp[i]
            #print i, pp[i]
            print i, total
 

def cumulant(t):
    """return the cumulant sum of t (forward)
    """
    s = t[:]
    sum = 0
    for i in range(len(t)):
        sum += t[i]
        s[i] = sum
    return s
        

def read_data(filename):
    t = []
    for line in file(filename):
        try:
            x = float(line)
            t.append(x)
        except ValueError:
            pass
    return t


def generate_data(mus=[0, 1], sigma=1.0, n=100):
    data = []
    for mu in mus:
        data.extend([normalvariate(mu, sigma) for i in xrange(n)])
    return data


def generate_var_data(mu=0.0, sigmas=[1.0, 0.5], n=40):
    data = []
    data.extend([normalvariate(1.0, 1.0) for i in xrange(n)])
    
    for sigma in sigmas:
        data.extend([normalvariate(mu, sigma) for i in xrange(n)])
    return data


def generate_trial(n, m=100, sigma=1.0):
    """generate a dataset with a changepoint after n.
    Before the changepoint, datapoints are N(mu0, sigma); after
    the changepoint, generate n datapoints with (mu1, sigma)
    """
    mu0 = 0.0
    mu1 = 1.0
    
    data = []
    data.extend([normalvariate(mu0, sigma) for i in xrange(n)])
    data.extend([normalvariate(mu1, sigma) for i in xrange(m)])

    return data


def generate_trend_data(sigma=1.0, n=100):
    t = []
    for i in range(n):
        mu = 1.0 * i / n
        x = normalvariate(mu,sigma)
        t.append(x)
    return t


def chi2variate(nu):
    """generate a single random value from the chi-square
    distribution with parameter nu by generating a value
    from the Gamma distribution with parameters (nu/2, 2).
    See http://en.wikipedia.org/wiki/Gamma_distribution
    """
    return gammavariate(nu/2, 2)


def scaled_inverse_chi2variate(nu, s2):
    """generate a single random value from the scaled
    inverse chi-square
    distribution with parameters (nu, s^2) by generating a value
    from the chi square distribution with parameter (nu).
    See http://en.wikipedia.org/wiki/Scale-inverse-chi-square_distribution
    """
    x = gammavariate(nu/2, 2)
    y = s2 * nu / x
    return y

def test_chi2():
    """confirm that the distribution generated by chi2variate(k)
    is the same as the distribution of the sum of square of 10
    normal variates with (0, 1)
    """
    
    k = 10
    n = 100000

    fp = open('sums', 'w')
    for i in range(n):
        t = [normalvariate(0,1) for j in range(k)]
        t = [x*x for x in t]
        s = sum(t)
        fp.write('%f\n' % (s,))
    fp.close()

    fp = open('chi2', 'w')
    for i in range(n):
        s = chi2variate(k)
        fp.write('%f\n' % (s,))
    fp.close()


def dump_list(t, prefix, sigma):
    filename = prefix + str(sigma)
    fp = open(filename, 'w')
    for i, x in enumerate(t):
        s = '%d\t%f\n' % (i, x)
        fp.write(s)
    fp.close()

def dump_lists(xs, ys, prefix, sigma):
    filename = prefix + str(sigma)
    fp = open(filename, 'w')
    for x, y in zip(xs, ys):
        s = '%f\t%f\n' % (x, y)
        fp.write(s)
    fp.close()


class GraphTransform(Transform):
    """the origin is in the middle of the leftmost column
    the positive y-axis is up, and range is the range of the y axis
    """
    def __init__(self, height, range=[-1, 1], xscale=5):
        self.xscale = xscale
        self.height = height
        self.min, self.max = range
        self.dy = 1.0 * (self.max - self.min)
    
    def trans(self, p):
        x, y = p
        x *= self.xscale
        y = self.height - self.height * (y - self.min) / self.dy
        return [x, y]



class ScrollableCanvas:
    def __init__(self, gui, width=500, height=500, **options):
        self.width, self.height = width, height
        
        self.gr = gui.gr(1, **options)


        coptions = dict(width=width, height=height,
                        scrollregion=(0, 0, 3*width, height),
                        bg='white', bd=2, relief=RIDGE, sticky=E+W)

        t1 = GraphTransform(height, [-3.1, 3.1])
        t2 = GraphTransform(height, [-0.1, 1.1])

        self.c1 = gui.ca(transforms=[t1], **coptions)
        self.c2 = gui.ca(transforms=[t2], **coptions)

        self.xb = gui.sb(command=self.xview, orient=HORIZONTAL,
                          sticky=E+W)

        self.c1.configure(xscrollcommand=self.xb.set)
        self.c2.configure(xscrollcommand=self.xb.set)
        gui.endgr()

        self.c1.axes(range(-3, 4))
        self.c2.axes(range(2))

    def xview(self, *args):
        self.c1.xview(*args)
        self.c2.xview(*args)


class GraphCanvas(GuiCanvas):
    def __init__(self, *args, **kwds):
        GuiCanvas.__init__(self, *args, **kwds)
        self.lines = []
        
    def axes(self, ys):
        options = dict(fill='gray50')
        dash = {True:'', False:'.'}
        for y in ys:
            self.line([[0,y], [3*self.width, y]], dash=dash[y==0], **options)

    def plot(self, data, **options):
        if len(data) < 2: return
        t = [x for x  in enumerate(data)]
        line = self.line(t, **options)
        self.lines.append(line)

    def clear(self):
        for line in self.lines:
            self.delete(line)
        

class World(Gui):

    def __init__(self):
        Gui.__init__(self)
        self.exists = True
        self.delay = 0.05
        self.setup()

    def setup(self):
        self.ca_width = 750
        self.ca_height = 300

        self.canvas = ScrollableCanvas(self, self.ca_width, self.ca_height)

        self.fr(TOP)
        self.bu(LEFT, text='Run', command=self.run)
        self.bu(LEFT, text='Stop', command=self.stop)
        self.bu(LEFT, text='Step', command=self.step)
        self.bu(RIGHT, text='Quit', command=self.quit)
        self.endfr()

    def ca(self, *args, **options):
        underride(options, fill=BOTH, expand=1)
        return self.widget(GraphCanvas, *args, **options)

        
    def quit(self):
        self.exists = False
        Gui.quit(self)

    def step(self):
        t = self.t
        i = self.i
        if i == len(t): return
        
        x = t[i]
        #print i, x
        
        self.i += 1
        
        #series.add(x)
        #series.process()

        gk = calc_gk(t, i, known_mu, known_sigma)
        #print i, x, gk

        series.make_cumulants()

        p = series.cump
        pp = series.cumpp
        winp = series.winp

        canvas = self.canvas
        canvas.c1.clear()
        canvas.c2.clear()

        canvas.c1.plot(t[:i+1], fill='red')
        canvas.c2.plot(p, fill='blue')
        canvas.c2.plot(pp, fill='green')
        #canvas.c2.plot(winp, fill='red')

        self.update()
        time.sleep(self.delay)
        
    def run(self):
        self.running = 1
        while self.exists and self.running:
            self.step()
            self.update()
            time.sleep(self.delay)

    def stop(self):
        self.running = 0


def display_series(t):
    global series
    world = World()
    world.t = t
    world.i = 0
    series = Series(fprior)
    world.series = series
    #world.run()
    world.mainloop()
    return series


def test_series(t, xs=None, prefix='', indices=[]):
    series = Series(fprior)

    for i, x in enumerate(t):
        series.add(x)
        series.process()
        n = len(series.data)
        p = series.estimates[n-1].p
        gk = sum(p)
        #print i, x, gk

        if i in indices:
            series.make_cumulants()
            if xs==None:
                dump_list(series.cump, prefix, i)
            else:
                dump_lists(xs, series.cump, prefix, i)

    return series
            

def calc_gk_ipp(t, k, mu0, sigma):
    """compute the probability that we have seen at least one
    changepoint in the first k data points in t.

    WARNING: before this function is called, series has to be
    initialized, and then this function has to be called
    consecutively with values of k from 0 to len(t)
    """
    x = t[k]
    series.add(x)
    series.process()
    p = series.estimates[k].p
    gk = sum(p)
    return gk


def cumulant(t):
    """return the cumulant sum of t (backward)
    """
    s = t[:]
    sum = 0
    for i in range(len(t)-1,-1,-1):
        sum += t[i]
        s[i] = sum
    return s
        
def calc_gk_glr0(t, k, mu0, sigma):
    """calculate gk with numin=0.0, using Equation 2.4.40 on
    page 53 of Basseville and Nikiforov
    """
    n = k+1
    t1 = t[:n]
    t2 = [y-mu0 for y in t1]
    t3 = cumulant(t2)

    t4 = t3[:]
    for j in range(n):
        t4[j] = t3[j]**2 / (k-j+1)

    gk = max(t4) / 2.0 / sigma / sigma
    return gk
    


def calc_gk_glr(t, k, mu0, sigma):
    """calculate gk with numin (global), using Equation 2.4.38 on
    page 53 of Basseville and Nikiforov
    """
    n = k+1
    t1 = t[:n]
    
    # compute the signum of the sum of centered observations
    co = cumulant([y-mu0 for y in t1])
    s1 = [signum(x) for x in co]

    # compute estimated nu for all j
    t2 = cumulant([abs(y-mu0) - numin for y in t1])
    t3 = [t2[j] / (k-j+1) for j in range(n)]
    t4 = [max(x, 0) + numin for x in t3]

    # sign of nuhat[j] should be the same as s1[j]
    s2 = [signum(x) for x in t4]
    nuhat = [s1[j] * abs(t4[j]) for j in range(n)]

    t5 = [(co[j] - nuhat[j] / 2.0 * (k-j+1)) * nuhat[j] for j in range(n)]

    gk = max(t5) / sigma**2
    return gk
    

def signum(x):
    if x<0: return -1
    if x>0: return 1
    return 0

series = None

def test_true_alarm(rseed, threshes, sigma=1.0, flag=False):
    """rseed is a random seed; threshes is a set of thresholds.
    Return a dictionary that maps from a
    threshold in threshes to the time the alarm triggered
    that threshold.
    """
    global series
    global known_sigma, known_mu

    # initialize series
    series = Series(fprior)
    mu0 = 0
    known_mu = mu0
    known_sigma = sigma

    #sys.stderr.write('%d\n' % (rseed,))

    # generate data
    seed(rseed)

    #mus = [mu0, 1, 1, 1]
    #n = 50
    #t = generate_data(mus, sigma, n)

    f = fprior
    n = int(ceil(expovariate(f)))
    m = 100

    t = generate_trial(n, m, sigma)
    #t = generate_data([0, 1, 1], sigma, n)

    # tba is time before alarm
    tba = {}
    for thresh in threshes:
        tba[thresh] = float('Inf')
            
    # iterate through the data until all the thresholds have
    # been reached
    for k in range(0, len(t)):
        gk = calc_gk(t, k, mu0, sigma)

        #sys.stderr.write('%d\t%d\t%f\t%f\t%d\n' %
        #                 (rseed, k, t[k], gk, len(threshes)))

        # we have to make a copy of threshes so we can modify
        # the original inside the loop
        for thresh in set(threshes):
            if gk > thresh:
                # print k, gk, thresh
                tba[thresh] = k - n
                threshes.remove(thresh)

        if not threshes: break

    #if threshes:
    #    sys.stderr.write('%d\t%d\t%d\n' % (rseed, k, len(threshes)))

    return tba


def alarm_time(threshes, runs=100, sigma=1.0):
    """run test_true_alarm the given number of times with
    the given sigma.
    Returns complete, farate, mat (see below).
    """

    # tbas maps from a threshold to a list of "time before alarm"
    tbas = {}
    for thresh in threshes:
        tbas[thresh] = []

    for i in range(runs):
        d = test_true_alarm(i, set(threshes), sigma)

        # d maps from a threshold to a single time before alarm
        for thresh, i in d.iteritems():
            tbas[thresh].append(i)

    complete = {}     # fraction of runs that generated an alarm
    farate = {}       # fraction of runs that had a false alarm
    mat = {}          # mean alarm time of runs that didn't fa

    for thresh in threshes:
        ncomp = len([tba for tba in tbas[thresh] if tba<1e300])
        complete[thresh] = 1.0 * ncomp / runs

        # count the number of false alarms and make a list
        # of time before alarm for runs that didn't fa
        fa = 0
        t = []
        for k in tbas[thresh]:
            if k <= 0:
                fa += 1
            else:
                t.append(k)

        # compute farate and mat
        farate[thresh] = 1.0 * fa / runs
        if len(t):
            mat[thresh] = trimmed_average(t)
        else:
            mat[thresh] = None

    return complete, farate, mat


def mta_vs_thresh(threshes, runs=100, sigma=1.0):
    complete, farate, mat = alarm_time(threshes, runs, sigma)

    for thresh in threshes:
        print thresh, complete[thresh], farate[thresh], mat[thresh]


def mta_vs_sigma(threshes, runs=100):
    """run alarm_times for a range of sigmas.
    For each sigma, interpolate the mean time until alarm
    the corresponds to a false alarm rate of fa
    """
    fas = [0.05, 0.10]
    sigmas = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]

    # matfa maps from sigma to mean alarm time, given fa
    matfas = {}
    
    for sigma in sigmas:
        complete, farate, mat = alarm_time(threshes, runs, sigma)
        #matfas[sigma] = interpolate(fa, complete, farate, mat)

        print sigma
        for thresh in threshes:
            print thresh, complete[thresh], farate[thresh], mat[thresh]

        for fa in fas:
            matfa = interpolate(fa, complete, farate, mat)
            print 'fa', fa, sigma, matfa 

    #for sigma in sigmas:
    #    print fa, sigma, matfas[sigma]


def interpolate(fa, complete, farate, mat):
    """return the mean alarm time that corresponds to the given
    false alarm rate (fa).
    """
    t = [(v, k) for k, v in farate.iteritems()]
    t.sort()
    
    from bisect import bisect
    i = bisect(t, (fa, 0.0))
    fa1, thresh1 = t[i-1]
    fa2, thresh2 = t[i]

    # if there were incomplete runs for the relevant
    # thresholds, then the data are unreliable
    #c1 = complete[thresh1]
    #c2 = complete[thresh2]
    #if c1 < 1.0 or c2 < 1.0:
    #    return None

    mat1 = mat[thresh1]
    mat2 = mat[thresh2]
    frac = (fa-fa1) / (fa2-fa1)
    res = mat1 + frac * (mat2-mat1)
    return res


def average(t):
    # compute the average of the elements in t
    return 1.0 * sum(t) / len(t)


def trimmed_average(t, percent=0.05):
    # compute the average of the elements in t, excluding the
    # top percent
    t = t[:]
    t.sort()
    trim = int(round(percent * len(t)))
    trim = max(trim, 1)
    t = t[trim:-trim]
    try:
        return 1.0 * sum(t) / len(t)
    except ZeroDivisionError:
        return None

numin = None
calc_gk = None
mu_known = False
sigma_known = False
known_sigma = None
known_mu = None
simple = False
pooled = False
fudge = None

def main(script, *args):
    #decay_problem()
    #test_dirichlet()
    #test_chi2()
    global calc_gk, numin, mu_known, sigma_known, numin, iters
    global known_mu, known_sigma, fprior, simple, pooled, fudge

    longopts = ['numin=', 'alg=', 'runs=', 'simple', 'trial', 'matsig', 'f=',
                'example', 'nile', 'pooled', 'fudge=', 'var', 'rseed=',
                'mu=', 'sigma=', 'iters=']
    
    (opts, args) = getopt.getopt(args, 'dms', longopts)
    optdict = {}
    for (key, val) in opts:
        optdict[key] = val

    try:
        fprior = float(optdict['--f'])
    except:
        fprior = 0.02

    try:
        known_mu = float(optdict['--mu'])
    except:
        known_mu = 0.0

    try:
        known_sigma = float(optdict['--sigma'])
    except:
        known_sigma = 1.0

    try:
        fudge = float(optdict['--fudge'])
    except:
        fudge = 1.0

    try:
        runs = int(optdict['--runs'])
    except:
        runs = 1000

    try:
        iters = int(optdict['--iters'])
    except:
        iters = 1

    try:
        rseed = int(optdict['--rseed'])
        seed(rseed)
    except:
        pass

    if '--simple' in optdict:
        simple = True

    if '--pooled' in optdict:
        pooled = True

    if '-m' in optdict:
        mu_known = True

    if '-s' in optdict:
        sigma_known = True

    try:
        numin = float(optdict['--numin'])
    except:
        numin = 0.5

    try:
        alg = optdict['--alg']
    except KeyError:
        alg = 'glr'

    name = 'calc_gk_' + alg
    calc_gk = eval(name)

    if alg == 'glr0' or alg == 'glr':
        mu_known = True
        sigma_known = True
        threshes = [4.8 + 0.2*i for i in range(20)]
    else:
        #threshes = [0.4 + 0.02*i for i in range(28)]
        threshes = [0.8 + 0.01*i for i in range(20)]

    if '--var' in optdict:
        mu = 0.0
        known_mu = 0.0
        t = generate_var_data(mu=0.0, sigmas=[1.0, 0.5])
        dump_list(t, 'var.', 0)

        if '-d' in optdict:
            series = display_series(t)
        else:
            series = test_series(t)

        series.make_cumulants()
        dump_list(series.cump, 'var.p', 0)
        dump_list(series.cumpp, 'var.pp', 0)
        
    elif '--nile' in optdict:
        t = read_data('nile.dat')
        xs = range(1871, 1971)
        t = [x/10 for x in t]
        dump_lists(xs, t, 'nile.', 0)

        if '-d' in optdict:        
            t = [x/100 for x in t]
            display_series(t)
        else:
            test_series(t, xs=xs, indices=[33, 66, 99], prefix='p.nile.')

    elif '--example' in optdict:
        mus = [-0.5, 0.5, 0.0]
        sigma = 1.0
        n = 50
        t = generate_data(mus, sigma, n)

        if '-d' in optdict:  
            display_series(t)
        else:
            series = test_series(t)
            
    elif '--trial' in optdict:
        mus = [0.0, 1.0]
        sigma = 0.7
        f = 0.02
        n = int(ceil(expovariate(f)))
        m = 50

        n, m = 20, 20
        t = generate_trial(n, m, sigma)

        if '-d' in optdict:  
            display_series(t)
        else:
            series = test_series(t)
            
      

    elif '--matsig' in optdict:
        mta_vs_sigma(threshes, runs)
    else:
        mta_vs_thresh(threshes, runs)


if __name__ == '__main__':
    try:
        import psyco
        psyco.full()
    except ImportError:
        pass

    profile = 0
    if profile:
        import profile
        profile.run('main(*sys.argv)')
    else:
        main(*sys.argv)

