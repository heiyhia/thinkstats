import descriptive
import Pmf
import matplotlib.pyplot as pyplot
import plot

def ConditionPmf(pmf, filter_func):
    """Computes a conditional PMF based on a filter function.
    
    Args:
        pmf: Pmf object
        
        filter_func: function that takes a value from the Pmf and returns
                     a boolean
        
    Returns:
        new Pmf object
    """
    d = pmf.GetDict()
    # make a copy before we mangle it
    d = dict(d)
    for x in d:
        if filter_func(x):
            d[x] = 0

    cond_pmf = Pmf.Pmf(d)
    cond_pmf.Normalize()
    return cond_pmf


def ConditionOnWeeks(pmf, week=39):
    """Computes a PMF conditioned on the given number of weeks.
    
    Args:
        pmf: Pmf object
        
        weeks: the current duration of the pregnancy
        
    Returns:
        new Pmf object
    """
    def filter_func(x):
        return x < week

    cond = ConditionPmf(pmf, filter_func)
    return cond


def main():
    MakeFigure()
    
def MakeFigure():
    pool, firsts, others = descriptive.MakeTables()

    weeks= range(35, 46)
    
    # probs is a map from table name to list of conditional probabilities
    probs = {}
    for table in [firsts, others]:
        name = table.pmf.name
        probs[name] = []
        for week in weeks:
            cond = ConditionOnWeeks(table.pmf, week)
            prob = cond.Prob(week)
            print week, prob, table.pmf.name
            probs[name].append(prob)
            
    # make a plot with one line for each table
    pyplot.clf()        
    for name, ps in probs.iteritems():
        pyplot.plot(weeks, ps, label=name)
        print name, ps
        
    plot.Plot('conditional',
              xlabel='weeks',
              ylabel=r'Prob{x $=$ weeks | x $\geq$ weeks}',
              title='Conditional Probability',
              show=True,
              )

if __name__ == "__main__":
    main()
