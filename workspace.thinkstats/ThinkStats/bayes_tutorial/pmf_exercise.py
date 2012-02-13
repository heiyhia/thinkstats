import Pmf

def PrintPmf(pmf):
    for value, prob in pmf.Items():
        print value, prob

def ProbBigger(pmf1, pmf2):
    new = Pmf.Pmf()
    return new

six_sided_die = Pmf.MakePmfFromList([1, 2, 3, 4, 5, 6])

PrintPmf(six_sided_die)
