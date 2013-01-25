import Pmf

def PrintPmf(pmf):
    for value, prob in pmf.Items():
        print value, prob

def ProbBigger(pmf1, pmf2):
    total = 0.0
    return total

six_sided_die = Pmf.MakePmfFromList([1, 2, 3, 4, 5, 6])

PrintPmf(six_sided_die)
