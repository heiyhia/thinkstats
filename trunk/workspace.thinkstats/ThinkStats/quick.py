import random

def Sample(n):
    return [random.random() for i in xrange(n)]

def Fraction(sample, p):
    success = len([x for x in sample if x<p])
    return float(success) / len(sample)

def Test(n, p):
    s1 = Sample(n)
    f1 = Fraction(s1, p)
    s2 = Sample(n)
    f2 = Fraction(s2, p)
    return f1 - f2
    
count = 0
for i in range(10000):
    diff = Test(1250, 0.07)
    if diff > 0.033:
        count += 1

print count
