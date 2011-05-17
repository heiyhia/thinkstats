import random

def Sample(n=2500):
    return [random.random() for i in xrange(n)]

def Count(sample, p):
    return len([x for x in sample if x<p])

def Test(n=1250):
    s1 = Sample(n)
    f1 = 1.0 * Count(s1, 0.07) / n
    s2 = Sample(n)
    f2 = 1.0 * Count(s2, 0.07) / n
    return f1 - f2
    
count = 0
for i in range(10000):
    diff = Test()
    if diff > 0.033:
        count += 1

print count
