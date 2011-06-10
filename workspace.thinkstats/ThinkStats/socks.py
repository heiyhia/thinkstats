"""A quick simulation in response to a question on reddit:
http://www.reddit.com/r/statistics/comments/hszzd/get_me_on_the_right_track/

There are 12 socks in a drawer, 10 white, 2 red. I draw 5 socks, what
is the probability that I get 2 reds?

My assumption was that i could get the probability of all white, which
is pretty easy (10/12)(9/11)(8/10)(7/9)(6/8) = 7/22 and add the
probability of 1 red. But I couldn't even figure out how to get that
probability.

"""

from random import shuffle

n = 10000
count = 0.0
t = list('w' * 10 + 'r' * 2)

for i in xrange(n):
    shuffle(t)
    if t[:5].count('r') == 2:
        count += 1

print count / n
