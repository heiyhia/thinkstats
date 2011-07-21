"""This file contains solutions to exercises in "Think Stats",
by Allen B. Downey, available from greenteapress.com

Copyright 2010 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html


"""

"""
If I roll two dice and the total is 8, what is the chance that
one of the dice is a 6?

Solution: there are 5 ways to make 8 (2-6, 3-5, 4-4, 5-3, 6-2)
and two of them have a 6, so the answer is 2/5.
"""

"""
If I roll 100 dice, what is the chance of getting all sixes?
"""
print (1.0/6)**100

"""
What is the chance of getting no sixes?
"""
print (5.0/6)**100


"""
The following questions are adapted from Mlodinow, The Drunkard's Walk.

* If a family has two children, what is the chance that they
  have two girls?

1/4

* If a family has two children and we know that at least one of
  them is a girl, what is the chance that they have two girls?

1/3

* If a family has two children and we know that the older one is a
  girl, what is the chance that they have two girls?

1/2

* If a family has two children and we know that at least one of
  them is a girl named Florida, what is the chance that they have
  two girls?

1/2
"""

"""
Write a program that simulates the Monty Hall problem and use
it to estimate the probability of winning if you stick and if
you switch.
"""



"""
To understand the Monty Hall problem, it is important to realize
that by deciding which door to open, Monty is giving you information.
To see why this matters, imagine the case where Monty doesn't
know where the prizes are, so he chooses Door B or C at random.

If he opens the door with the car, the game is over, you lose, and
you don't get to choose whether to switch or stick.

Otherwise, are you better off switching or sticking?

Solution: In this case, Monty is contributing no information, so
the chance of winning is 1/2 whether you stick or switch.

"""


"""
If you go to a dance where partners are paired up randomly, what
percentage of opposite sex couples will you see where the woman is
taller than the man?

In the BRFSS (see Section~\ref{lognormal}), the distribution of
heights is roughly normal with parameters \mymu~=~178 cm and
\sigmasq~=~59.4 cm for men, and \mymu~=~163 cm and \sigmasq~=~52.8 cm for
women.
"""
import random
import math
count = 0.0
for i in range(1000):
    male_height = random.normalvariate(178, math.sqrt(59.4))
    female_height = random.normalvariate(163, math.sqrt(52.8))
    if female_height > male_height:
        count += 1
print count / 1000

"""
If I roll two dice, what is the chance of rolling at least one 6?
\index{dice}
"""
print 1.0/6 + 1.0/6 - 1.0/36

"""
What is the general formula for the probability of \A~or \B~but not both?

Solution: P(A XOR B) = P(A) + P(B) - 2 P(A \AND B)
"""

"""
If you flip a coin 100 times, you expect about 50 heads, but what
is the probability of getting exactly 50 heads?

Solution: 0.079589
"""

"""
This exercise is from http://wikipedia.org/wiki/Bayesian_inference

``Suppose there are two full bowls of cookies. Bowl 1 has 10 chocolate
  chip and 30 plain cookies, while Bowl 2 has 20 of each. Our friend
  Fred picks a bowl at random, and then picks a cookie at random. The
  cookie turns out to be a plain one. How probable is it that Fred
  picked it out of Bowl 1?''

"""

"""
The blue M&M was introduced in 1995.  Before then, the color mix in
a bag of plain M&Ms was (30% Brown, 20% Yellow, 20% Red, 10%
Green, 10% Orange, 10% Tan).  Afterward it was (24% Blue , 20%
Green, 16% Orange, 14% Yellow, 13% Red, 13% Brown).


A friend of mine has two bags of M&Ms, and he tells me
that one is from 1994 and one from 1996.  He won't tell me which is
which, but he gives me one M&M from each bag.  One is yellow and
one is green.  What is the probability that the yellow M&M came
from the 1994 bag?
"""

"""
This exercise is adapted from MacKay, Information
  Theory, Inference, and Learning Algorithms:

Elvis Presley had a twin brother who died at birth.  According to the
Wikipedia article on twins:

``Twins are estimated to be approximately 1.9% of the world population,
with monozygotic twins making up 0.2% of the total---and 8% of all
twins.''

What is the probability that Elvis was an identical twin?
"""
