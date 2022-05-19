import numpy as np
import collections
import itertools
def pf(k):
    i = 2
    while i * i <= k:
        if k % i == 0:
            k /= i
            yield i
        else:
            i += 1
    if k > 1:
        yield k

def product(s):
    result = 1
    for i in s:
        result *= i
    return result

def get_divisors(k):
    factors = pf(k)
    factors = collections.Counter(factors)
    powers = [[factor**i for i in range(count + 1)] for factor, count in factors.items()]
    for combs in itertools.product(*powers):
        yield product(combs)

N = 128
Lambda = 24
N = int(input("How many processors?\n"))
Lambda = int(input("number of lattice points per dimension in cubic lattice?\n"))

divisors = list(get_divisors(Lambda))
divisors = [i-1  for i in divisors]
combinations = itertools.combinations_with_replacement(divisors,4)
result = [[x, (x[0]+1)*(x[1]+1)*(x[2]+1)*(x[3]+1)] for x in combinations if (x[0]+1)*(x[1]+1)*(x[2]+1)*(x[3]+1)<=N]
sorted_list = sorted(result,key = lambda x: x[1])
for i in range(len(sorted_list)):
    print(sorted_list[i])
