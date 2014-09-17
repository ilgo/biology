#!/usr/bin/env python
import sys, math
from bio1 import pp
from collections import defaultdict


with open(sys.argv[1], 'r') as f:
    k, n = [int(i) for i in f.readline().strip().split()]
    half = math.floor(n/2)
    res = []
    for i in range(k):
        a = (int(a) for a in f.readline().strip().split(' '))
        data = defaultdict(int)
        done = False
        for n in a:
            data[n] += 1
            if data[n] > half:
                res.append(n)
                done = True
                break
        if not done:
            res.append(-1)
print(pp(res)) 
