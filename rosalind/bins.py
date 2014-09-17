#!/usr/bin/env python
import sys, bisect
from bio1 import pp
from graph import GraphFactory, connected_components

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect.bisect_left(a, x)
    j = bisect.bisect_right(a, x)
    return -1 if i==j else i+1
    #print(a, x, i, j)
    #try:
    #    return a[i:j].index(x) + i
    #except:
    #    return -1

with open(sys.argv[1], 'r') as f:
    n = int(f.readline().strip())
    m = int(f.readline().strip())
    A = [int(a) for a in f.readline().strip().split(' ')]
    B = [int(a) for a in f.readline().strip().split(' ')]
   
assert n == len(A)
assert m == len(B)

res = [index(A,b) for b in B] 
print(pp(res))
