#!/usr/bin/env python
import sys

n = int(sys.argv[1])

f1 = 0
f2 = 1
c = 2

while c <= n:
    tmp = f1 + f2
    f1 = f2
    f2 = tmp
    c += 1
print(f2)
