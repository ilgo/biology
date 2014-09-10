#!/usr/bin/env python
import sys
import bio1, bio2

with open(sys.argv[1], 'r') as f:

    pep = f.readline().strip()
    #amino = f.readline().strip()
    masses = [int(p) for p in pep.split(' ')]

    res = bio2.branch_bound_1(masses)
    print(bio1.pp(res))
