#!/usr/bin/env python
import sys
import util
import bio1, bio2, bio3,bio4, bio5, bio6
import graph


import cProfile, pstats, io


with open(sys.argv[1], 'r') as f:
    #dna = f.readline().strip()
    sequence = [int(n) for n in f.readline().strip().split(' ')]
    res = sorted(bio2.branch_bound_1(sequence), reverse=True)
    #print(res)
    print(util.pp(res, join_char=' '))
