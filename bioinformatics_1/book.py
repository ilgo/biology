#!/usr/bin/env python
import sys
import util
import turnpike
import bio1, bio2, bio3,bio4, bio5, bio6, bio7
import graph

from datetime import datetime
from collections import defaultdict
import cProfile, pstats, io


k = 3 
t = 5
dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
res = bio3.greedy_motif(k, t, dnas, pseudo=True)
print( res)

sys.exit()

#now = datetime.now()

f_name = 'data/test.txt' if len(sys.argv) == 1 else sys.argv[1]
with open(f_name, 'r') as f:
   
    #k = int(f.readline().strip())
    dnas = [line.strip() for line in f.readlines()] 
    print(dnas)
    res = bio4.contigs(dnas)
    #print(res)
    #res = [bio6.signed_permutation_string(seq) for seq in seqs]
    print(util.pp(res, join_char='\n'))

    #print(datetime.now()-now)

    #src = int(f.readline().strip())
    #sink = int(f.readline().strip())
    #edges = [bio5.make_edge(l.strip()) for l in f.readlines()]
    #dag = bio5.DAG(sink, edges)
    #res = bio5.longest_path(src, sink, dag)
    #res = bio5.longest_path_fmt(res)
    #print(res)
    #print(util.pp(res, join_char='\n'))

