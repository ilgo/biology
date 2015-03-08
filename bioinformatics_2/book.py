#!/usr/bin/env python
import sys
import util
import bio2_1

from datetime import datetime
from collections import defaultdict



f_name = 'data/test.txt' if len(sys.argv) == 1 else sys.argv[1]
with open(f_name, 'r') as f:
  
    #k = int(f.readline().strip())
    s1 = f.readline().strip()
    s2 = f.readline().strip()
    #perm = [int(k) for k in f.readline().split(' ')]
    #res = bio5.middle_edge(s1, s2)
    #print(res)
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

