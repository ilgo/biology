#!/usr/bin/env python
import sys
import util
import turnpike
import bio1, bio2, bio3,bio4, bio5, bio6, bio7
import graph

from datetime import datetime
from collections import defaultdict
import cProfile, pstats, io



#now = datetime.now()

f_name = 'data/test.txt' if len(sys.argv) == 1 else sys.argv[1]
with open(f_name, 'r') as f:
   
     
    s1 = f.readline().strip()
    s2 = f.readline().strip()
    pathLen, path = bio5.overlap_align(s1, s2)
    print(pathLen, path)
    res = bio5.align_strings(s1, s2, path)
    res.insert(0, pathLen)

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

