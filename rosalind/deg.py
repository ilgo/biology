#!/usr/bin/env python
import sys
from graph import GraphFactory, pp

fp = sys.argv[1]

g = GraphFactory(fp)

res = []
for s in sorted(g.keys()):
    res.append(sum([len(g[t]) for t in g[s]]))

print(pp(res))

