#!/usr/bin/env python
import sys
from graph import GraphFactory, connected_components

   
graph = GraphFactor(sys.argv[1]) 
cc = connected_components(graph)
print(cc)
