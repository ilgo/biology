#!/usr/bin/env python
import sys
from itertools import permutations, combinations_with_replacement, product
from bio1 import pp
from bio2 import cyclic_cutter

with open(sys.argv[1], 'r') as f:
    dna = f.readline().strip()
    #length = int(f.readline().strip())
#
#res = [n+1 for n, txt in enumerate(cyclic_cutter(dna, len(pattern))) if txt == pattern]    
#print(pp(res))

#n = int(sys.argv[1])
#
#c = 0
#res = []
#for perm in permutations(range(1, n+1)):
#    c +=1
#    res.append(pp(perm))
#print(c)
#print('\n'.join(res))

#t = {
#    'A' :   71.03711,
#    'C' :   103.00919,
#    'D' :   115.02694,
#    'E' :   129.04259,
#    'F' :   147.06841,
#    'G' :   57.02146,
#    'H' :   137.05891,
#    'I' :   113.08406,
#    'K' :   128.09496,
#    'L' :   113.08406,
#    'M' :   131.04049,
#    'N' :   114.04293,
#    'P' :   97.05276,
#    'Q' :   128.05858,
#    'R' :   156.10111,
#    'S' :   87.03203,
#    'T' :   101.04768,
#    'V' :   99.06841,
#    'W' :   186.07931,
#    'Y' :   163.06333 ,
#}
#
#aa = sys.argv[1]
#res = sum([t[c] for c in aa])
#print("%0.3f" %res)

#chars.insert(0, '')
#res = [''.join(p) for p in product(*[chars for _ in range(length)])]
#
#print('\n'.join(res))

res = [dna.count(c) for c in 'ACGT']
print(pp(res))
