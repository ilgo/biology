#!/usr/bin/env python
import sys
import bio1, bio2, bio4
import graph


with open(sys.argv[1], 'r') as f:

    kmers = [kmer.strip() for kmer in f.readlines()]
    res = bio4.str_reconstruct(kmers)
    print(res)
