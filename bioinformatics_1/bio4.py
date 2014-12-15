#!/usr/bin/env python

import argparse
import math
import util
from bio2 import cyclic_cutter
from itertools import product
from collections import defaultdict
from graph import GraphFactory, euler_cycle, is_balanced

from Bio import SeqIO, Alphabet

concat = lambda ps: ','.join([p for p in ps]) if isinstance(ps, list) else ps 


def str_composition(text, n):
    '''
    get  all substrings of length n

    retuen a list of substrings len(n)

    :param text: the string to analyse
    :type text: `str`
    :param n: lenght of substrings
    :type n: `int`

    :rtype: [`str`, ...]
    '''
    pos = 0
    while pos < len(text) - n + 1:
        yield text[pos:pos+n]
        pos += 1


def format_overlaps(adjacency_dict):
    '''
    format an overlap dictionary for Rosalind

    returns a formateted list of lexicographically sorted overlaps. eg:'ATCGG' -> 'TCGGT'

    :param adjacency_dict: a dict of overlaps
    :type adjacency_dict: {'source:'target', ...}
    
    :rtype: `str`
    '''
    out = ['%s -> %s' %(src, concat(adjacency_dict[src])) for src in sorted(adjacency_dict.keys())]
    return '\n'.join(out)


def str_overlap(kmers):
    '''
    create an adjacency graph for the given overlapping kmers

    return a dict of overlapping kmers. eg:'ATCGG' -> 'TCGGT'
    
    :param kmers: the overlapping kmers
    :type kmers: [`str`, ...]

    :rtype: {'source:'target', ...}
    '''
    overlaps = {}
    for kmer in kmers:
        suffix = kmer[1:]
        if suffix not in overlaps:
            overlaps[suffix] = [kmer, []]   
        else:
            overlaps[suffix][0] = kmer

        prefix = kmer[:-1]
        if prefix not in overlaps:
            overlaps[prefix] = ['', [kmer]]
        else:
            overlaps[prefix][1].append(kmer)
                
    return {k:v for k, v in overlaps.values() if k and v}


def deBruijn(kmers, pairs=False):
    '''
    constructs a deBruijn graph from a list of kmers

    returns an adjacency deBruijn graph 

    :param kmers: an iterable returning kmers
    :type k: `iter`
    :param pairs: if we got read-pairs
    :type pairs: `bool`

    :rtype: {'source: ['target'...], ...}
    '''
    overlaps = {}
    for kmer in kmers:
        if pairs:
            suffix = (kmer[0][1:],  kmer[1][1:])
            prefix = (kmer[0][:-1], kmer[1][:-1])
        else:
            suffix = kmer[1:]
            prefix = kmer[:-1]
        if prefix not in overlaps:
            overlaps[prefix] = [suffix]   
        else:
            overlaps[prefix].append(suffix)
    #return {k:sorted(list(set(v))) for k, v in overlaps.items()}
    return overlaps


def eulerCycle(edge_iter, start_node=0):
    '''
    find an euler cycle from a directed graph

    returns a list of nodes that form an euler cycle

    :param edge_iter: an iterable for edges
    :type edge_iter: `iter` -> (`int, [`int`, ...])
    :param start_node: for debugging purpose the node to start the tour
    :type start_node: `int`

    :rtype: [`int`, ...]
    '''
    graph = GraphFactory(edge_iter, directed=True)
    if not is_balanced(graph):
        return []
    return euler_cycle(graph, start_node)
    

def eulerPath(edge_iter):
    '''
    find an euler path in a directed graph

    returns a list of nodes that form an euler path

    :param edge_iter: an iterable for edges
    :type edge_iter: `iter` -> (`int, [`int`, ...])

    :rtype: [`int`, ...]
    '''
    graph = GraphFactory(edge_iter, directed=True)
    in_outs = {}
    for k,vs in graph.items():
        entry = in_outs[k] if k in in_outs else [0,0]
        entry[1] += len(vs)
        in_outs[k] = entry
        for v in vs:
            entry2 = in_outs[v] if v in in_outs else [0,0]
            entry2[0] += 1
            in_outs[v] = entry2
   
    start, end = None, None
    for k, inout in in_outs.items():
        #diff = inout[0] - inout[1]
        #if diff < 0:
        #    start = k
        #elif diff > 0:
        #    end = k
        
        if inout[0] == 0:
            start = k
        elif inout[1] == 0:
            end = k
        if start and end:
            break
    graph[end].append(start)
    path = euler_cycle(graph)
    start_idx = path.index(end) + 1
    return path[start_idx:-1] + path[:start_idx]
    

def str_reconstruct(kmers, pairs=False, distance=0):
    '''

    '''
    edge_iter = ((k,v) for k, v in deBruijn(kmers, pairs=pairs).items())
    kmers = eulerPath(edge_iter)
    #print(kmers)
    if pairs:
        s1 = kmers[0][0][:-1] + ''.join([kmer[0][-1] for kmer in kmers])
        s2 = kmers[0][1][:-1] + ''.join([kmer[1][-1] for kmer in kmers])
        return s1 + s2[-(len(kmers[0][0])+distance)-1:]
    else:
        return kmers[0][:-1] + ''.join([kmer[-1] for kmer in kmers])
   

def complete(kmers, f):
    '''
    try to get a more complete debrujin graph by cutting exisitng kmers into smaller pieces

    :param kmers: the incomplete list of kmers 
    :type kmers: [`str`, ...]
    :param `int` f: ratio to cut kmers
    :returns: [`str`, ...]
    '''
    small_kmers = []
    lenght = math.ceil(len(kmers[0]) / int(f))
    for kmer in kmers:
        small_kmers.extend([k for k in util.kmer_gen(kmer, lenght)])
    return small_kmers

def contigs(kmers):
    '''
    generate the contigs for a set of kmers

    :param [`str`, ...] kmers: a list of kmers
    :rtype: [`str`, ...]
    :returns: a list of contigs
    '''
    kmers = complete(kmers, 2)
    overlaps = {}
    cons = []

    for kmer in kmers:
        suffix = kmer[1:]
        prefix = kmer[:-1]
        if prefix not in overlaps:
            overlaps[prefix] = [suffix]   
        else:
            overlaps[prefix].append(suffix)
    keys = set(overlaps.keys())
    for vs in overlaps.values():
        keys.update(vs)

    incoming = {k:0 for k in keys}
    outgoing = {k:0 for k in keys}
    for k, vs in overlaps.items():
        outgoing[k] = len(vs)        
        for v in vs:
            incoming[v] +=1

    starts = [k for k,vs in overlaps.items() for _ in vs]
    for start in starts:
        val = incoming[start]
        if val == 1:
            continue
        con = [start]
        ends = overlaps.get(start, None)
        while ends and len(ends) > 0:
            end = ends.pop()
            con.append(end[-1])
            if incoming[end] != 1:
                break
            start = end
            ends = overlaps.get(start, None)
        cons.append(''.join(con))    
    
    return sorted(cons)
    


def universal_circular(n):
    '''
    
    '''
    kmers = [''.join(i) for i in product('01', repeat=n)]
    edge_iter = ((k,v) for k, v in deBruijn(kmers).items())
    graph = GraphFactory(edge_iter, directed=True)
    kmers = euler_cycle(graph, start_node=kmers[0][:-1])
    #print(kmers)
    return kmers[0][:-1] + ''.join([kmer[-1] for kmer in kmers[:-n+1]])


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', required=True, help='fasta file')
    parser.add_argument('-k', type=int, required=True, help='kmer size')
    parser.add_argument('-d', required=False, help='pair distance')

    args = parser.parse_args()

    seqio = SeqIO.parse(args.f, 'fasta', Alphabet.DNAAlphabet())
    seq = next(seqio).seq
    reads = [str(s) for s in str_composition(seq, args.k)]
    
    restruct = str_reconstruct(reads)
    assert restruct == str(seq)
    print('done')
