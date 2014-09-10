#!/usr/bin/env python
import sys
import re
from random import choice
from numpy import array

def seq(n, alphabeth=['A', 'T', 'C', 'G']):
    '''
    Random sequence of N chars from the set(['A', 'T', 'C', 'G'])

    :param n: the size of the sequence
    :type n: `int`
    :param alphabeth: the letters in the sequence
    :type alphabeth: [`str`, ... ]

    :rtype: `str`
    '''
    seq = (choice(alphabeth) for _ in range(n))
    return ''.join(seq)


def find_all_re(segment, gene):
    '''
    Find the segment in the gene with a regex
    
    returns a list of positions where the segmen is found in the gene

    :param segment: the substring to be found
    :type segment: `str`
    :param gene: the sequence to be searched
    :type gene: `str`

    :rtype: [`int`, ... ]
    '''
    pat = re.compile(segment)
    _iter = pat.finditer(gene)
    return [m.start() for m in _iter]

