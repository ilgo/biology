from itertools import product, chain, permutations, combinations, filterfalse

peptide = {
    'A' : ('GCU', 'GCC', 'GCA', 'GCG'),
 	'L' : ('UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'),
    'R' : ('CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'K' : ('AAA', 'AAG'),
    'N' : ('AAU', 'AAC'),
    'M' : ('AUG',),
    'D' : ('GAU', 'GAC'),
    'F' : ('UUU', 'UUC'),
    'C' : ('UGU', 'UGC'),
    'P' : ('CCU', 'CCC', 'CCA', 'CCG'),
    'Q' : ('CAA', 'CAG'),
    'S' : ('UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'),
    'E' : ('GAA', 'GAG'),
    'T' : ('ACU', 'ACC', 'ACA', 'ACG'),
    'G' : ('GGU', 'GGC', 'GGA', 'GGG'),
    'W' : ('UGG',),
    'H' : ('CAU', 'CAC'),
    'Y' : ('UAU', 'UAC'),
    'I' : ('AUU', 'AUC', 'AUA'),
    'V' : ('GUU', 'GUC', 'GUA', 'GUG'),
#START 	AUG 	
    'X' : ('UAA', 'UGA', 'UAG')
}

masses = {
    'A' : 71,
    'R' : 156,
    'N' : 114,
    'D' : 115,
    'C' : 103,
    'E' : 129,
    'Q' : 128,
    'G' : 57,
    'H' : 137,
    'I' : 113,
    'L' : 113,
    'K' : 128,
    'M' : 131,
    'F' : 147,
    'P' : 97,
    'S' : 87,
    'T' : 101,
    'W' : 186,
    'Y' : 163,
    'V' : 99,
}


to_dna = lambda rna : rna.replace('U', 'T')
to_rna = lambda dna : dna.replace('T', 'U')
mass = lambda pep : sum([masses[p] for p in pep])


def _rna_translate(rna, stop=True):
    '''
    translates an rna into its amino acid

    returns a list of amino-acid characters

    :param rna: a RNA string
    :type rna: `str`
    :param stop: return if a STOP codon is found
    :type stop: `bool`

    :rtype: [`str`, ...]
    '''
    amino = []
    for codon in _codon_cutter(rna):
        for k,v in peptide.items():
            if codon in v:
                if k == 'X' and stop:
                    return amino
                amino.append(k)
    return amino


def _codon_cutter(rna):
    '''
    cuts an RNA string into its triple letter amino parts

    returns a sequence of triple letters

    :param rna: a RNA string
    :type rna: `str`

    :rtype: `str`
    '''
    pos = 0
    while pos < len(rna):
        yield rna[pos:pos+3]
        pos += 3
    

def _dna_frames(dna):
    '''

    '''
    yield dna
    #yield dna[1:-2]
    #yield dna[2:-1]

def peptide_to_rna(aminos):
    '''
    translate a peptide back into the various rna represenations

    returns a list of all possible translations

    :param aminos: a peptide 
    :type aminos: `str`

    :rtype: [`str`, ...]
    '''
    return [to_dna(''.join(pep)) for pep in product(*[peptide[amino] for amino in aminos])]

def protein_translate(rna, stop=True):
    '''
    translate a RNA string to its amino-acid representation

    returns an amino-acid string

    :param rna: a RNA string
    :type rna: `str`
    :param stop: return if a STOP codon is found
    :type stop: `bool`

    :rtype: `str`
    '''
    return ''.join(_rna_translate(rna, stop=stop))


#def peptide_encode(dna, amino):
#    '''
#
#    '''
#    rna = to_rna(dna)
#    valid_codes = []
#    for pep, geno in product(*[_amino_perms(amino), _dna_frames(rna)]):
#        pep = ''.join(pep)
#        if pep in geno:
#            print(pep, geno)
#            valid_codes.append(pep)
#    print(valid_codes)
#    for pep, geno in product(*[_amino_perms(amino[::-1]), _dna_frames(rna)]):
#        pep = ''.join(pep)
#        dna = dna[::-1]
#        #print(pep, geno)
#        if pep in geno:
#            valid_codes.append(pep)
#    return valid_codes  


def cyclic_cutter(text, cycle=1):
    '''
    gets all unique and sequential pieces of lenght cycle from text

    returns a string generator

    :param text: a dna, rna, peptide....
    :type text: `str`
    :param cycle: the lenght of the pieces to cut out
    :type cycle: `int`

    :rtype: generator -> `str`
    '''
    pos = 0
    _text = text + text[:cycle-1]
    while pos < len(text):
        yield _text[pos:pos+cycle]
        pos += 1


def theoretical_spectrum(peptide):
    '''
    calculates the theoretical spectrum of a given peptide

    returns a sorted list of all the molecular masses, including 0 and the full mass

    :param peptide: a peptide string
    :type peptide: `str`
    
    :rtype: [0, `int`, ...]
    '''
    peps = chain(*[cyclic_cutter(peptide, n) for n in range(1, len(peptide))])
    pep_masses = [mass(pep) for pep in peps]
    pep_masses.insert(0,0)
    pep_masses.append(mass(peptide))
    return sorted(pep_masses)

def branch_bound_1(theoretical_masses):
    '''
    finds all peptides that are possible matches for the theoretical mass spectrum

    returns a formatted list of peptides

    :param theoretical_masses: a sorted list of peptide masses
    :type theoretical_masses: [0, `int`, ... ]

    :rtype: [ `str`, ... ]
    '''
    single_peps = [[v] for k,v in masses.items()]
    peps = [p for p in single_peps]
    while True:
        peps = [p for p in filter(lambda x: sum(x) in theoretical_masses, peps)]
        matches = [p for p in peps if sum(p) == theoretical_masses[-1]]
        if matches:
            peps = list(set(tuple(m) for m in matches))
            break
        peps = [p+sp for p in peps for sp in single_peps]
     
    return ['-'.join([str(i) for i in p]) for p in peps]
