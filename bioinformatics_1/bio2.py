from itertools import product, chain, permutations, combinations, filterfalse
from collections import defaultdict
import bio1, util

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
rev_mass = {v:k for k,v in masses.items()}


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


def peptide_encode(dna, amino):
    '''

    '''
    valid_codes = []

    amino_gen = product(*[peptide[acid] for acid in amino])
    for pep in amino_gen:
        amino_acid = to_dna(''.join([p for p in pep]))
        for kmer in util.kmer_gen(dna, len(amino_acid)): 
            if kmer == amino_acid:
                valid_codes.append(amino_acid)
        for kmer in util.kmer_gen(bio1.complement(dna), len(amino_acid)): 
            if kmer == amino_acid:
                valid_codes.append(bio1.complement(amino_acid))

    return valid_codes


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


def theoretical_spectrum(peptide, cyclic=True):
    '''
    calculates the theoretical spectrum of a given peptide

    returns a sorted list of all the molecular masses, including 0 and the full mass

    :param peptide: a peptide string
    :type peptide: `str`
    :param `boolb` cyclic: True if the theoretical spectrum is a cyclic one.
    
    :rtype: [0, `int`, ...]
    '''
    peptide_gen = cyclic_cutter if cyclic else util.kmer_gen
    peps = chain(*[peptide_gen(peptide, n) for n in range(1, len(peptide))])
    pep_masses = [mass(pep) for pep in peps]
    pep_masses.append(mass(peptide))
    pep_masses.insert(0,0)
    return sorted(pep_masses)


def cyclo_peptide_sequence(theoretical_masses):
    '''
    finds all peptides that are possible matches for the theoretical mass spectrum

    returns a formatted list of peptides

    :param theoretical_masses: a sorted list of peptide masses
    :type theoretical_masses: [0, `int`, ... ]

    :rtype: [ `str`, ... ]
    '''
    aminos = tuple(v for v in set(masses.values()))
    outer_spectrum = spectrum_dict(theoretical_masses)

    results = []
    peps = [[0]]
    while peps:
        peps = _branch(peps, aminos)
        peps, res = _bound(peps, outer_spectrum)
        results.extend(res)        
        #print(res, peps)
    
    return ['-'.join([str(i) for i in p]) for p in res]


#-----------------------------------------------------------------------
#
# functions used by cyclo_peptide_sequenc
#

def _branch(peps, aminos):
    if peps == [[0]]:
        return [[amino] for amino in aminos] 

    next_peps = []
    for pep in peps:
        for amino in aminos:
            next_pep = pep[:]
            next_pep.append(amino)
            next_peps.append(next_pep)
    return next_peps


def _bound(peps, outer_spectrum):
    next_peps, res = [],[]
    parent_mass = max(outer_spectrum.keys())

    _peps = []
    for pep in peps:
        _pep = cyclic_spectrum(pep)
        inner_spectrum = spectrum_dict(_pep)
        if is_inner_spectrum(inner_spectrum, outer_spectrum):
            if sum(pep) == parent_mass:
                res.append(pep)
            else:
                next_peps.append(pep)
    return next_peps, res


def is_inner_spectrum(inner_spec_dict, outer_spec_dict):
    '''
    determines if the inner_dict is contained within the outer dict.
    containemtn means that all masses in inner must be in outer, and no inner mass_count can be greater than the respective outer mass_count
    ''' 
    #print(inner_spec_dict)
    #print(outer_spec_dict)
    #print('>>>>>')
    for k, v in inner_spec_dict.items():
        if k not in outer_spec_dict:
            return False
        if v > outer_spec_dict[k]:
            return False
    return True


#-----------------------------------------------------------------------

def number_of_peptides_with_mass(mass, in_data=None):
    '''
    dynamic programming solution to find all possible ways peptides add up to a given mass
    
    returns the number of ways peptides can have a given mass

    :param mass: the mass to be achieved
    :type mass: `int`
    :param in_data: the masses of the components
    :type in_data: [`int, ...]

    :rtype: `int`
    '''
    out = 0
    data = set(in_data)
    sources = {k:1 for k in data}
    while sources:
        targets = {}
        for src, add in product(*[sources.keys(), data]):
            val = src + add
            targets[val] = targets.get(val, 0) + sources.get(src,1)
        if mass in targets:
            out += targets.pop(mass)
        sources = {k:v for k,v in targets.items() if k <= mass}

    return out


def cyclic_spectrum(masses, zero=False):
    '''
    calculates the  cyclopeptide sequence for the given masses.
    
    the result will not include a zero.

    :param [`int`, ...] masses: the masses to be added    
    :returns: the list of all masses in the cycle formed by the masses, including the mass
    :rtype: [`int`, ...]
    '''
    spectrum = masses[:]
    if len(masses) <= 1:
        return spectrum

    for d in range(2, len(masses)):
        cy_mass = masses + masses[:d-1]
        for p in range(len(cy_mass)-d+1):
            spectrum.append(sum(cy_mass[p:p+d]))

    spectrum.append(sum(masses))
    if zero:
        spectrum.insert(0,0)
    return spectrum
   



def spectrum_dict(masses):
    '''
    generates the spectrum dictionary. a mapping from mass to mass_count
 
    '''
    spec_dict = {}
    for m in masses:
        if m in spec_dict:
            spec_dict[m] += 1
        else:
            spec_dict[m] = 1
    return spec_dict


def leaderboard_cyclopeptide_sequence(sequence, N, allow=None, cyclic=True):

    if not allow:
        aminos = tuple(v for v in set(masses.values()))
    else:
        aminos = tuple(set(allow))
    spectrum = spectrum_dict(sequence)
    mass = max(sequence)

    best_peptide = None
    best_score = 0
    leaders = [[0]]
    while leaders:
        leaders = _branch(leaders, aminos)
        leaders, res = _bound2(leaders, mass)
        for leader in res:
            score = cyclopeptide_scoring(leader, spectrum,cyclic=cyclic)
            if score > best_score:
                best_peptide = leader
                best_score = score
        leaders = trim(leaders, spectrum, N, cyclic=cyclic)
         
    return '-'.join([str(i) for i in best_peptide])


def _bound2(leaders, mass):
    new_leaders = []
    res = []
    for leader in leaders:
        lead_mass = sum(leader)
        if lead_mass == mass:
            res.append(leader)
        elif lead_mass < mass:
            new_leaders.append(leader)
    return new_leaders, res    


def cyclopeptide_scoring(peptide, spec_dict, cyclic=True):
    '''
    what is the score of the peptide against the spectrum
    spec_dict parameter is generated with the spectrum_dict(mases) function    

    :param `str` peptide: the peptide to get the score from
    :param {`int`: `int`} spectrum: an value:count dictionary for masses and their occurrence counts
    :param `bool` cyclic: True if the peptide is cyclic
    :rtype: `int`
    :returns the score
    '''
    #peptide_spectrum = theoretical_spectrum(peptide, cyclic=cyclic)
    peptide_spectrum = cyclic_spectrum(peptide, zero=True)
    score = 0
    for k,v in spectrum_dict(peptide_spectrum).items():
        if k in spec_dict:
            score += spec_dict[k] if spec_dict[k] <= v else v
    return score


def spectral_convolution(seq):
    data = {}
    seq = sorted(seq)
    #print(seq)
    for i in range(len(seq)-1):
        for n,m in zip(seq[i+1:], seq[:-1]):
            #print(i, n, m)
            diff = n-m
            if diff in data:
                data[diff] += 1
            else:
                data[diff] = 1
    return sorted([(v,k) for k,v in data.items() if k > 0], reverse=True)


def expand_convolution(convolution_tuples):
    '''

    '''
    masses = []
    for n, m in convolution_tuples:
        for i in range(n):
            masses.append(m)
    return masses


def convolution_sequencing(m, n, seq):
   
    if 0 not in seq:
        seq.append(0)
    slots = [slot for slot in spectral_convolution(seq) if slot[1] >= 57 and slot[1] < 200]
    top_m = slots[m-1] if len(slots) > m-1 else slots[-1]
    mass = [v for c, v in slots if c >= top_m[0]]
    #print(slots)    
    #print(mass)    

    return leaderboard_cyclopeptide_sequence(seq, n, allow=mass, cyclic=True)
    

  
def subpeptides(n):
    '''
    calculate the number of subpeptides for a cyclic peptide 
    '''  
    return (n-1)*n 


def masses2peptide(mass):
    pep = [rev_mass[m] for m in mass] 
    return ''.join(pep)
    

def trim(leaderboard, spectrum, N, cyclic=True):
    '''
    pick the N highest scoring peptides from leaderboard with respect to spectrum

    :param [`str`, ...] leaderboard: the list of peptides to choose from 
    :param [`int`, ...] spectrum: the spectrum to compare the peptides to.
    :param `int` N: the top N to return
    :returns: the N highest scoring leaders 
    '''
    if isinstance(spectrum, dict):
        spec_dict = spectrum
    else:
        spec_dict = spectrum_dict(spectrum)
    scores = [(cyclopeptide_scoring(peptide, spec_dict, cyclic=cyclic),peptide) for peptide in leaderboard]
    top_N = sorted(scores, reverse=True)[:N]
    return [top[1] for top in top_N]
    
