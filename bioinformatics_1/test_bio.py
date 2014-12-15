#!/usr/bin/env python
import unittest
from random import shuffle
import bio1, bio2, bio3, bio4, bio5, bio6, util
from graph import *

class TestBioFunctions(unittest.TestCase):

    def _test_pp(self):
        genome = 'GATATATGCATATACTT'
        pattern = 'ATAT'
        res = bio1.pattern_start(genome, pattern)
        ss = bio1.pp(res)
        self.assertEqual(ss, "1 3 9")

    def _test_kmer(self):
        t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        res = sorted(bio1.kmer(t, k))
        self.assertEqual(res, ['CATG', 'GCAT'])

    def _test_kmer_min_count(self):
        t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        res = sorted(bio1.kmer(t, k, 2))
        self.assertEqual(res, ['ATGA', 'CATG', 'GCAT', 'TGCA'])

    def _test_complement(self):
        t = 'AAAACCCGGT'
        res = bio1.complement(t)
        self.assertEqual(res, 'ACCGGGTTTT')

    def _test_pattern_start(self):
        genome = 'GATATATGCATATACTT'
        pattern = 'ATAT'
        res = bio1.pattern_start(genome, pattern)
        self.assertEqual(res, [1, 3, 9])

    def _test_clumpsO(self):
        genome = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
        k = 5
        clump = 75
        t = 4
        res = bio1.clumps(genome, k, clump, t)
        self.assertEqual(sorted(res), ['AATGT', 'CGACA', 'GAAGA'])

    def _test_skew(self):
        genome = 'CATGGGCATCGGCCATACGCC'
        expected = [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]
        res = bio1.skew(genome)
        self.assertEqual(res, expected)

    def _test_min_skew(self):
        genome = 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
        res = bio1.min_skew(genome)
        self.assertEqual(res, [53, 97])

    def _test_approximate_match(self):
        genome = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
        pattern = 'ATTCTGGA'
        mismatch = 3
        res = bio1.approximate_pattern(genome, pattern, mismatch)
        self.assertEqual(res, [6, 7, 26, 27])

    def _test_pat2num(self):
        pat = 'GT'
        res = bio1.pat2num(pat)
        self.assertEqual(res, 11)

    def _test_num2pat(self):
        n = 5437
        k = 7
        res = bio1.num2pat(n, k)
        self.assertEqual(res, 'CCCATTC')
        
        k = 8
        res = bio1.num2pat(n, k)
        self.assertEqual(res, 'ACCCATTC')

    def _test_freq_array(self):
        dna = 'ACGCGGCTCTGAAA'
        k = 2
        res = bio1.frequency_array(dna, k)
        self.assertEqual(res, [2, 1, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 0])

    def _test_pattern_count(self):
        text = 'GCGCG'
        pat  = 'GCG' 
        res = bio1.pattern_count(text, pat)
        self.assertEqual(res, 2)    

        text = 'AACAAGCTGATAAACATTTAAAGAG'
        pat = 'AAAAA'
        d = 2
        res = bio1.pattern_count(text, pat, d)
        self.assertEqual(res, 11)

        text = 'TTTAGAGCCTTCAGAGG'
        pat = 'GAGG'
        d = 2
        res = bio1.pattern_count(text, pat, d)
        self.assertEqual(res, 4)

    def _test_most_freq_mismatch(self):
        text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        d = 1
        res = sorted(bio1.mismatch_kmers(text, k, d))
        self.assertEqual(res, sorted(['GATG', 'ATGC', 'ATGT']))


#---------------------------------------------------------------

    def _test_rna_translate(self):
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        amino = bio2.protein_translate(rna)
        self.assertEqual('MAMAPRTEINSTRING', amino)
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGAGCCAUG'
        amino = bio2.protein_translate(rna)
        self.assertEqual('MAMAPRTEINSTRING', amino)
        amino = bio2.protein_translate(rna, stop=False)
        self.assertEqual('MAMAPRTEINSTRINGXAM', amino)


    def _test_to_rna_to_dna(self):
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        dna = bio2.to_dna(rna)
        self.assertEqual(dna.find('U'), -1)
        rna2 = bio2.to_rna(dna)
        self.assertEqual(rna, rna2)
        dna2 = bio2.to_dna(rna2)
        self.assertEqual(dna, dna2)

    def _test_peptide_to_rna(self):
        amino = 'MA'
        res = bio2.peptide_to_rna(amino)
        self.assertEqual(res, ['ATGGCT', 'ATGGCC', 'ATGGCA', 'ATGGCG'])


    def _test_peptide_encoding(self):
        dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
        amino = 'MA'
        res = sorted(bio2.peptide_encode(dna, amino))
        self.assertEqual(res, ['ATGGCC', 'ATGGCC', 'GGCCAT'])

    def _test_cyclic_cutter(self):
        peptide = 'LEQN'
        res = [c for c in bio2.cyclic_cutter(peptide)]
        self.assertEqual(res, ['L', 'E', 'Q', 'N']) 
        res = [cs for cs in bio2.cyclic_cutter(peptide, cycle=2)]
        self.assertEqual(res, ['LE', 'EQ', 'QN', 'NL']) 

    def _test_peptide_sequence_mass(self):
        peptide = 'LEQN'
        res = bio2.theoretical_spectrum(peptide)
        self.assertEqual(res, [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484])

    def _test_branch_bound_1(self):
        masses = [0, 113, 128, 186, 241, 299, 314, 427]
        res = bio2.cyclo_peptide_sequence(masses)
        expect = sorted(['186-128-113', '186-113-128', '128-186-113', '128-113-186', '113-186-128', '113-128-186'])
        self.assertEqual(sorted(res), expect)

    def _test_count_masses(self):
        in_data = [1,1,2]
        value = 3
        res = bio2.number_of_peptides_with_mass(value, in_data=in_data)
        self.assertEqual(res, 3)

        in_data = [1,5,10]
        value =  16
        res = bio2.number_of_peptides_with_mass(value, in_data=in_data)
        self.assertEqual(res, 58)

        in_data = bio2.masses.values()
        value = 1024
        res = bio2.number_of_peptides_with_mass(value, in_data=in_data)
        self.assertEqual(res, 14712706211)

    def _test_cyclic_spectrum(self):
        self.fail()

        res = bio2.cyclic_spectrum([])
        self.assertEqual(res, [])

        res = bio2.cyclic_spectrum([4])
        self.assertEqual(res, [4])

        res = bio2.cyclic_spectrum([3,4])
        self.assertEqual(res, [3,4,7])

        masses = [97, 101, 97]
        res = bio2.cyclic_spectrum(masses)
        expected = [97, 97, 101, 194, 198, 198, 295]
        self.assertEqual(sorted(res), expected)    

        masses = [97, 99, 103, 97, 101]
        res = sorted(bio2.cyclic_spectrum(masses))
        expected = [97, 97, 99, 101, 103, 196, 198, 200, 202, 299, 299, 301, 396, 400, 497]
        self.assertEqual(res, expected)   
         
    def _test_sepctrum_dict(self):
        masses = [97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295, 297, 299, 299, 301, 394, 396, 398, 400, 400, 497]
        res = bio2.spectrum_dict(masses)
        expected = {97:2, 99:1, 101:1, 103:1, 196:1, 198:2, 200:1, 202:1,
                    295:1, 297:1, 299:2, 301:1, 394:1, 396:1, 398:1, 400:2, 497:1}
        self.assertEqual(res, expected)

    def _test_inner_spectrum(self):
        outer = {97:2, 99:1, 101:1, 103:2, 196:1, 198:2, 200:1, 202:1,
                    295:1, 297:1, 299:2, 301:1, 394:1, 396:1, 398:1, 400:2, 497:1}
        inner = {121:3}
        self.assertFalse(bio2.is_inner_spectrum(inner, outer))

        inner = {103:3}
        self.assertFalse(bio2.is_inner_spectrum(inner, outer))

        inner = {103:2, 295:1, 297:1}
        self.assertTrue(bio2.is_inner_spectrum(inner, outer))

    def _test_cyclopetic_sequence(self):
        self.fail()

        theoretical_masses = [0, 113, 128, 186, 241, 299, 314, 427]
        res = sorted(bio2.cyclo_peptide_sequence(theoretical_masses))
        expect = sorted(['186-128-113', '186-113-128', '128-186-113', '128-113-186', '113-186-128', '113-128-186'])
        self.assertEqual(res, expect)

        theoretical_masses = [0, 97, 97, 99, 101, 103, 196, 198, 198, 200, 202, 295, 297, 299, 299, 301, 394, 396, 398, 400, 400, 497]
        res = sorted(bio2.cyclo_peptide_sequence(theoretical_masses))
        expect = ['101-97-103-99-97', '101-97-99-103-97', '103-97-101-97-99', '103-99-97-101-97', '97-101-97-103-99', '97-101-97-99-103', '97-103-99-97-101', '97-99-103-97-101', '99-103-97-101-97', '99-97-101-97-103']
        self.assertEqual(res, expect)

    def _test_spectral_convolution(self):
        seq = [0, 137, 186, 323]
        res = bio2.spectral_convolution(seq)
        expect = [(2,186), (2,137), (1,323), (1,49)]
        self.assertEqual(res, expect)

    def _test_ConvolutionCyclopeptideSequencing(self):
        M, N = 20, 60
        seq = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]
        res = bio2.convolution_sequencing(M, N, seq)
        expect = '99-71-137-57-72-57'
        self.assertEqual(res, expect)

    def _test_subpeptide_count(self):
        peptide = 'NQEL'
        res = bio2.subpeptides(len(peptide))
        self.assertEqual(res, 12)

        res = bio2.subpeptides(31315)
        self.assertEqual(res, 980597910)

    def _test_theoretical_spectrum(self):
        peptide = 'LEQN'
        res = bio2.theoretical_spectrum(peptide)
        expect = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        self.assertEqual(res, expect)

        peptide = 'NQEL'
        res = bio2.theoretical_spectrum(peptide, cyclic=False)
        expect = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
        self.assertEqual(res, expect)

    def _test_cyclopeptide_scoring(self):
        peptide = 'NQEL'
        spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
        spec_dict = bio2.spectrum_dict(spectrum)
        masses = [bio2.masses[c] for c in peptide]
        res = bio2.cyclopeptide_scoring(masses, spec_dict)
        expect =  11
        self.assertEqual(res, expect)

    def _test_leaderboard_cyclopeptide_sequencing(self):
        N = 10
        seq = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
        res = bio2.leaderboard_cyclopeptide_sequence(seq, N)
        expect = '147-71-129-113'
        self.assertEqual(res, expect)
    
    def _test_liner_peptide_scoring(self):
        peptide = 'NQEL'
        spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
        spec_dict = bio2.spectrum_dict(spectrum)
        peptide = [bio2.masses[c] for c in peptide]
        res = bio2.cyclopeptide_scoring(peptide, spec_dict, cyclic=False)
        expect = 11 
        self.assertEqual(res, expect)

    def _test_leaderboard_trim(self):
        leaderboard = ['LAST', 'ALST', 'TLLT', 'TQAS']
        spectrum = [0, 71, 87, 101, 113, 158, 184, 188, 259, 271, 372]
        N = 2
        leaderboard = [[bio2.masses[c]for c in leader] for leader in leaderboard]
        res = bio2.trim(leaderboard, spectrum, N, cyclic=False)
        expect = [[71, 113, 87, 101], [113, 71, 87, 101]]
        self.assertEqual(sorted(res), expect)

#---------------------------------------------------
# bio 3

    def _test_motif_enumerate(self):
        k = 3
        d = 1
        dnas = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
        expected = ['ATA', 'ATT', 'GTT', 'TTT']
        res = sorted(bio3.motif_enumerate(dnas, k, d))
        self.assertEqual(res, expected)

    def _test_motif_generator(self):
        d = 1
        kmer = 'ATT'
        res = sorted([motif for motif in util.motif_generator(kmer, d)])
        self.assertEqual(res, ['AAT', 'ACT', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CTT', 'GTT', 'TTT'])

    def _test_motif_in_dna(self):
        motifs = ['ATT', 'CTT', 'TTT', 'GTA']
        dnas = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAGTTT']
        self.assertTrue(bio3.motif_in_dna(motifs, dnas))

        motifs = ['ATT', 'CTT', 'TTT', 'GAT']
        dnas = ['TTTGGC', 'TGCCGTA']
        self.assertFalse(bio3.motif_in_dna(motifs, dnas))

    def _test_median_string(self):
        k = 3
        dnas = ['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG','GCTGAGCACCGG','AGTACGGGACAG']
        res = bio3.median_string(dnas, k)
        self.assertEqual(res , 'GAC')
    
    def _test_most_probable(self):
        dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
        k = 5
        profile = [[0.2, 0.2, 0.3, 0.2, 0.3],
                  [0.4, 0.3, 0.1, 0.5, 0.1],
                  [0.3, 0.3, 0.5, 0.2, 0.4],
                  [0.1, 0.2, 0.1, 0.1, 0.2]]
        res = bio3.most_probable(dna, k, profile)
        self.assertEqual(res, 'CCGAG')

    def _test_motif_matrix(self):
        dnas = ['GGC','AAG','CAA','CAC','CAA']
        res = bio3.motif_matrix(dnas)
        expect = [[0.2, 0.8, 0.4],[0.6, 0.0, 0.4],[0.2, 0.2, 0.2],[0.0, 0.0, 0.0]]
        self.assertEqual(res, expect)

    
    def _test_greedy_motif(self):
        self.fail()

        k = 3 
        t = 5
        dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
        res = bio3.greedy_motif(k, t, dnas, pseudo=True)
        excpect = ['CAG','CAG','CAA','CAA','CAA']
        self.assertEqual(res, excpect)

    def _test_greedy_motif_pseudocount(self):
        k = 3 
        t = 5
        dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
        res = bio3.greedy_motif(k, t, dnas, pseudo=True)
        excpect = ['TTC','ATC','TTC','ATC','TTC']
        self.assertEqual(res, excpect)

    def _test_wrong_greedy(self):
        with open('data/rosalind_3d.txt', 'r') as f:
            k, t = [int(n) for n in f.readline().strip().split(' ')]
            dnas = [line.strip() for line in f.readlines()]
            res = bio3.greedy_motif(k, t, dnas)
            not_expect = ['CGTCTAATACGA', 'TGTCTAGCCGTA', 'TAACTGGCGGGA', 'TAAATAAAGAAA',
                          'TGGATTCCGACA', 'TGTAAAAGGGTT', 'TGATCAGTGGGA', 'TGTTTTGTGTCA',
                          'TAACTGGAAACA', 'CTACTAACAATA', 'CGAGTATTGCTA', 'TGGATTACGCGA',
                          'TGATTACTAGTT', 'CATATAACGCTA', 'GGTATGATGGAA', 'TGTTTAGTGTTT',
                          'TGACTGTACCTT', 'TCTTTGGTGTTA', 'GGACAGAGGCCA', 'TATTTAGTGTTA',
                          'CCTTTGAAGGGA', 'TGTATTTCGTGT', 'TATTTGGTGTCA', 'TTTTCGACGCTA',
                          'CGTTTTTCGTGT']
            self.assertNotEqual(res, not_expect)

    def _test_consensus(self):
        motifs = [
            'TCGGGGGTTTTT',
            'CCGGTGACTTAC',
            'ACGGGGATTTTC',
            'TTGGGGACTTTT',
            'AAGGGGACTTCC',
            'TTGGGGACTTCC',
            'TCGGGGATTCAT',
            'TCGGGGATTCCT',
            'TAGGGGAACTAC',
            'TCGGGTATAACC'
        ]
        res = bio3.consensus(motifs)
        self.assertEqual(res, 'TCGGGGATTTCC')

    def _test_score(self):
        motifs = [
            'TCGGGGGTTTTT',
            'CCGGTGACTTAC',
            'ACGGGGATTTTC',
            'TTGGGGACTTTT',
            'AAGGGGACTTCC',
            'TTGGGGACTTCC',
            'TCGGGGATTCAT',
            'TCGGGGATTCCT',
            'TAGGGGAACTAC',
            'TCGGGTATAACC'
        ]
        consensus = 'TCGGGGATTTCC'
        res = bio3.motif_ham_score(consensus, motifs)
        self.assertEqual(res, 30)

    def _test_randomized_motifs(self):
        self.fail()

        k = 8
        t = 5
        dna = [
            'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
            'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
            'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
            'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
            'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
        ]
        N = 3000
        res, score = bio3.randomized_motif_search(dna, k, t, N)
        expected = sorted(['TCTCGGGG','CCAAGGTG','TACAGGCG','TTCAGGTG','TCCACGTG'])
        self.assertEqual(sorted(res), expected)


    def _test_profile_random_motif(self):
        dna = 'CCGGCGTTAG'
        k = 4
        profile = [[3/8, 2/8, 2/8, 2/8],
                   [1/8, 2/8, 2/8, 2/8],
                   [2/8, 2/8, 2/8, 1/8],
                   [2/8, 2/8, 2/8, 3/8]]
        res = defaultdict(int)
        for _ in range(1000):
            motif = bio3.profile_random_motif(dna, k, profile)
            res[motif] += 1
        max_val = max(res.values())
        motif = [k for k,v in res.items() if v == max_val][0]
        self.assertEqual(motif, 'GCGT')

    def _test_gibbs_sampler(self):
        k = 4
        t = 5
        N = 2000
        dna = [
            'TTACCTTAAC',
            'GATGTCTGTC',
            'CCGGCGTTAG',
            'CACTAACGAG',
            'CGTCAGAGGT',
        ]
        res, score = bio3.gibbs_sampler(dna, k, t, N)
        expected = ['ACCT', 'ATGT', 'GCGT', 'ACGA', 'AGGT']
        self.assertEqual(res, expected)

#---------------------------------------------------
# graphs
    
    def _test_balanced(self):
        edge_iter = ((0,[1,2]), (1,[3]), (2,[1,3]), (3,[0]))
        g = GraphFactory(edge_iter, directed=True)
        self.assertFalse(is_balanced(g))

        edge_iter = ((0,[1,2]), (1,[3]), (2,[3]), (3,[0, 4]), (4,[0]))
        g = GraphFactory(edge_iter, directed=True)
        self.assertTrue(is_balanced(g))

    def _test_kmer_gen(self):
        dna = 'CAATCCAAC'
        k = 3
        res = [kmer for kmer in util.kmer_gen(dna, k)]
        expected = ['CAA', 'AAT', 'ATC','TCC', 'CCA','CAA', 'AAC']
        self.assertEqual(res, expected) 
       

#----------------------------------------------------
#lesson 4

    def _test_string_complement(self):
        text = 'CAATCCAAC'
        length = 5
        res = sorted([comp for comp in bio4.str_composition(text, length)])
        expect = ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']
        self.assertEqual(res, expect)

    def _test_overlap(self):
        kmers = ['ATGCG','GCATG','CATGC','AGGCA','GGCAT']
        res = bio4.str_overlap(kmers)
        expect = {'AGGCA': ['GGCAT'], 'CATGC': ['ATGCG'], 'GCATG' : ['CATGC'], 'GGCAT' : ['GCATG']}
        self.assertEqual(res, expect)

        res = bio4.format_overlaps(res)
        expect = 'AGGCA -> GGCAT\nCATGC -> ATGCG\nGCATG -> CATGC\nGGCAT -> GCATG'
        self.assertEqual(res, expect)

    def _test_deBrujin(self):
        self.fail()

        text = 'AAGATTCTCTAC'
        kmer = 4
        _iter = bio4.str_composition(text, kmer)
        res = bio4.deBruijn(_iter)
        expect = {'AGA': ['GAT'], 'GAT': ['ATT'], 'AAG': ['AGA'], 'TTC': ['TCT'], 'CTC': ['TCT'], 'CTA': ['TAC'], 'ATT': ['TTC'], 'TCT': ['CTA', 'CTC']}
        self.assertEqual(res, expect) 

    def _test_euler_cycle(self):
        with open('data/euler.dat', 'r') as f:
            lines = f.readlines()
            idxs = [i for i in range(len(lines))]
            shuffle(idxs)
            for n in idxs:
                edge_iter = g_read(lines)
                res = '->'.join([str(node) for node in bio4.eulerCycle(edge_iter, start_node=n)])
                res += res[1:]
                expect = '6->8->7->9->6->5->4->2->1->0->3->2->6'
                self.assertTrue(expect in res)

    def _test_euler_path(self):
        self.fail()

        with open('data/euler_2.dat', 'r') as f:
            edge_iter = g_read(f.readlines())
            res = '->'.join([str(node) for node in bio4.eulerPath(edge_iter)])
            expect = '6->7->8->9->6->3->0->2->1->3->4'
            self.assertEqual(res, expect)

    def _test_str_reconstruct(self):
        kmers = ['CTTA','ACCA','TACC','GGCT','GCTT','TTAC']
        expect = 'GGCTTACCA'
        res = bio4.str_reconstruct(kmers)
        self.assertEqual(res, expect)

    def _test_str_pair_reconstruct(self):
        kmers = [('GAGA','TTGA'),('TCGT','GATG'),('CGTG','ATGT'),('TGGT','TGAG'),('GTGA','TGTT'),('GTGG','GTGA'),('TGAG','GTTG'),('GGTC','GAGA'),('GTCG','AGAT')]
        d = 2
        expect = 'GTGGTCGTGAGATGTTGA'
        res = bio4.str_reconstruct(kmers, pairs=True, distance=d)
        self.assertEqual(res, expect)

    def _test_circular(self):
        n = 4
        expect = '0001111011001010'
        res = bio4.universal_circular(n)
        self.assertEqual(res, expect)

    def _test_contigs(self):
        self.fail()

        kmers = ['ATG','ATG','TGT','TGG','CAT','GGA','GAT','AGA']
        res = sorted(bio4.contigs(kmers))
        expect = sorted(['AGA', 'ATG', 'ATG', 'CAT', 'GAT', 'TGGA', 'TGT'])
        self.assertEqual(res, expect)

#----------------------------------------------------
#lesson 5

    def test_change(self):
        value = 40
        coins = [1,5,10,20,25,50]
        res = bio5.change(value, coins)
        self.assertEqual(res, 2)

        value = 49
        res = bio5.change(value, coins)
        self.assertEqual(res, 6)


    def test_read_manhatten_data(self):
        path = 'data/manhatten_data.txt'
        g,d,r = bio5.read_manhatten_data(path)
        self.assertEqual(g, [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]])
        self.assertNotEqual(id(g[0]), id(g[1]))
        self.assertEqual(d, ([1, 0, 2, 4, 3], [4, 6, 5, 2, 1], [4, 4, 5, 2, 1], [5, 6, 8, 5, 3]))
        self.assertEqual(r, ([3, 2, 4, 0], [3, 2, 4, 2], [0, 7, 3, 3], [3, 3, 0, 2], [1, 3, 2, 2]))

        path = 'data/rosalind_5b_1_dataset.txt'
        g,d,r = bio5.read_manhatten_data(path)
        self.assertEqual(len(d), 11)
        self.assertEqual(len(d[0]), 8)
        self.assertEqual(len(r), 12)
        self.assertEqual(len(r[0]), 7)


    def test_manhatten_path(self):
        path = 'data/manhatten_data.txt'
        grid, down, right = bio5.read_manhatten_data(path)
        res = bio5.manhatten_path(grid, down, right) 
        self.assertEqual(res, 34)

    def test_LCS(self):
        s1 = 'AACCTTGG'
        s2 = 'ACACTGTGA'
        expected = 'AACTTG'
        res = bio5.LCS(s1, s2)
        self.assertEqual(res, expected)

    def test_longestPath(self):
        f_name = 'data/dataset_245_7.txt'
        res = None
        with open(f_name, 'r') as f:
            src = int(f.readline().strip())
            sink = int(f.readline().strip())
            edges_str = [line.strip() for line in f.readlines()]
            edges = [bio5.make_edge(edge) for edge in edges_str]
            res = bio5.longest_path(src, sink, bio5.DAG(sink, edges))
        expected = (57, [0, 3, 29, 31, 39])
        self.assertEqual(res, expected)

    def test_global_align(self):
        s1 = 'PLEASANTLY'
        s2 = 'MEANLY'
        res = bio5.global_align(s1, s2)
        #expected = (8, [0, 1, 13, 25, 37, 38, 39, 51, 52, 64, 76])
        #self.assertEqual(res, expected)
        res = bio5.align_strings(s1, s2, res[1])
        expected = ['PLEASANTLY', '-ME--AN-LY']
        self.assertEqual(res, expected)

    def test_local_align(self):
        s1 = 'AAAAACMEANLYHRTWWWWWWWW'
        s2 = 'HHHHLPENALTYAAAAA'
        path_len, path = bio5.local_align(s1, s2)
        res = bio5.align_strings(s1, s2, path)
        self.assertEqual(res, ['EANL-Y', 'ENALTY'])


    def test_global_align_strings(self):
        s1 = 'PLEASANTLY'
        s2 = 'MEANLY'
        expected = ['PLEASANTLY', '-MEA--N-LY']
        path = [0, 1, 13, 25, 37, 38, 39, 51, 52, 64, 76]
        res = bio5.align_strings(s1, s2, path)
        self.assertEqual(res, expected)

        s1 = 'MEANLY'
        s2 = 'PENALTY'
        expected = ['EANL-Y','ENALTY']
        path = [8, 16, 24, 32, 40, 47, 55]
        res = bio5.align_strings(s1, s2, path)
        self.assertEqual(res, expected)

    def test_levenstein(self):
        s1 = 'PLEASANTLY'
        s2 = 'MEANLY'
        res = bio5.levenstein(s1, s2)
        self.assertEqual(res, 5)

    def test_fit_align(self):
        s1 = 'GTAGGCTTAAGGTTA'
        s2 = 'TAGATA'
        expected = (2, [7, 24, 25, 42, 43, 60, 77, 94, 111])
        res = bio5.fit_align(s1, s2)
        self.assertEqual(res, expected)

    def test_overlap_align(self):
        s1 = 'PAWHEAE'
        s2 = 'HEAGAWGHEE'
        expected = (1, ['HEAE', 'HEAG'])
        path_len, path = bio5.overlap_align(s1, s2)
        res = bio5.align_affine_strings(s1, s2, path)
        self.assertEqual((path_len, res), expected)

    def test_gap_align(self):
        s1 = 'PRTEINS'
        s2 = 'PRTWPSEIN'
        expected = (8, [0, 9, 18, 27, 115, 123, 131, 51, 60, 69, 78, 239, 79])
        res = bio5.affine_gap(s1, s2)
        self.assertEqual(res, expected)

        res = bio5.align_affine_strings(s1, s2, res[1])
        expected = ['PRT---EINS', 'PRTWPSEIN-']
        self.assertEqual(res, expected)
        
#----------------------------------------------------
#lesson 6

    def _test_greedy_revesal_sort(self):
        data = [-3, +4, +1, +5, -2]
        res = bio6.greedy_reversal(data)
        res = bio6.nest_list_print(res)
        expect = '(-1 -4 +3 +5 -2)\n(+1 -4 +3 +5 -2)\n(+1 +2 -5 -3 +4)\n(+1 +2 +3 +5 +4)\n(+1 +2 +3 -4 -5)\n(+1 +2 +3 +4 -5)\n(+1 +2 +3 +4 +5)'
        self.assertEqual(res, expect)

    def _test_breakpoint_count(self):
        data = [+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]
        res = bio6.breakpoint_count(data)
        self.assertEqual(res, 8)

if __name__ == '__main__':  
    unittest.main()
