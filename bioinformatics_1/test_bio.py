#!/usr/bin/env python
import unittest
from random import shuffle
import bio1, bio2, bio3, bio4, bio5, bio6, util
from graph import *

class TestBioFunctions(unittest.TestCase):

    def test_pp(self):
        genome = 'GATATATGCATATACTT'
        pattern = 'ATAT'
        res = bio1.pattern_start(genome, pattern)
        ss = bio1.pp(res)
        self.assertEqual(ss, "1 3 9")

    def test_kmer(self):
        t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        res = sorted(bio1.kmer(t, k))
        self.assertEqual(res, ['CATG', 'GCAT'])

    def test_kmer_min_count(self):
        t = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        res = sorted(bio1.kmer(t, k, 2))
        self.assertEqual(res, ['ATGA', 'CATG', 'GCAT', 'TGCA'])

    def test_complement(self):
        t = 'AAAACCCGGT'
        res = bio1.complement(t)
        self.assertEqual(res, 'ACCGGGTTTT')

    def test_pattern_start(self):
        genome = 'GATATATGCATATACTT'
        pattern = 'ATAT'
        res = bio1.pattern_start(genome, pattern)
        self.assertEqual(res, [1, 3, 9])

    def test_clumpsO(self):
        genome = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
        k = 5
        clump = 75
        t = 4
        res = bio1.clumps(genome, k, clump, t)
        self.assertEqual(sorted(res), ['AATGT', 'CGACA', 'GAAGA'])

    def test_skew(self):
        genome = 'CATGGGCATCGGCCATACGCC'
        expected = [0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2]
        res = bio1.skew(genome)
        self.assertEqual(res, expected)

    def test_min_skew(self):
        genome = 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
        res = bio1.min_skew(genome)
        self.assertEqual(res, [53, 97])

    def test_approximate_match(self):
        genome = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
        pattern = 'ATTCTGGA'
        mismatch = 3
        res = bio1.approximate_pattern(genome, pattern, mismatch)
        self.assertEqual(res, [6, 7, 26, 27])

    def test_pat2num(self):
        pat = 'GT'
        res = bio1.pat2num(pat)
        self.assertEqual(res, 11)

    def test_num2pat(self):
        n = 5437
        k = 7
        res = bio1.num2pat(n, k)
        self.assertEqual(res, 'CCCATTC')
        
        k = 8
        res = bio1.num2pat(n, k)
        self.assertEqual(res, 'ACCCATTC')

    def test_freq_array(self):
        dna = 'ACGCGGCTCTGAAA'
        k = 2
        res = bio1.frequency_array(dna, k)
        self.assertEqual(res, [2, 1, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 0])

    def test_pattern_count(self):
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

    def test_most_freq_mismatch(self):
        text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        k = 4
        d = 1
        res = sorted(bio1.mismatch_kmers(text, k, d))
        self.assertEqual(res, sorted(['GATG', 'ATGC', 'ATGT']))


#---------------------------------------------------------------

    def test_rna_translate(self):
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        amino = bio2.protein_translate(rna)
        self.assertEqual('MAMAPRTEINSTRING', amino)
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGAGCCAUG'
        amino = bio2.protein_translate(rna)
        self.assertEqual('MAMAPRTEINSTRING', amino)
        amino = bio2.protein_translate(rna, stop=False)
        self.assertEqual('MAMAPRTEINSTRINGXAM', amino)


    def test_to_rna_to_dna(self):
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        dna = bio2.to_dna(rna)
        self.assertEqual(dna.find('U'), -1)
        rna2 = bio2.to_rna(dna)
        self.assertEqual(rna, rna2)
        dna2 = bio2.to_dna(rna2)
        self.assertEqual(dna, dna2)

    def test_peptide_to_rna(self):
        amino = 'MA'
        res = bio2.peptide_to_rna(amino)
        self.assertEqual(res, ['ATGGCT', 'ATGGCC', 'ATGGCA', 'ATGGCG'])


    def test_peptide_encoding(self):
        dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
        amino = 'MA'
        res = sorted(bio2.peptide_encode(dna, amino))
        self.assertEqual(res, ['ATGGCC', 'ATGGCC', 'GGCCAT'])

    def test_cyclic_cutter(self):
        peptide = 'LEQN'
        res = [c for c in bio2.cyclic_cutter(peptide)]
        self.assertEqual(res, ['L', 'E', 'Q', 'N']) 
        res = [cs for cs in bio2.cyclic_cutter(peptide, cycle=2)]
        self.assertEqual(res, ['LE', 'EQ', 'QN', 'NL']) 

    def test_peptide_sequence_mass(self):
        peptide = 'LEQN'
        res = bio2.theoretical_spectrum(peptide)
        self.assertEqual(res, [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484])

    def test_branch_bound_1(self):
        masses = [0, 113, 128, 186, 241, 299, 314, 427]
        res = bio2.branch_bound_1(masses)
        expect = sorted(['186-128-113', '186-113-128', '128-186-113', '128-113-186', '113-186-128', '113-128-186'])
        self.assertEqual(sorted(res), expect)

    def test_count_masses(self):
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

#---------------------------------------------------
# bio 3

    def test_motif_enumerate(self):
        k = 3
        d = 1
        dnas = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
        expected = ['ATA', 'ATT', 'GTT', 'TTT']
        res = sorted(bio3.motif_enumerate(dnas, k, d))
        self.assertEqual(res, expected)

    def test_motif_generator(self):
        d = 1
        kmer = 'ATT'
        res = sorted([motif for motif in util.motif_generator(kmer, d)])
        self.assertEqual(res, ['AAT', 'ACT', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CTT', 'GTT', 'TTT'])

    def test_motif_in_dna(self):
        motifs = ['ATT', 'CTT', 'TTT', 'GTA']
        dnas = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAGTTT']
        self.assertTrue(bio3.motif_in_dna(motifs, dnas))

        motifs = ['ATT', 'CTT', 'TTT', 'GAT']
        dnas = ['TTTGGC', 'TGCCGTA']
        self.assertFalse(bio3.motif_in_dna(motifs, dnas))

    def test_median_string(self):
        k = 3
        dnas = ['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG','GCTGAGCACCGG','AGTACGGGACAG']
        res = bio3.median_string(dnas, k)
        self.assertEqual(res , 'GAC')
    
    def test_most_probable(self):
        dna = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
        k = 5
        profile = [[0.2, 0.2, 0.3, 0.2, 0.3],
                  [0.4, 0.3, 0.1, 0.5, 0.1],
                  [0.3, 0.3, 0.5, 0.2, 0.4],
                  [0.1, 0.2, 0.1, 0.1, 0.2]]
        res = bio3.most_probable(dna, k, profile)
        self.assertEqual(res, 'CCGAG')

    def test_motif_matrix(self):
        dnas = ['GGC','AAG','CAA','CAC','CAA']
        res = bio3.motif_matrix(dnas)
        expect = [[0.2, 0.8, 0.4],[0.6, 0.0, 0.4],[0.2, 0.2, 0.2],[0.0, 0.0, 0.0]]
        self.assertEqual(res, expect)


    def test_greedy_motif(self):
        k = 3 
        t = 5
        dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
        res = bio3.greedy_motif(k, t, dnas)
        excpect = ['CAG','CAG','CAA','CAA','CAA']
        self.assertEqual(res, excpect)

    def test_greedy_motif_pseudocount(self):
        k = 3 
        t = 5
        dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
        res = bio3.greedy_motif(k, t, dnas, pseudo=True)
        excpect = ['TTC','ATC','TTC','ATC','TTC']
        self.assertEqual(res, excpect)

    def test_wrong_greedy(self):
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

    def test_consensus(self):
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

    def test_score(self):
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

    def __test_randomized_motifs(self):
        k = 8
        t = 5
        dna = [
            'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
            'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
            'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
            'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
            'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
        ]
        res = bio3.randomized_motif_search(dna, k, t)
        expected = ['TCTCGGGG','CCAAGGTG','TACAGGCG','TTCAGGTG','TCCACGTG']
        self.assertEqual(res, expected)


#---------------------------------------------------
# graphs
    
    def test_balanced(self):
        edge_iter = ((0,[1,2]), (1,[3]), (2,[1,3]), (3,[0]))
        g = GraphFactory(edge_iter, directed=True)
        self.assertFalse(is_balanced(g))

        edge_iter = ((0,[1,2]), (1,[3]), (2,[3]), (3,[0, 4]), (4,[0]))
        g = GraphFactory(edge_iter, directed=True)
        self.assertTrue(is_balanced(g))

    def test_kmer_gen(self):
        dna = 'CAATCCAAC'
        k = 3
        res = [kmer for kmer in util.kmer_gen(dna, k)]
        expected = ['CAA', 'AAT', 'ATC','TCC', 'CCA','CAA', 'AAC']
        self.assertEqual(res, expected) 
       

#----------------------------------------------------
#lesson 4

    def test_string_complement(self):
        text = 'CAATCCAAC'
        length = 5
        res = sorted([comp for comp in bio4.str_composition(text, length)])
        expect = ['AATCC', 'ATCCA', 'CAATC', 'CCAAC', 'TCCAA']
        self.assertEqual(res, expect)

    def test_overlap(self):
        kmers = ['ATGCG','GCATG','CATGC','AGGCA','GGCAT']
        res = bio4.str_overlap(kmers)
        expect = {'AGGCA': ['GGCAT'], 'CATGC': ['ATGCG'], 'GCATG' : ['CATGC'], 'GGCAT' : ['GCATG']}
        self.assertEqual(res, expect)

        res = bio4.format_overlaps(res)
        expect = 'AGGCA -> GGCAT\nCATGC -> ATGCG\nGCATG -> CATGC\nGGCAT -> GCATG'
        self.assertEqual(res, expect)

    def test_deBrujin(self):
        text = 'AAGATTCTCTAC'
        kmer = 4
        _iter = bio4.str_composition(text, kmer)
        res = bio4.deBruijn(_iter)
        expect = {'AGA': ['GAT'], 'GAT': ['ATT'], 'AAG': ['AGA'], 'TTC': ['TCT'], 'CTC': ['TCT'], 'CTA': ['TAC'], 'ATT': ['TTC'], 'TCT': ['CTA', 'CTC']}
        self.assertEqual(res, expect) 

    def test_euler_cycle(self):
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

    def test_euler_path(self):
        with open('data/euler_2.dat', 'r') as f:
            edge_iter = g_read(f.readlines())
            res = '->'.join([str(node) for node in bio4.eulerPath(edge_iter)])
            expect = '6->7->8->9->6->3->0->2->1->3->4'
            self.assertEqual(res, expect)

    def test_str_reconstruct(self):
        kmers = ['CTTA','ACCA','TACC','GGCT','GCTT','TTAC']
        expect = 'GGCTTACCA'
        res = bio4.str_reconstruct(kmers)
        self.assertEqual(res, expect)

    def test_str_pair_reconstruct(self):
        kmers = [('GAGA','TTGA'),('TCGT','GATG'),('CGTG','ATGT'),('TGGT','TGAG'),('GTGA','TGTT'),('GTGG','GTGA'),('TGAG','GTTG'),('GGTC','GAGA'),('GTCG','AGAT')]
        d = 2
        expect = 'GTGGTCGTGAGATGTTGA'
        res = bio4.str_reconstruct(kmers, pairs=True, distance=d)
        self.assertEqual(res, expect)

    def test_circular(self):
        n = 4
        expect = '0001111011001010'
        res = bio4.universal_circular(n)
        self.assertEqual(res, expect)

    def __test_contigs(self):
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

    def __test_LCS(self):
        data = ['AACCTTGG', 'ACACTGTGA']
        expected = 'AACTGG'
        res = bio5.LCS(data)
        self.assertEqual(res, expected)

#----------------------------------------------------
#lesson 6

    def test_greedy_revesal_sort(self):
        data = [-3, +4, +1, +5, -2]
        res = bio6.greedy_reversal(data)
        res = bio6.nest_list_print(res)
        expect = '(-1 -4 +3 +5 -2)\n(+1 -4 +3 +5 -2)\n(+1 +2 -5 -3 +4)\n(+1 +2 +3 +5 +4)\n(+1 +2 +3 -4 -5)\n(+1 +2 +3 +4 -5)\n(+1 +2 +3 +4 +5)'
        self.assertEqual(res, expect)

    def test_breakpoint_count(self):
        data = [+3, +4, +5, -12, -8, -7, -6, +1, +2, +10, +9, -11, +13, +14]
        res = bio6.breakpoint_count(data)
        self.assertEqual(res, 8)

if __name__ == '__main__':  
    unittest.main()
