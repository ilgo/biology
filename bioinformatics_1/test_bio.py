#!/usr/bin/env python

import unittest
from random import shuffle
import bio1, bio2, bio4, bio5
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
        genome = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC'
        pattern = 'ATTCTGGA'
        mismatch = 3
        res = bio1.pattern_start(genome, pattern, mismatch=mismatch)
        self.assertEqual(res, [6, 7, 26, 27, 78])

    def test_kmer_approximate(self):
        genome = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        kmer = 4
        mismatch = 1
        res = bio1.approximate_kmer(genome, kmer, mismatch=mismatch)
        self.assertEqual(sorted(res), ['ATGC', 'ATGT', 'GATG'])

        #genome = 'AAAACAAACGTGGCGGACGTGGCGGACGTGGCGGACGTGGCGGACCATTTGACGCTCTAAAAAACAAAAAACAAAAAACAAACCATTTGAAAAACAACGCTCTAACGCTCTAACTACTGAGACGTGGCGGCTACTGAGACCATTTGAACCATTTGAACCATTTGAACCATTTGACTACTGAGCGCTCTAAAAAACAAAAAACAAACGTGGCGGCTACTGAGAAAACAACTACTGAGAAAACAACGCTCTAACTACTGAGACCATTTGAACGTGGCGGAAAACAACTACTGAGACCATTTGACTACTGAGACCATTTGAAAAACAACGCTCTAACTACTGAGACCATTTGACTACTGAGAAAACAACGCTCTAAACCATTTGAACCATTTGACGCTCTAAACCATTTGACTACTGAGACGTGGCGGACCATTTGAACCATTTGACGCTCTAAACCATTTGACGCTCTAAACGTGGCGGAAAACAAACCATTTGAACGTGGCGGAAAACAAACGTGGCGGCTACTGAGACGTGGCGGCTACTGAGACGTGGCGGCTACTGAGACCATTTGAAAAACAACTACTGAGCGCTCTAACGCTCTAACGCTCTAACGCTCTAAACCATTTGACGCTCTAACGCTCTAAAAAACAACGCTCTAACGCTCTAAACGTGGCGGACCATTTGACTACTGAGACGTGGCGGACCATTTGAACGTGGCGGAAAACAAAAAACAACGCTCTAAACGTGGCGGACGTGGCGGACCATTTGACGCTCTAAACCATTTGAACGTGGCGGACCATTTGA'
        #kmer = 10
        #mismatch = 3
        #res = bio1.approximate_kmer(genome, kmer, mismatch=mismatch)
        #self.assertNotEqual(res, ['AAAAAACAAA'])

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

    def _test_peptide_to_rna(self):
        amino = 'MA'
        res = bio2.peptide_to_rna(amino)
        self.assertEqual(res, [])


    def _test_peptide_encoding(self):
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

#---------------------------------------------------
# graphs
    
    def test_balanced(self):
        edge_iter = ((0,[1,2]), (1,[3]), (2,[1,3]), (3,[0]))
        g = GraphFactory(edge_iter, directed=True)
        self.assertFalse(is_balanced(g))

        edge_iter = ((0,[1,2]), (1,[3]), (2,[3]), (3,[0, 4]), (4,[0]))
        g = GraphFactory(edge_iter, directed=True)
        self.assertTrue(is_balanced(g))

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


#----------------------------------------------------
#lesson 5

    def _test_lcs(self):

        s2 = 'AACCTTGG'
        s1 = 'ACACTGTGA'
        res = bio5.lcs(s1, s2)
        self.assertEqual(res, 'AACTGG')


if __name__ == '__main__':  
    unittest.main()
