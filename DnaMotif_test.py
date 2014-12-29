#!/usr/bin/python

from DnaMotif import MotifEnumeration, DistancePatternStrings, MedianString
from DnaMotif import ProfileProbability, ProfileMostProbableKmer, Score
from DnaMotif import GreedyMotifSearch, GetProfileList, RepeatRandomizedMotifSearch
from DnaMotif import GibbsSampler
from SequenceAntiBiotics_test import ReadFile, WriteFile
import unittest


class FindMotifTesting(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        k, d = [int(i) for i in input[0].split()]
        return input[1:], k, d

    def ReadExpect(self, filename):
        expect = ReadFile(filename)
        return set(expect[0].split())

    def testData1(self):
        dna_list = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
        k = 3
        d = 1
        output = MotifEnumeration(dna_list, k, d)
        expect = set(['ATA', 'ATT', 'GTT', 'TTT'])
        self.assertEqual(expect, output)

    def testData2(self):
        dna_list, k, d = self.ReadInput('motif_enumeration_input.txt')
        output = MotifEnumeration(dna_list, k, d)
        expect = self.ReadExpect('motif_enumeration_output.txt')
        self.assertEqual(expect, output)

    def testData3(self):
        dna_list, k, d = self.ReadInput('dataset_156_7.txt')
        output = MotifEnumeration(dna_list, k, d)
        WriteFile('test_output.txt', output)
        expect = self.ReadExpect('dataset_156_7_output.txt')
        self.assertEqual(expect, output)


class DistancePatternDnaTest(unittest.TestCase):

    def testData1(self):
        pattern = 'AAA'
        dna_list = ['TTACCTTAAC', 'GATATCTGTC', 'ACGGCGTTCG', 'CCCTAAAGAG', 'CGTCAGAGGT']
        self.assertEqual(5, DistancePatternStrings(pattern, dna_list))

    def testData1(self):
        pattern = 'CAAAACG'
        dna_list = 'CATCTTGTCTTCAAGGGTAACCGTTTGTTGCTTCGTTTGATGTCAGAGGCGAGAACCAGAACGATGAAGGGTATAAAACCCTGGATCAATCC AGCCATGGTGTAGCGTCAAGCTACCCATGGGAAACGTAGGTAGAAGAGCATATGCAAACACGTCTTCATAGACTATTAGTCCACTCTGTGGT TGTGGATGTACAATACGGTTAACTCTAAAAAGCTGCCAGAATGGCATAAAAGACGGGTAAGGATAGCCTACGGACTAGAGTAGAACTAGGCT GTAACAACCGGTTTCCAGTGAGCCCGCGGGCTCGAATATGGGAGAAGGCGGCCGCGGGCGACAACCAAAGAGATGTGCACAAATGATTCAAA TTAAGCACTGCGGTATCGGCCTACGCTGGTGGCAAACGTAACGCGCCCATTCCGTTGTCCAGTGTTTGAGCATGAACTCCAAATTGGCAAGT CAAGCTCCGGCAACCGAACTAATGTCTTACCGGCGACCCTAATGATGCGGAAAGTTCCTTATGAGACCTATTCTTTGACGCCAGGGCTATCT TTCCATATGCAATTTATGGTACGCACTGATTGGGCTTGTCCGACGCTTGAAACCATGGCCTTGTCGGTTCTAGCCCTCACTAAGAGGAGCCG TGCGAATTAAGACTTGAGATCTAGCTTAAAGGATTGTAATTTAAAGGGGTAAATCGGTCTATAGCAGAGGAGATTAACTGTGAACAGCCTGG AAAACCTTAATCCCGCGGATGGCGCGCGCCGTAACAGACAGGAGCCATTTACTTAAGCTGATGTATACAGAGAGAATAGAGGCCGTGGGGGT CGGCCTTCTGAACTGGGGGGTACCGCAGGCACCCACCTGTACATGGAGCTTGGACGAACTCGCGCCCGGAAGAATTCATAGCAGGACATTAC TCTCGTCTCACTCTACTCTTGTCGCACGTCCATTGGACAGTATCAAAGTATGGAGATTAGCCAGGCACAGTGGTCTAAAACGGGGCACATTA CTGAAACCTAACAATGTCTATACTTCTAGTATAAAATTCAGTCCATAGTACTCGCTACTTAGGTATACAATGATAGCGCTGCAGCAACCGAA GCCGATTGTCACTTACGTTGCAATATCCTAACCGATCACCGCGTGATCACATAGCGGACGCGGTTCCTATGCGTTTTACAACTGGCGTGAAC TGGACAAGATGCATTACGTAGGTCGCTAAGGTTATCAGCGTTTTCGGGTAATTTTGCTGCTATGAAGAGTCTGTGGTGCTCGAAGGCACGAA AGCGAAATTACTCCCATAGTAGCTTTAGTCTTCAGGTACCACCTGTGGATCGATCGGGGGACCGACACTTAGTTTCATTCAGCATATATCAC GGCGTTACTGGCGTACAACTGAAGCACCGAAATGGCCGCTCGGGATCCGCGGGCTTATCAGAACAGCGTGCATCATGGAGCAGTTAGTGTGA ATCGACATAGTTCCTCCAGTCGGGACTCCGAGCAATATCCCGCACTCCAGCGGACCCGCAAGGGCCCCTGTTGCTTCCTTCGTGTTGTCGAG GCTTGACCGGTCTGAACGGCTTTCTAACAACCCGCGAATGCGAGGCACCCTTTAGGGCATGGATTCCACCCGTATCTTCGGCTCCTGTCCAA CTCACTCTTCGGATAGTAGCTGAGTTGCTCACATTGGTCTCCGAGCTTACTCCAAAAGGAAACGAAGTAGACTTTACAAATTAGCGCACTTT ATACTTCGATAAGCATTCCAGAATAGCACCCTGGTCTGGACTGGCATGCTTTATCGCGTAGACCCACTCTTTGACAAGCGTGTACGGGGACA ACCTTATTCTTTGCCCTCTGTCCTTTTCGCCAGTGGGGGAGGACAGGCTGCCGGAACAAGAGTACCAGGAGGCGACATTGATTTATAAAGAG CTCCCCGGGTAGGATGCGTATGATTAAGATACCTCATGTGGGGGGAATCACTCAACCAAACTTTGATACCAGGTTGACGTTCGCGAGCCTAC TCATCAAGTTATCGCTATATATCCTCCGCCTGTGAGGTGTCCTCGTACGCAACCTATTTGTATAACCCGAATGGCAGCTGTGTCCTCCCACA TTGATTCCTCTAACGGTCTAGCACTCCTGGGGGTATCGGTAAACCCAGTAACGAGGCACGAATGCTGGGTAACTTGTTCGCCACCTTAAATA TGCCCGTAAGGTTACGTGACATGACGGACGCCTCCCTTCCAACGCCAGTTTCAGTTATTAGATCAGCGTACGTTGGGGCCGGGATCATTTCA ATAGTGTCAACATACCTTAGTGCTTTGTAGTGCACCACCTGTCCGGGGCCCCTCTGCGATTTGGATTGCTGTTTCAGAGCGCTCGTCGCGGC CCGAAAGCATGATGTGCGCGTGCCAGGATCATTTGGGACCGACCTGGGCGCAGAAGAGAATTACGACAGGCAATGTGCGTTGTCAGGCATTA CGTGGGATTGTGATTAGCTCTGTGGAGTGTCGGGTCACGCTCGCGTTGAAGACTAATATTGTAGCGTGACATTGGTATGGTTTGCAAAGTCT CTTCTTTCGTTCTGGTCATCGCGAAAGGTGGCTCTCGACCACGGACACGGCATGATTCGCTCCACGGAGACTTACGGTGCAGCGTATGAAGC GTGGCAAAAAGTCAACCAACCTTTTCCTGTTTTTGCCGCACGGGATCCCCGCTTGAAAAAGTCTGTACGAAGACCGTTGTCAGTGCCACGTA'.split()
        self.assertEqual(68, DistancePatternStrings(pattern, dna_list))


class MedianStringTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return int(input[0]), input[1:]

    def testData1(self):
        k = 3
        dna_list = ['AAATTGACGCAT', 'GACGACCACGTT', 'CGTCAGCGCCTG', 'GCTGAGCACCGG', 'AGTACGGGACAG']
        output = MedianString(k, dna_list)
        self.assertIn('GAC', output)

    def testData2(self):
        k, dna_list = self.ReadInput('medium_string_input.txt')
        output = MedianString(k, dna_list)
        self.assertIn('CGGCGA', output)

    def testData3(self):
        k, dna_list = self.ReadInput('dataset_158_9.txt')
        output = MedianString(k, dna_list)
        self.assertIn('ACGAAC', output)

    def testData4(self):
        k = 7
        dna_list = ['CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC',
                    'GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC',
                    'GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG']
        output = MedianString(k, dna_list)
        expect = ['AATCCTA', 'AATCCTA', 'GAACCAC', 'GTAGGAA', 'TAGTTTC']
        self.assertEqual(expect, output)


class ProfileProbabilityTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        profile_list = [[float(p) for p in line.split()]
                        for line in input[2:6]]
        return input[0], int(input[1]), profile_list

    def testData1(self):
        kmer = 'TCGGGGATTTCC'
        profile_list = [[0.2, 0.2, 0, 0, 0, 0, 0.9, 0.1, 0.1, 0.1, 0.3, 0],
                        [0.1, 0.6, 0, 0, 0, 0, 0, 0.4, 0.1, 0.2, 0.4, 0.6],
                        [0, 0, 1, 1, 0.9, 0.9, 0.1, 0, 0, 0, 0, 0],
                        [0.7, 0.2, 0, 0, 0.1, 0.1, 0, 0.5, 0.8, 0.7, 0.3, 0.4]]
        self.assertEqual(0.020575296, ProfileProbability(kmer, profile_list))

    def testData2(self):
        kmer = 'TCGTGGATTTCC'
        profile_list = [[0.2, 0.2, 0, 0, 0, 0, 0.9, 0.1, 0.1, 0.1, 0.3, 0],
                        [0.1, 0.6, 0, 0, 0, 0, 0, 0.4, 0.1, 0.2, 0.4, 0.6],
                        [0, 0, 1, 1, 0.9, 0.9, 0.1, 0, 0, 0, 0, 0],
                        [0.7, 0.2, 0, 0, 0.1, 0.1, 0, 0.5, 0.8, 0.7, 0.3, 0.4]]
        self.assertEqual(0, ProfileProbability(kmer, profile_list))

    def testData3(self):
        kmer = 'ACGGGGATTACC'
        profile_list = [[0.2, 0.2, 0, 0, 0, 0, 0.9, 0.1, 0.1, 0.1, 0.3, 0],
                        [0.1, 0.6, 0, 0, 0, 0, 0, 0.4, 0.1, 0.2, 0.4, 0.6],
                        [0, 0, 1, 1, 0.9, 0.9, 0.1, 0, 0, 0, 0, 0],
                        [0.7, 0.2, 0, 0, 0.1, 0.1, 0, 0.5, 0.8, 0.7, 0.3, 0.4]]
        self.assertEqual(0.0008398080000000002, ProfileProbability(kmer, profile_list))

    def testData4(self):
        text = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
        k = 5
        profile_list = [[0.2, 0.2, 0.3, 0.2, 0.3],
                        [0.4, 0.3, 0.1, 0.5, 0.1],
                        [0.3, 0.3, 0.5, 0.2, 0.4],
                        [0.1, 0.2, 0.1, 0.1, 0.2]]
        self.assertEqual('CCGAG', ProfileMostProbableKmer(text, k, profile_list))

    def testData5(self):
        text, k, profile_list = self.ReadInput('profile_most_1.txt')
        self.assertEqual('TGTCGC', ProfileMostProbableKmer(text, k, profile_list))

    def testData6(self):
        text, k, profile_list = self.ReadInput('dataset_159_3.txt')
        self.assertEqual('GTGTGGCTTAACAAG', ProfileMostProbableKmer(text, k, profile_list))

    def testData7(self):
        kmer = 'CAGTGA'
        profile_list = [[0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
                        [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
                        [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
                        [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]]
        self.assertEqual(0.0108, ProfileProbability(kmer, profile_list))

    def testData8(self):
        kmer = 'TCGGTA'
        profile_list = [[0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
                        [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
                        [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
                        [0.3, 0.1, 0.0, 0.4, 0.5, 0.0]]
        self.assertEqual(0.00405, ProfileProbability(kmer, profile_list))


class ScoreTest(unittest.TestCase):

    def testData1(self):
        motifs = ['TCGGGGgTTTtt',
                  'cCGGtGAcTTaC',
                  'aCGGGGATTTtC',
                  'TtGGGGAcTTtt',
                  'aaGGGGAcTTCC',
                  'TtGGGGAcTTCC',
                  'TCGGGGATTcat',
                  'TCGGGGATTcCt',
                  'TaGGGGAacTaC',
                  'TCGGGtATaaCC']
        self.assertEqual(30, Score(motifs))


class GetProfileListTest(unittest.TestCase):

    def testData1(self):
        motifs = ['TCGGGGGTTTTT',
                  'CCGGTGACTTAC',
                  'ACGGGGATTTTC',
                  'TTGGGGACTTTT',
                  'AAGGGGACTTCC',
                  'TTGGGGACTTCC',
                  'TCGGGGATTCAT',
                  'TCGGGGATTCCT',
                  'TAGGGGAACTAC',
                  'TCGGGTATAACC']
        expect = [(0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0),
                  (0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6),
                  (0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0),
                  (0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4)]
        self.assertEqual(expect, GetProfileList(motifs))


class GreedyMotifSearchTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        k, t = [int(s) for s in input[0].split()]
        return k, t, input[1:]

    def testData1(self):
        k = 3
        t = 5
        dna_list = ['GGCGTTCAGGCA',
                    'AAGAATCAGTCA',
                    'CAAGGAGTTCGC',
                    'CACGTCAATCAC',
                    'CAATAATATTCG']
        output = GreedyMotifSearch(k, t, dna_list)
        expect = ['CAG', 'CAG', 'CAA', 'CAA', 'CAA']
        self.assertEqual(expect, output)

    def testData2(self):
        k, t, dna_list = self.ReadInput('greedy_data_input.txt')
        output = GreedyMotifSearch(k, t, dna_list)
        expect = ReadFile('greedy_data_output.txt')
        self.assertEqual(expect, output)

    def testData3(self):
        k, t, dna_list = self.ReadInput('dataset_159_5.txt')
        output = GreedyMotifSearch(k, t, dna_list)
        expect = ReadFile('dataset_159_5_output.txt')
        self.assertEqual(expect, output)

    def testData4(self):
        k = 3
        t = 5
        dna_list = ['GGCGTTCAGGCA',
                    'AAGAATCAGTCA',
                    'CAAGGAGTTCGC',
                    'CACGTCAATCAC',
                    'CAATAATATTCG']
        output = GreedyMotifSearch(k, t, dna_list, laplace=True)
        expect = ['TTC', 'ATC', 'TTC', 'ATC', 'TTC']
        self.assertEqual(expect, output)

    def testData5(self):
        k, t, dna_list = self.ReadInput('greedy_pseudo_input.txt')
        output = GreedyMotifSearch(k, t, dna_list, laplace=True)
        expect = ReadFile('greedy_pseudo_output.txt')
        self.assertEqual(expect, output)

    def testData6(self):
        k, t, dna_list = self.ReadInput('dataset_160_9.txt')
        output = GreedyMotifSearch(k, t, dna_list, laplace=True)
        expect = ReadFile('dataset_160_9_output.txt')
        self.assertEqual(expect, output)


class RandomizedMotifSearchTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        k, t = [int(s) for s in input[0].split()]
        return k, t, input[1:]

    def testData1(self):
        k = 8
        t = 5
        dna_list = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
        output = RepeatRandomizedMotifSearch(k, t, dna_list, 1000)
        expect = ['TCTCGGGG',
                  'CCAAGGTG',
                  'TACAGGCG',
                  'TTCAGGTG',
                  'TCCACGTG']
        self.assertEqual(expect, output)

    def testData2(self):
        k, t, dna_list = self.ReadInput('randomized_input.txt')
        output = RepeatRandomizedMotifSearch(k, t, dna_list, 1000)
        expect = ReadFile('randomized_output.txt')
        self.assertEqual(expect, output)

    def testData3(self):
        k, t, dna_list = self.ReadInput('dataset_161_5.txt')
        output = RepeatRandomizedMotifSearch(k, t, dna_list, 1000)
        expect = ReadFile('dataset_161_5_output.txt')
        self.assertEqual(expect, output)


class GibbsSamplerTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        k, t, N = [int(s) for s in input[0].split()]
        return k, t, N, input[1:]

    def testData1(self):
        k, t, N = 8, 5, 100
        dna_list = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
                    'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
                    'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
                    'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
                    'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
        # output = GibbsSampler(dna_list, k, t, N)
        output = RepeatRandomizedMotifSearch(k, t, dna_list, 1000)
        expect = ['TCTCGGGG',
                  'CCAAGGTG',
                  'TACAGGCG',
                  'TTCAGGTG',
                  'TCCACGTG']
        self.assertEqual(expect, output)

    def testData2(self):
        k, t, N, dna_list = self.ReadInput('dataset_163_4.txt')
        # output = GibbsSampler(dna_list, k, t, N)
        output = RepeatRandomizedMotifSearch(k, t, dna_list, 1000)
        expect = ReadFile('dataset_163_4_output.txt')
        self.assertEqual(expect, output)


if __name__ == '__main__':
    unittest.main()
