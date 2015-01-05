#!/usr/bin/python

from SignPermutations import BreakPoints, Chrom2Cycle, Cycle2Chrom, GreedySort
from SignPermutations import ColorEdges, Graph2Genome, TwoBreakDistance, TwoBreakSort
from SignPermutations import SharedKmers
from SequenceAntiBiotics_test import ReadFile, WriteFile
import re
import unittest


class GreedySortTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return re.search(r'\((.*)\)', input[0]).group(1).split()

    def testData1(self):
        input = self.ReadInput('greedy_sort0_input.txt')
        gsort = GreedySort(input)
        output = gsort.Sort()
        expect = ReadFile('greedy_sort0_output.txt')
        self.assertEqual(expect, output)

    def testData2(self):
        input = self.ReadInput('greedy_sort1_input.txt')
        gsort = GreedySort(input)
        output = gsort.Sort()
        expect = ReadFile('greedy_sort1_output.txt')
        self.assertEqual(expect, output)

    def testData3(self):
        input = self.ReadInput('dataset_286_3.txt')
        gsort = GreedySort(input)
        output = gsort.Sort()
        # WriteFile('test_output.txt', output)
        expect = ReadFile('dataset_286_3_output.txt')
        self.assertEqual(expect, output)

    def testData4(self):
        input = self.ReadInput('greedy_sort2_input.txt')
        gsort = GreedySort(input)
        output = gsort.Sort()
        self.assertEqual(24, len(output))


class BreakPointsTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return re.search(r'\((.*)\)', input[0]).group(1).split()

    def testData1(self):
        input = self.ReadInput('breakpoint0_input.txt')
        self.assertEqual(8, BreakPoints(input))

    def testData2(self):
        input = self.ReadInput('breakpoint1_input.txt')
        self.assertEqual(178, BreakPoints(input))

    def testData3(self):
        input = self.ReadInput('dataset_287_4.txt')
        self.assertEqual(162, BreakPoints(input))

    def testData4(self):
        input = self.ReadInput('breakpoint2_input.txt')
        self.assertEqual(16, BreakPoints(input))


class Chrom2CycleTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return re.search(r'\((.*)\)', input[0]).group(1).split()

    def ReadOutput(self, filename):
        return tuple(int(c) for c in self.ReadInput(filename))

    def testData1(self):
        input = self.ReadInput('chrom_cycle0_input.txt')
        output = Chrom2Cycle(input)
        expect = self.ReadOutput('chrom_cycle0_output.txt')
        self.assertEqual(expect, output)

    def testData2(self):
        input = self.ReadInput('chrom_cycle1_input.txt')
        output = Chrom2Cycle(input)
        expect = self.ReadOutput('chrom_cycle1_output.txt')
        self.assertEqual(expect, output)


class Cycle2ChromTest(unittest.TestCase):

    def ReadInput(self, filename):
        return tuple(int(c) for c in self.ReadOutput(filename))

    def ReadOutput(self, filename):
        input = ReadFile(filename)
        return tuple(re.search(r'\((.*)\)', input[0]).group(1).split())

    def testData1(self):
        input = self.ReadInput('cycle_chrom0_input.txt')
        output = Cycle2Chrom(input)
        expect = self.ReadOutput('cycle_chrom0_output.txt')
        self.assertEqual(expect, output)

    def testData2(self):
        input = self.ReadInput('cycle_chrom1_input.txt')
        output = Cycle2Chrom(input)
        expect = self.ReadOutput('cycle_chrom1_output.txt')
        self.assertEqual(expect, output)


class ColorEdgesTest(unittest.TestCase):

    def ReadInput(self, filename):
        line = ReadFile(filename)[0]
        return [tuple(g.split())
                for g in re.findall(r'\((.*?)\)', line)]

    def ReadOutput(self, filename):
        line = ReadFile(filename)[0]
        return [tuple(int(e) for e in edge)
                for edge in re.findall(r'(\d+), (\d+)', line)]

    def testData1(self):
        genome = self.ReadInput('color_edge0_input.txt')
        output = ColorEdges(genome)
        expect = self.ReadOutput('color_edge0_output.txt')
        self.assertEqual(expect, output)

    def testData2(self):
        genome = self.ReadInput('color_edge1_input.txt')
        output = ColorEdges(genome)
        expect = self.ReadOutput('color_edge1_output.txt')
        self.assertEqual(expect, output)


class Graph2GenomeTest(unittest.TestCase):

    def Check(self, expect, output):
        expect = [set(g) for g in expect]
        output = [set(g) for g in expect]
        for s in expect:
            self.assertIn(s, output)

    def ReadInput(self, filename):
        line = ReadFile(filename)[0]
        return [tuple(int(e) for e in edge)
                for edge in re.findall(r'(\d+), (\d+)', line)]

    def ReadOutput(self, filename):
        line = ReadFile(filename)[0]
        return [tuple(g.split())
                for g in re.findall(r'\((.*?)\)', line)]

    def testData1(self):
        edges = self.ReadInput('graph_genome0_input.txt')
        output = Graph2Genome(edges)
        expect = self.ReadOutput('graph_genome0_output.txt')
        self.Check(expect, output)

    def testData2(self):
        edges = self.ReadInput('graph_genome1_input.txt')
        output = Graph2Genome(edges)
        expect = self.ReadOutput('graph_genome1_output.txt')
        self.Check(expect, output)


class TwoBreakDistanceTest(unittest.TestCase):

    def ReadGenome(self, line):
        return [tuple(gs.split())
                for gs in re.findall(r'\((.*?)\)', line)]

    def ReadInput(self, filename):
        lines = ReadFile(filename)
        return [self.ReadGenome(line) for line in lines]

    def testData1(self):
        genomes = self.ReadInput('two_break_dist0_input.txt')
        tb = TwoBreakDistance(*genomes)
        self.assertEqual(3, tb.GetDistance())

    def testData2(self):
        genomes = self.ReadInput('two_break_dist1_input.txt')
        tb = TwoBreakDistance(*genomes)
        self.assertEqual(9671, tb.GetDistance())

    def testData3(self):
        genomes = self.ReadInput('dataset_288_4.txt')
        tb = TwoBreakDistance(*genomes)
        self.assertEqual(9951, tb.GetDistance())


class TwoBreakSortTest(unittest.TestCase):

    def ReadGenome(self, line):
        return [tuple(gs.split())
                for gs in re.findall(r'\((.*?)\)', line)]

    def ReadInput(self, filename):
        lines = ReadFile(filename)
        return [self.ReadGenome(line) for line in lines]

    def WriteGenome(self, genome):
        return ''.join('({})'.format(' '.join(chrom))
                       for chrom in genome)

    def WriteOutput(self, filename, genomes):
        output = [self.WriteGenome(g) for g in genomes]
        WriteFile(filename, output)

    def testData1(self):
        genomes = self.ReadInput('two_break_sort0_input.txt')
        tb = TwoBreakSort(*genomes)
        output = tb.Sort()
        expect = self.ReadInput('two_break_sort0_output.txt')
        self.assertEqual(expect, output)

    def testData2(self):
        genomes = self.ReadInput('two_break_sort1_input.txt')
        tb = TwoBreakSort(*genomes)
        output = tb.Sort()
        expect = self.ReadInput('two_break_sort1_output.txt')
        self.assertEqual(expect, output)

    def testData3(self):
        genomes = self.ReadInput('dataset_288_5.txt')
        tb = TwoBreakSort(*genomes)
        output = tb.Sort()
        # self.WriteOutput('test_output.txt', output)
        expect = self.ReadInput('dataset_288_5_output.txt')
        self.assertEqual(expect, output)

    def testData4(self):
        genomes = self.ReadInput('dataset_288_5_1.txt')
        tb = TwoBreakSort(*genomes)
        output = tb.Sort()
        self.WriteOutput('test_output.txt', output)
        expect = self.ReadInput('dataset_288_5_1_output.txt')
        self.assertEqual(expect, output)


class SharedKmersTest(unittest.TestCase):
        
    def ReadInput(self, filename):
        input = ReadFile(filename)
        return int(input[0]), input[1], input[2]

    def ReadOutput(self, filename):
        lines = ReadFile(filename)
        return [tuple(int(s) for s in re.findall(r'\((\d+), (\d+)\)', line)[0])
                for line in lines]

    def WriteOutput(self, filename, output):
        output = ['({}, {})'.format(*tu) for tu in output]
        WriteFile(filename, output)

    def testData1(self):
        k, string1, string2 = self.ReadInput('shared_kmers0_input.txt')
        output = SharedKmers(k, string1, string2)
        expect = self.ReadOutput('shared_kmers0_output.txt')
        self.assertEqual(expect, output)

    def testData2(self):
        k, string1, string2 = self.ReadInput('shared_kmers1_input.txt')
        output = SharedKmers(k, string1, string2)
        expect = self.ReadOutput('shared_kmers1_output.txt')
        self.assertEqual(set(expect), set(output))

    def testData3(self):
        k, string1, string2 = self.ReadInput('dataset_289_5.txt')
        output = SharedKmers(k, string1, string2)
        # self.WriteOutput('test_output.txt', output)
        expect = self.ReadOutput('dataset_289_5_output.txt')
        self.assertEqual(set(expect), set(output))

    def testData4(self):
        k, string1, string2 = self.ReadInput('shared_kmers2_input.txt')
        output = SharedKmers(k, string1, string2)
        self.assertEqual(8, len(output))


if __name__ == '__main__':
    unittest.main()
