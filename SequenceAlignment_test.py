#!/usr/bin/python

from SequenceAlignment import DPChange, Manhattan, OutputLcs, LpDag
from SequenceAlignment import GlobalAlignment, LocalAlignment, EditDistance
from SequenceAlignment import FitAlignment, OverlapAlignment, GapAlignment
from SequenceAlignment import LinearSpace, LinearPeptides, MultiAlignment
from SequenceAntiBiotics_test import ReadFile, WriteFile
import itertools
import sys
import unittest


sys.setrecursionlimit(300000)


class DPChangeTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        coins = [int(c) for c in input[1].split(',')]
        return int(input[0]), coins

    def testData1(self):
        money = 40
        coins = [50, 25, 20, 10, 5, 1]
        change = DPChange(coins)
        expect = 2
        self.assertEqual(expect, change.GetMinNumCoins(money))

    def testData2(self):
        money, coins = self.ReadInput('change_problem_input.txt')
        change = DPChange(coins)
        expect = 338
        self.assertEqual(expect, change.GetMinNumCoins(money))

    def testData3(self):
        money, coins = self.ReadInput('dataset_243_9.txt')
        change = DPChange(coins)
        expect = 2105
        self.assertEqual(expect, change.GetMinNumCoins(money))

    def testData4(self):
        money = 23
        coins = [3, 2]
        change = DPChange(coins)
        expect = 8
        self.assertEqual(expect, change.GetMinNumCoins(money))


class ManhattanTest(unittest.TestCase):

    def ReadDown(self, n, m, down_list):
        return {(row, column): int(weight)
                for row, weight_line in enumerate(down_list)
                for column, weight in itertools.izip(xrange(n+1),
                                                     weight_line.split())}

    def ReadInput(self, filename):
        input = ReadFile(filename)
        n, m = [int(i) for i in input[0].split()]
        separator = input.index('-')
        down_map = self.ReadDown(n, m, input[1:separator])
        right_map = self.ReadRight(n, m, input[separator+1:])
        return n, m, down_map, right_map

    def ReadRight(self, n, m, right_list):
        return {(row, column): int(weight)
                for row, weight_line in enumerate(right_list)
                for column, weight in itertools.izip(xrange(n),
                                                     weight_line.split())}

    def testData1(self):
        n, m, down_map, right_map = self.ReadInput('manhattan_input.txt')
        man = Manhattan(n, m, down_map, right_map)
        expect = 34
        self.assertEqual(expect, man.GetLPLength())

    def testData2(self):
        n, m, down_map, right_map = self.ReadInput('longest_path_1_input.txt')
        man = Manhattan(n, m, down_map, right_map)
        expect = 84
        self.assertEqual(expect, man.GetLPLength())

    def testData3(self):
        n, m, down_map, right_map = self.ReadInput('dataset_261_9.txt')
        man = Manhattan(n, m, down_map, right_map)
        expect = 62
        self.assertEqual(expect, man.GetLPLength())


class OutputLcsTest(unittest.TestCase):

    def testData1(self):
        strings = ReadFile('lcs_input.txt')
        lcs = OutputLcs(strings)
        expect = 'AACTTG'
        self.assertEqual(expect, lcs.GetLcs())

    def testData2(self):
        strings = ReadFile('longest_common_subsequence_input.txt')
        lcs = OutputLcs(strings)
        output = lcs.GetLcs()
        expect = ReadFile('longest_common_subsequence_output.txt')[0]
        self.assertEqual(expect, output)

    def testData3(self):
        strings = ReadFile('dataset_245_5.txt')
        lcs = OutputLcs(strings)
        output = lcs.GetLcs()
        expect = ReadFile('dataset_245_5_output.txt')[0]
        self.assertEqual(expect, output)

    def testData4(self):
        lcs = OutputLcs(['AGACTG', 'GTACGA'])
        output = lcs.GetLcs()
        expect = 'GACG'
        self.assertEqual(expect, output)


class LpDagTest(unittest.TestCase):

    def ReadInput(self, filename):
        input = ReadFile(filename)
        return input[0], input[1], input[2:]

    def testData1(self):
        source, dest, adj = self.ReadInput('lp_dag_input.txt')
        dag = LpDag(adj)
        output = dag.GetLP(source, dest)
        self.assertEqual(9, output[0])
        self.assertEqual('0->2->3->4', output[1])

    def testData2(self):
        source, dest, adj = self.ReadInput('longest_path_in_DAG_input.txt')
        dag = LpDag(adj)
        output = dag.GetLP(source, dest)
        expect = ReadFile('longest_path_in_DAG_output.txt')
        self.assertEqual(int(expect[0]), output[0])
        self.assertEqual(expect[1], output[1])

    def testData3(self):
        source, dest, adj = self.ReadInput('dataset_245_7.txt')
        dag = LpDag(adj)
        output = dag.GetLP(source, dest)
        expect = ReadFile('dataset_245_7_output.txt')
        self.assertEqual(int(expect[0]), output[0])
        self.assertEqual(expect[1], output[1])

    def testData4(self):
        source, dest, adj = self.ReadInput('lp_dag1_input.txt')
        dag = LpDag(adj)
        output = dag.GetLP(source, dest)
        self.assertEqual(14, output[0])
        self.assertEqual('a->d->e->g', output[1])


class GlobalAlignmentTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('global_alignment_input.txt')
        alignment = GlobalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        expect = ReadFile('global_alignment_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertListEqual(expect[1:], aligned_strings)

    def testData2(self):
        input = ReadFile('global_alignment1_input.txt')
        alignment = GlobalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        expect = ReadFile('global_alignment1_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertListEqual(expect[1:], aligned_strings)

    def testData3(self):
        input = ReadFile('dataset_247_3.txt')
        alignment = GlobalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        expect = ReadFile('dataset_247_3_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertListEqual(expect[1:], aligned_strings)

    def testData4(self):
        input = ReadFile('linear_space3_input.txt')
        alignment = GlobalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        # WriteFile('test_output.txt', [str(score)] + aligned_strings)
        expect = ReadFile('linear_space3_output.txt')
        self.assertEqual(int(expect[0]), score)
        # self.assertListEqual(expect[1:], aligned_strings)

    def testData5(self):
        input = ReadFile('dataset_250_14.txt')
        alignment = GlobalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        # WriteFile('test_output.txt', [str(score)] + aligned_strings)
        expect = ReadFile('dataset_250_14_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertListEqual(expect[1:], aligned_strings)

    def testData6(self):
        input = ReadFile('linear_space0_input.txt')
        alignment = GlobalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        print 'score: {}'.format(score)
        print 'aligned_strings: {}'.format(aligned_strings)


class LocalAlignmentTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('local_alignment0_input.txt')
        alignment = LocalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        expect = ReadFile('local_alignment0_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertListEqual(expect[1:], aligned_strings)

    def testData2(self):
        input = ReadFile('local_alignment_input.txt')
        alignment = LocalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        expect = ReadFile('local_alignment_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertListEqual(expect[1:], aligned_strings)

    def testData3(self):
        input = ReadFile('dataset_247_9.txt')
        alignment = LocalAlignment(input)
        score, aligned_strings = alignment.GetAlignment()
        expect = ReadFile('dataset_247_9_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertListEqual(expect[1:], aligned_strings)


class EditDistanceTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('edit_distance0_input.txt')
        distance = EditDistance(input)
        self.assertEqual(5, distance.GetDistance())

    def testData2(self):
        input = ['ATGTTAT', 'ATCGTAC']
        distance = EditDistance(input)
        self.assertEqual(3, distance.GetDistance())

    def testData3(self):
        input = ReadFile('edit_distance1_input.txt')
        distance = EditDistance(input)
        self.assertEqual(400, distance.GetDistance())

    def testData4(self):
        input = ReadFile('dataset_248_3.txt')
        distance = EditDistance(input)
        self.assertEqual(426, distance.GetDistance())


class FitAlignmentTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('fit_alignment0_input.txt')
        fit = FitAlignment(input)
        score, aligned_strings = fit.GetAlignment()
        expect = ReadFile('fit_alignment0_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)

    def testData2(self):
        input = ReadFile('fit_alignment1_input.txt')
        fit = FitAlignment(input)
        score, aligned_strings = fit.GetAlignment()
        expect = ReadFile('fit_alignment1_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)

    def testData3(self):
        input = ReadFile('dataset_248_5.txt')
        fit = FitAlignment(input)
        score, aligned_strings = fit.GetAlignment()
        # WriteFile('test_output.txt', [str(score)] + aligned_strings)
        expect = ReadFile('dataset_248_5_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)


class OverlapAlignmentTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('overlap_alignment0_input.txt')
        overlap = OverlapAlignment(input)
        score, aligned_strings = overlap.GetAlignment()
        expect = ReadFile('overlap_alignment0_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)

    def testData2(self):
        input = ReadFile('overlap_alignment1_input.txt')
        overlap = OverlapAlignment(input)
        score, aligned_strings = overlap.GetAlignment()
        expect = ReadFile('overlap_alignment1_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)

    def testData3(self):
        input = ReadFile('dataset_248_7.txt')
        overlap = OverlapAlignment(input)
        score, aligned_strings = overlap.GetAlignment()
        expect = ReadFile('dataset_248_7_output.txt')
        # WriteFile('test_output.txt', [str(score)] + aligned_strings)
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)


class GapAlignmentTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('gap_alignment0_input.txt')
        gap = GapAlignment(input)
        score, aligned_strings = gap.GetAlignment()
        expect = ReadFile('gap_alignment0_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)

    def testData2(self):
        input = ReadFile('affine_gap_input.txt')
        gap = GapAlignment(input)
        score, aligned_strings = gap.GetAlignment()
        expect = ReadFile('affine_gap_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)

    def testData3(self):
        input = ReadFile('dataset_249_8.txt')
        gap = GapAlignment(input)
        score, aligned_strings = gap.GetAlignment()
        # WriteFile('test_output.txt', [str(score)] + aligned_strings)
        expect = ReadFile('dataset_249_8_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)


class LinearSpaceTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('linear_space0_input.txt')
        linear = LinearSpace(input)
        edge = linear.GetMiddleEdge()
        output = ' '.join(str(node) for node in edge)
        expect = '(4, 3) (5, 4)'
        self.assertEqual(expect, output)

    def testData2(self):
        input = ReadFile('linear_space1_input.txt')
        linear = LinearSpace(input)
        edge = linear.GetMiddleEdge()
        output = ' '.join(str(node) for node in edge)
        expect = '(512, 510) (513, 511)'
        self.assertEqual(expect, output)

    def testData3(self):
        input = ReadFile('dataset_250_12.txt')
        linear = LinearSpace(input)
        edge = linear.GetMiddleEdge()
        output = ' '.join(str(node) for node in edge)
        expect = '(533, 515) (534, 516)'
        self.assertEqual(expect, output)


class LinearPeptidesTest(unittest.TestCase):

    def testData1(self):
        lp = LinearPeptides([2, 3])
        self.assertEqual(1, lp.GetTotal(2))
        self.assertEqual(1, lp.GetTotal(3))
        self.assertEqual(1, lp.GetTotal(4))
        self.assertEqual(2, lp.GetTotal(5))
        self.assertEqual(2, lp.GetTotal(6))
        self.assertEqual(3, lp.GetTotal(7))
        self.assertEqual(4, lp.GetTotal(8))
        self.assertEqual(5, lp.GetTotal(9))
        self.assertEqual(265, lp.GetTotal(23))


class MultiAlignmentTest(unittest.TestCase):

    def testData1(self):
        input = ReadFile('multi_alignment0_input.txt')
        multi = MultiAlignment(input)
        score, aligned_strings = multi.GetAlignment()
        expect = ReadFile('multi_alignment0_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)

    def testData2(self):
        input = ReadFile('dataset_251_5.txt')
        multi = MultiAlignment(input)
        score, aligned_strings = multi.GetAlignment()
        expect = ReadFile('dataset_251_5_output.txt')
        self.assertEqual(int(expect[0]), score)
        self.assertEqual(expect[1:], aligned_strings)


if __name__ == '__main__':
    unittest.main()
