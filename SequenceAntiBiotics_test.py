#!/usr/bin/python

import collections
from SequenceAntiBiotics import Translation, FindSubstringsPeptideEncoding
from SequenceAntiBiotics import LinearSpectrum, CycloPeptideSequencing
from SequenceAntiBiotics import TheoreticalSpectrum, CycloPeptideScoring
from SequenceAntiBiotics import LinearScore, Trim, Convolution
from SequenceAntiBiotics import LeaderboardCycloPeptideSequencing
from SequenceAntiBiotics import ConvolutionCycloPeptideSequencing
import unittest


def ReadFile(filename):
    with open(filename, 'r') as fd:
        return [line.strip() for line in fd]


def WriteFile(filename, result):
    with open(filename, 'w') as fd:
        for res in result:
            fd.write(res+'\n')


class TranslationTest(unittest.TestCase):

    def testData1(self):
        pattern = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        expect = 'MAMAPRTEINSTRING'
        self.assertEqual(expect, Translation(pattern))

    def testData2(self):
        pattern = ReadFile('protein_translate_input.txt')[0]
        expect = ReadFile('protein_translate_output.txt')[0]
        self.assertEqual(expect, Translation(pattern))

    def testData3(self):
        pattern = ReadFile('dataset_96_5.txt')[0]
        expect = ReadFile('ouput_96_5.txt')[0]
        self.assertEqual(expect, Translation(pattern))

    def testData4(self):
        pattern = 'CCACGUACUGAAAUUAAC'
        expect = 'PRTEIN'
        self.assertEqual(expect, Translation(pattern))


class PeptideEncodingTest(unittest.TestCase):

    def testData1(self):
        dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
        peptide = 'MA'
        expect = ['ATGGCC', 'ATGGCC', 'GGCCAT']
        output = sorted(FindSubstringsPeptideEncoding(dna, peptide))
        self.assertEqual(expect, output)

    def testData2(self):
        dna = ReadFile('peptide_encoding_input.txt')[0]
        peptide = 'KEVFEPHYY'
        expect = sorted(ReadFile('peptide_encoding_output.txt'))
        output = sorted(FindSubstringsPeptideEncoding(dna, peptide))
        self.assertEqual(expect, output)

    def testData3(self):
        dna = ReadFile('dataset_96_8_input.txt')[0]
        peptide = 'MKQCTWHL'
        expect = sorted(ReadFile('dataset_96_8_output.txt'))
        output = sorted(FindSubstringsPeptideEncoding(dna, peptide))
        self.assertEqual(expect, output)


class LinearSpectrumTest(unittest.TestCase):

    def GetExpect(self, filename):
        contents = ReadFile(filename)
        return [int(c) for line in contents for c in line.split()]

    def WriteFile(self, filename, output):
        WriteFile(filename, [' '.join(str(i) for i in output)])

    def testData1(self):
        peptide = 'NQEL'
        expect = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
        self.assertEqual(expect, LinearSpectrum(peptide))

    def testData2(self):
        peptide = 'MPYENCCCWMFNIRKGQPDFFRKGAVPYVVPMNCIRWS'
        expect = self.GetExpect('LinearSpectrum_output.txt')
        output = LinearSpectrum(peptide)
        self.WriteFile('test_output.txt', output)
        self.assertEqual(expect, output)

    def testData3(self):
        peptide = 'ATFFLQWTMPSICWYNVTSHWDWKFASYRSMFAGDTPTHNY'
        expect = self.GetExpect('dataset_4912_2_output.txt')
        output = LinearSpectrum(peptide)
        self.assertEqual(expect, output)


class CycloPeptideSequencingTest(unittest.TestCase):

    def testData1(self):
        spectrum = [0, 113, 128, 186, 241, 299, 314, 427]
        output = set(CycloPeptideSequencing(spectrum))
        expect = set(['186-128-113', '186-113-128', '128-186-113',
                      '128-113-186', '113-186-128', '113-128-186'])
        self.assertEqual(output, expect)

    def testData2(self):
        spectrum = ReadFile('cycloseq_input.txt')[0]
        spectrum = [int(m) for m in spectrum.split()]
        expect = set(ReadFile('cycloseq_output.txt')[0].split())
        output = CycloPeptideSequencing(spectrum)
        WriteFile('test_output.txt', output)
        self.assertEqual(set(output), expect)

    def testData3(self):
        spectrum = ReadFile('dataset_100_5.txt')[0]
        spectrum = [int(m) for m in spectrum.split()]
        expect = set(ReadFile('dataset_100_5_output.txt'))
        output = CycloPeptideSequencing(spectrum)
        self.assertEqual(set(output), expect)


class TheoreticalSpectrumTest(unittest.TestCase):

    def testData1(self):
        peptide = 'NQEL'
        expect = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        output = TheoreticalSpectrum(peptide)
        self.assertEqual(output, expect)

    def testData2(self):
        peptide = 'LEQN'
        expect = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
        output = TheoreticalSpectrum(peptide)
        self.assertEqual(output, expect)

    def testData3(self):
        peptide = 'IAQMLFYCKVATN'
        expect = '0 71 71 99 101 103 113 113 114 128 128 131 147 163 170 172 184 199 215 227 227 231 244 259 260 266 271 286 298 298 310 312 328 330 330 372 385 391 394 399 399 399 401 413 423 426 443 443 470 493 498 502 513 519 526 527 541 554 556 557 564 569 590 598 616 626 640 654 657 658 665 670 682 697 697 703 711 729 729 753 753 771 779 785 785 800 812 817 824 825 828 842 856 866 884 892 913 918 925 926 928 941 955 956 963 969 980 984 989 1012 1039 1039 1056 1059 1069 1081 1083 1083 1083 1088 1091 1097 1110 1152 1152 1154 1170 1172 1184 1184 1196 1211 1216 1222 1223 1238 1251 1255 1255 1267 1283 1298 1310 1312 1319 1335 1351 1354 1354 1368 1369 1369 1379 1381 1383 1411 1411 1482'
        expect = [int(m) for m in expect.split()]
        output = TheoreticalSpectrum(peptide)
        self.assertEqual(output, expect)

    def testData4(self):
        peptide = 'QRQLGNTLMHANCR'
        expect = '0 57 71 101 103 113 113 114 114 128 128 131 137 156 156 170 171 185 208 214 215 217 241 244 259 268 272 284 284 284 284 288 298 322 328 339 345 373 381 385 385 387 397 412 412 425 440 444 452 453 454 459 482 498 501 513 516 525 543 553 556 566 568 568 572 581 582 596 626 629 653 657 667 667 669 669 671 681 696 709 712 724 728 738 757 766 770 781 782 784 785 797 825 837 838 840 841 852 856 865 884 894 898 910 913 926 941 951 953 953 955 955 965 969 993 996 1026 1040 1041 1050 1054 1054 1056 1066 1069 1079 1097 1106 1109 1121 1124 1140 1163 1168 1169 1170 1178 1182 1197 1210 1210 1225 1235 1237 1237 1241 1249 1277 1283 1294 1300 1324 1334 1338 1338 1338 1338 1350 1354 1363 1378 1381 1405 1407 1408 1414 1437 1451 1452 1466 1466 1485 1491 1494 1494 1508 1508 1509 1509 1519 1521 1551 1565 1622'
        expect = [int(m) for m in expect.split()]
        output = TheoreticalSpectrum(peptide)
        self.assertEqual(output, expect)

    def testData5(self):
        peptide = 'ALTM'
        expect = '0 71 113 101 131 184 202 214 232 285 303 315 345 416'
        expect = sorted(int(m) for m in expect.split())
        output = TheoreticalSpectrum(peptide)
        self.assertEqual(output, expect)


class CycloPeptideScoringTest(unittest.TestCase):

    def ReadFile(self, filename):
        result = ReadFile(filename)
        peptide = result[0]
        spectrum = [int(m) for m in result[1].split()]
        return peptide, spectrum

    def testData1(self):
        peptide = 'NQEL'
        spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
        output = CycloPeptideScoring(peptide, spectrum)
        expect = 11
        self.assertEqual(expect, output)

    def testData2(self):
        peptide, spectrum = self.ReadFile('CyclicScoring_input.txt')
        output = CycloPeptideScoring(peptide, spectrum)
        expect = 521
        self.assertEqual(expect, output)

    def testData3(self):
        peptide, spectrum = self.ReadFile('dataset_102_3.txt')
        output = CycloPeptideScoring(peptide, spectrum)
        expect = 713
        self.assertEqual(expect, output)


class LinearScoreTest(unittest.TestCase):

    def ReadFile(self, filename):
        result = ReadFile(filename)
        peptide = result[0]
        spectrum = [int(m) for m in result[1].split()]
        return peptide, spectrum

    def testData1(self):
        peptide = 'NQEL'
        spectrum = [0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484]
        output = LinearScore(peptide, spectrum)
        expect = 8
        self.assertEqual(expect, output)

    def testData2(self):
        peptide, spectrum = self.ReadFile('CyclicScoring_input.txt')
        output = LinearScore(peptide, spectrum)
        expect = 344
        self.assertEqual(expect, output)

    def testData3(self):
        peptide, spectrum = self.ReadFile('dataset_4913_1.txt')
        output = LinearScore(peptide, spectrum)
        expect = 232
        self.assertEqual(expect, output)


class TrimTest(unittest.TestCase):

    def ReadFile(self, filename):
        result = ReadFile(filename)
        leaderboard = result[0].split()
        spectrum = [int(m) for m in result[1].split()]
        return leaderboard, spectrum, int(result[2])

    def ReadExpectFile(self, filename):
        result = ReadFile(filename)
        return set(result)

    def testData1(self):
        leaderboard = ['LAST', 'ALST', 'TLLT', 'TQAS']
        spectrum = [0, 71, 87, 101, 113, 158, 184, 188, 259, 271, 372]
        N = 2
        output = Trim(leaderboard, spectrum, N)
        expect = set(['LAST', 'ALST'])
        self.assertEqual(expect, set(output))

    def testData2(self):
        leaderboard, spectrum, N = self.ReadFile('Trim_input.txt')
        output = Trim(leaderboard, spectrum, N)
        expect = set(['WASIGAIMRSAKDMYESLEFHKTHCTYFVYMVCKEARPGWTFFIEWV',
                      'YYGYRQCSWCQRWTVRRMLCWIDVLHKALHWHVCLLFHQALYGFSHE',
                      'WDTDTFFQKAMLKKDETADQIFNLRPYSLTCHNENILGNDNQEKQAG',
                      'YDSPTTYLSTHCHRLTNRMVHENPVICPPQDFAKYLIQSGWEFPLVA',
                      'MAPRDIRMYFDKYHETAALDSQWIIQQIYHLMNVRKLNRTNRFTSVG',
                      'FEKYHQQQILIDAQRVRLVHTVARAGPGWVQTGGWQQTCPRYKPYAW'])
        self.assertEqual(expect, set(output))

    def testData3(self):
        leaderboard, spectrum, N = self.ReadFile('dataset_4913_3.txt')
        output = Trim(leaderboard, spectrum, N)
        expect = self.ReadExpectFile('dataset_4913_3_output.txt')
        self.assertEqual(expect, set(output))


class LeaderboardCycloPeptideSequencingTest(unittest.TestCase):

    def ReadFile(self, filename):
        result = ReadFile(filename)
        spectrum = [int(m) for m in result[1].split()]
        return int(result[0]), spectrum

    def testData1(self):
        N = 10
        spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
        output = LeaderboardCycloPeptideSequencing(spectrum, N)
        expect = '113-147-71-129'
        self.assertEqual(expect, output)

    def testData2(self):
        N, spectrum = self.ReadFile('dataset_102_7.txt')
        output = LeaderboardCycloPeptideSequencing(spectrum, N)
        expect = '115-101-71-186-137-103-128-129-186-115-87-71-87-163-137-131-163-99-87-71-115-57-71'
        self.assertEqual(expect, output)


class ConvolutionTest(unittest.TestCase):

    def testData1(self):
        spectrum = [0, 137, 186, 323]
        output = Convolution(spectrum)
        expect = [137, 137, 186, 186, 323, 49]
        expect_counter = collections.Counter(expect)
        observe_counter = collections.Counter(output)
        self.assertEqual(expect_counter, observe_counter)

    def testData2(self):
        spectrum = ReadFile('spectral_convolution_input.txt')[0]
        spectrum = [int(m) for m in spectrum.split()]
        expect = ReadFile('spectral_convolution_output.txt')[0]
        expect = [int(m) for m in expect.split()]
        output = Convolution(spectrum)
        WriteFile('test_output.txt',
                  [' '.join(str(m) for m in output)])
        expect_counter = collections.Counter(expect)
        observe_counter = collections.Counter(output)
        self.assertEqual(expect_counter, observe_counter)

    def testData3(self):
        spectrum = ReadFile('dataset_104_4.txt')[0]
        spectrum = [int(m) for m in spectrum.split()]
        expect = ReadFile('dataset_104_4_output.txt')[0]
        expect = [int(m) for m in expect.split()]
        output = Convolution(spectrum)
        expect_counter = collections.Counter(expect)
        observe_counter = collections.Counter(output)
        self.assertEqual(expect_counter, observe_counter)


class ConvolutionCycloPeptideSequencingTest(unittest.TestCase):

    def ReadFile(self, filename):
        result = ReadFile(filename)
        spectrum = [int(m) for m in result[2].split()]
        return int(result[0]), int(result[1]), spectrum

    def testData1(self):
        M = 20
        N = 60
        spectrum = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285,
                    299, 307, 323, 356, 364, 394, 422, 493]
        output = ConvolutionCycloPeptideSequencing(M, N, spectrum)
        expect = '99-71-137-57-72-57'
        self.assertIn(expect, output)

    def testData2(self):
        M, N, spectrum = self.ReadFile('convolution_input.txt')
        output = ConvolutionCycloPeptideSequencing(M, N, spectrum)
        expect = '113-115-114-128-97-163-131-129-129-147-57-57-129'
        self.assertIn(expect, output)

    def testData3(self):
        M, N, spectrum = self.ReadFile('dataset_104_7.txt')
        output = ConvolutionCycloPeptideSequencing(M, N, spectrum)
        expect = '115-128-57-137-129-114-97-97-156-99-128-186'
        self.assertIn(expect, output)


if __name__ == '__main__':
    unittest.main()
