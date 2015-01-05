#!/opt/local/bin/pypy


import collections
import functools
import itertools
import sys


sys.setrecursionlimit(10000)


BLOSUM = {'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1,
                'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1,
                'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2},
          'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1,
                'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3,
                'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2},
          'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3,
                'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1,
                'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
          'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3,
                'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1,
                'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3},
          'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4,
                'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2,
                'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3},
          'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0,
                'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4,
                'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3},
          'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4,
                'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3,
                'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1},
          'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3,
                'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2,
                'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2},
          'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3,
                'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1,
                'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},
          'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1,
                'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2,
                'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1},
          'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2,
                'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3,
                'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1},
          'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3,
                'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2,
                'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2},
          'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3,
                'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1,
                'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1},
          'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3,
                'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7,
                'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3},
          'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2,
                'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1,
                'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2},
          'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3,
                'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2,
                'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2},
          'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1,
                'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1,
                'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2},
          'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3,
                'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4,
                'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2},
          'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3,
                'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2,
                'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1},
          'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1,
                'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3,
                'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}}


class Error(Exception):
    """Base error class."""


class DPChange(object):
    INVALID = 10000

    def __init__(self, coins):
        self.coins = set(coins)
        self.num_coins = {coin: 1 for coin in coins}

    def xGetMinNumCoins(self, money):
        # Recursive version.
        if money in self.num_coins:
            return self.num_coins[money]
        self.num_coins[money] = min(
            1+self.GetMinNumCoins(money-coin)
            for coin in self.coins if money > coin)
        return self.num_coins[money]

    def GetMinNumCoins(self, money):
        # Bottom-up, iterative approach.
        for i in xrange(2, money+1):
            if i in self.num_coins:
                continue
            self.num_coins[i] = min(
                1+self.num_coins.get(i-coin, self.INVALID)
                for coin in self.coins if i > coin)
        return self.num_coins[money]


class LinearPeptides(object):

    def __init__(self, mass_list):
        self.mass_list = mass_list
        self.mass = {1: {m: 1 for m in mass_list}}
        self.max_amount = 1

    def GetTotal(self, mass):
        return sum(self.GetMass(mass, i) for i in xrange(1, mass+1))
            
    def GetMass(self, mass, amount):
        if self.mass.get(amount):
            return self.mass[amount].get(mass, 0)
        for i in xrange(self.max_amount+1, amount+1):
            self.UpdateMassAmount(i)
        self.max_amount = amount
        return self.mass[amount].get(mass, 0)

    def UpdateMassAmount(self, i):
        mass_counter = collections.Counter()
        for mass, amount in self.mass[i-1].iteritems():
            for m in self.mass_list:
                mass_counter[mass+m] += amount
        self.mass[i] = mass_counter


class LpDag(object):

    def __init__(self, adj):
        self.GetGraph(adj)

    def GetBackTrack(self, source, orders):
        # Limit the node of interest to those included in the topological
        # sort orders.
        valid = set(orders)
        scores, backtrack = {source: 0}, {}
        for i, start in enumerate(orders[:-1]):
            node = orders[i+1]
            scores[node], backtrack[node] = max(
                (scores[p]+self.graph[p][node], p)
                for p in self.prevs[node]
                if p in valid)
        return backtrack

    def GetGraph(self, adj):
        self.graph = collections.defaultdict(dict)
        self.prevs = collections.defaultdict(list)
        counter = collections.Counter()
        for ad in adj:
            start, end_weight = ad.split('->')
            end, weight = end_weight.split(':')
            self.graph[start][end] = int(weight)
            self.prevs[end].append(start)
            counter[start] += 1
            counter[end] -= 1
        self.sources = [node for node, freq in counter.iteritems()
                        if freq > 0]

    def GetLP(self, source, dest):
        orders = self.SortNodes(source)
        backtrack = self.GetBackTrack(source, orders)
        return self.GetLPFromBackTrack(source, dest, backtrack)

    def GetLPFromBackTrack(self, source, dest, backtrack):
        total = 0
        path = [dest]
        while dest != source:
            path.append(backtrack[dest])
            total += self.graph[path[-1]][dest]
            dest = path[-1]
        return total, '->'.join(reversed(path))

    def SortNodes(self, source):
        orders, stack, visited = [], [source], set()
        # DFS.
        while stack:
            neighbors = [node for node in self.graph[stack[-1]]
                         if node not in visited]
            if neighbors:
                stack.extend(neighbors)
                continue
            node = stack[-1]
            stack = stack[:-1]
            if node in visited: continue
            orders.append(node)
            visited.add(node)
        return list(reversed(orders))


class Alignment(object):
    INDEL = -5
    INVALID = -100
    STOP = (-1, -1)

    def __init__(self, strings):
        self.s1, self.s2 = strings
        self.scores = {(0, 0): 0}
        self.backtrack = {}
        self.max_node = None
        self.max_score = 0

    def GetMaxScore(self, this, c1, c2, left, left_score, top, top_score,
                    diag, diag_score):
        return max((self.GetWeight(left, this, c1, c2) + left_score, left),
                   (self.GetWeight(top, this, c1, c2) + top_score, top),
                   (self.GetWeight(diag, this, c1, c2) + diag_score, diag))
            
    def GetScores(self):
        for i in xrange(len(self.s1) + 1):
            ci = self.s1[i-1] if i >= 1 else None
            for j in xrange(len(self.s2) + 1):
                if (not i) and (not j): continue
                cj = self.s2[j-1] if j >= 1 else None
                left = (i, j-1)
                left_score = self.scores.get(left, self.INVALID)
                top = (i-1, j)
                top_score = self.scores.get(top, self.INVALID)
                diag = (i-1, j-1)
                diag_score = self.scores.get(diag, self.INVALID)
                self.scores[i,j], self.backtrack[i,j] = self.GetMaxScore(
                    (i, j), ci, cj, left, left_score, top, top_score, diag, diag_score)
                if self.scores[i,j] > self.max_score:
                    self.max_score = self.scores[i,j]
                    self.max_node = (i, j)

    def GetWeight(self, unused_node1, unused_node2, unused_c1, unused_c2):
        raise NotImplementedError

    def ShowScores(self):
        for i in xrange(len(self.s1)+1):
            for j in xrange(len(self.s2)+1):
                print '{:>3}'.format(self.scores.get((i, j), '-')),
            print


class GlobalAlignment(Alignment):

    def GetAlignedStrings(self):
        char1_list, char2_list = [], []
        this = (len(self.s1), len(self.s2))
        while this != (0, 0):
            prev = self.backtrack[this]
            prev0, prev1 = prev
            if ((prev0 + 1 == this[0]) and (prev1 + 1 == this[1])):
                char1_list.append(self.s1[prev0])
                char2_list.append(self.s2[prev1])
            elif (prev0 == this[0]):
                char1_list.append('-')
                char2_list.append(self.s2[prev1])
            else:
                char1_list.append(self.s1[prev0])
                char2_list.append('-')
            this = prev
        return [''.join(reversed(char1_list)), ''.join(reversed(char2_list))]
        
    def GetAlignment(self):
        if len(self.scores) <= 1:
            self.GetScores()
        scores = self.scores[len(self.s1), len(self.s2)]
        aligned_strings = self.GetAlignedStrings()
        return scores, aligned_strings
                    
    def GetWeight(self, node1, node2, c1, c2):
        if (c1 is None) or (c2 is None):
            return self.INDEL
        if (node1[0] == node2[0]) or (node1[1] == node2[1]):
            return self.INDEL
        return BLOSUM[c1][c2]


class EditDistance(GlobalAlignment):

    def FindDistance(self):
        distance = 0
        this = (len(self.s1), len(self.s2))
        while this != (0, 0):
            prev = self.backtrack[this]
            prev0, prev1 = prev
            if ((prev0 + 1 != this[0]) or (prev1 + 1 != this[1]) or
                self.s1[prev0] != self.s2[prev1]):
                distance += 1
            this = prev
        return distance

    def GetDistance(self):
        if len(self.scores) <= 1:
            self.GetScores()
        return self.FindDistance()
                    
    def GetWeight(self, node1, node2, c1, c2):
        if (c1 is None) or (c2 is None):
            return 0
        if (node1[0] == node2[0]) or (node1[1] == node2[1]):
            return 0
        # Encourage substitution.
        return 2 if c1 == c2 else 1


class LocalAlignment(Alignment):
    PAM = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1,
                 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1,
                 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3},
           'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2,
                 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3,
                 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0},
           'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2,
                 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1,
                 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4},
           'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2,
                 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1,
                 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4},
           'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3,
                 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0,
                 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5},
           'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1,
                 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5,
                 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7},
           'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5,
                 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2,
                 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1},
           'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2,
                 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0,
                 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0},
           'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2,
                 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1,
                 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4},
           'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2,
                 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2,
                 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2},
           'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2,
                 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3,
                 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1},
           'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2,
                 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0,
                 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2},
           'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2,
                 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0,
                 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4},
           'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2,
                 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6,
                 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5},
           'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1,
                 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1,
                 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3},
           'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2,
                 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0,
                 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4},
           'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0,
                 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0,
                 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3},
           'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5,
                 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6,
                 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0},
           'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4,
                 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1,
                 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2},
           'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1,
                 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5,
                 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}

    def GetAlignedStrings(self):
        char1_list, char2_list = [], []
        this = self.max_node
        score = 0
        while this != (0, 0) and score < self.max_score:
            prev = self.backtrack[this]
            if prev == self.STOP: break
            prev0, prev1 = prev
            if ((prev0 + 1 == this[0]) and (prev1 + 1 == this[1])):
                char1_list.append(self.s1[prev0])
                char2_list.append(self.s2[prev1])
                score += self.PAM[self.s1[prev0]][self.s2[prev1]]
            elif (prev0 == this[0]):
                char1_list.append('-')
                char2_list.append(self.s2[prev1])
                score += self.INDEL
            else:
                char1_list.append(self.s1[prev0])
                char2_list.append('-')
                score += self.INDEL
            this = prev
        return [''.join(reversed(char1_list)), ''.join(reversed(char2_list))]
        
    def GetAlignment(self):
        if len(self.scores) <= 1:
            self.GetScores()
        aligned_strings = self.GetAlignedStrings()
        return self.max_score, aligned_strings

    def GetMaxScore(self, this, c1, c2, left, left_score, top, top_score,
                    diag, diag_score):
        return max((self.GetWeight(left, this, c1, c2) + left_score, left),
                   (self.GetWeight(top, this, c1, c2) + top_score, top),
                   (self.GetWeight(diag, this, c1, c2) + diag_score, diag),
                   (0, self.STOP))
                    
    def GetWeight(self, node1, node2, c1, c2):
        if (c1 is None) or (c2 is None):
            return self.INDEL
        if (node1[0] == node2[0]) or (node1[1] == node2[1]):
            return self.INDEL
        return self.PAM[c1][c2]


class FitAlignment(GlobalAlignment):
    # Introduction to Computational Biology: Maps, Sequences and Genomes
    # By Michael S. Waterman
    # Sec. 9.5 Fitting One Sequence Into Another
    INDEL = -1

    def GetAlignedStrings(self):
        char1_list, char2_list = [], []
        # Along the last column, find the row where the score is the max.
        # This should give us the last character of s1.
        last_s2 = len(self.s2)
        _, row = max((self.scores[i, last_s2], i) for i in xrange(len(self.s1)+1))
        this = (row, last_s2)
        s2_counter = 0
        while this != (0, 0) and s2_counter < last_s2:
            prev = self.backtrack[this]
            prev0, prev1 = prev
            if ((prev0 + 1 == this[0]) and (prev1 + 1 == this[1])):
                char1_list.append(self.s1[prev0])
                char2_list.append(self.s2[prev1])
                s2_counter += 1
            elif (prev0 == this[0]):
                char1_list.append('-')
                char2_list.append(self.s2[prev1])
                s2_counter += 1
            else:
                char1_list.append(self.s1[prev0])
                char2_list.append('-')
            this = prev
        return [''.join(reversed(char1_list)), ''.join(reversed(char2_list))]
        
    def GetAlignment(self):
        if len(self.scores) <= 1:
            self.GetScores()
        aligned_strings = self.GetAlignedStrings()
        return self.GetScore(aligned_strings), aligned_strings

    def GetMaxScore(self, this, c1, c2, left, left_score, top, top_score,
                    diag, diag_score):
        if not this[0]:
            return (this[1] * self.INDEL, top)
        if not this[1]:
            return (0, left)
        return max((self.GetWeight(left, this, c1, c2) + left_score, left),
                   (self.GetWeight(top, this, c1, c2) + top_score, top),
                   (self.GetWeight(diag, this, c1, c2) + diag_score, diag))
        
    def GetScore(self, aligned_strings):
        _score = lambda c1, c2: 1 if c1 == c2 else self.INDEL
        return sum(_score(c1, c2) for c1, c2 in zip(*aligned_strings))
                    
    def GetWeight(self, node1, node2, c1, c2):
        if (c1 is None) or (c2 is None):
            return self.INDEL
        if (node1[0] == node2[0]) or (node1[1] == node2[1]):
            return self.INDEL
        return 1 if c1 == c2 else self.INDEL


class OverlapAlignment(FitAlignment):
    INDEL = -2

    def GetAlignedStrings(self):
        char1_list, char2_list = [], []
        # Along the last row, find the column where the score is the max.
        # This should give us the last character of s1.
        last_s1 = len(self.s1)
        _, col = max((self.scores[last_s1, i], i) for i in xrange(len(self.s2)+1))
        this = (last_s1, col)
        while this[1]:
            prev = self.backtrack[this]
            prev0, prev1 = prev
            if ((prev0 + 1 == this[0]) and (prev1 + 1 == this[1])):
                char1_list.append(self.s1[prev0])
                char2_list.append(self.s2[prev1])
            elif (prev0 == this[0]):
                char1_list.append('-')
                char2_list.append(self.s2[prev1])
            else:
                char1_list.append(self.s1[prev0])
                char2_list.append('-')
            this = prev
        return [''.join(reversed(char1_list)), ''.join(reversed(char2_list))]

    def GetMaxScore(self, this, c1, c2, left, left_score, top, top_score,
                    diag, diag_score):
        if not this[0]:
            return (0, top)
        if not this[1]:
            return (0, left)
        return max((self.GetWeight(left, this, c1, c2) + left_score, left),
                   (self.GetWeight(top, this, c1, c2) + top_score, top),
                   (self.GetWeight(diag, this, c1, c2) + diag_score, diag))


class OutputLcs(Alignment):

    def GetBackTrack(self):
        char_list = []
        this = (len(self.s1), len(self.s2))
        while this != (0, 0):
            prev = self.backtrack[this]
            if ((prev[0] + 1 == this[0]) and (prev[1] + 1 == this[1])):
                char_list.append(self.s1[prev[0]])
            this = prev
        return ''.join(reversed(char_list))
        
    def GetLcs(self):
        if len(self.scores) <= 1:
            self.GetScores()
        return self.GetBackTrack()
                    
    def GetWeight(self, node1, node2, c1, c2):
        if ((node1[0] + 1 == node2[0]) and
            (node1[1] + 1 == node2[1]) and
            c1 == c2):
            return 1
        return 0


class GapAlignment(object):
    EXT_PENALTY = -1
    OPEN_PENALTY = -11
    FROM_LOWER = 0
    FROM_MIDDLE = 1
    FROM_UPPER = 2

    def __init__(self, strings):
        self.s1, self.s2 = strings
        self.lower, self.middle, self.upper = {}, {}, {}
        self.bt_lower, self.bt_middle, self.bt_upper = {}, {}, {}

    def GetAlignedStrings(self):
        chr1_list, chr2_list = [], []
        this = (len(self.s1)-1, len(self.s2)-1)
        bt = self.bt_middle
        bt_map = {self.FROM_LOWER: self.bt_lower,
                  self.FROM_MIDDLE: self.bt_middle,
                  self.FROM_UPPER: self.bt_upper}
        while True:
            if this == (0, 0):
                chr1_list.append(self.s1[this[0]])
                chr2_list.append(self.s2[this[1]])
                break
            bt_key, prev0, prev1 = bt[this]
            if prev0 + 1 == this[0] and prev1 + 1 == this[1]:
                chr1_list.append(self.s1[this[0]])
                chr2_list.append(self.s2[this[1]])
            elif prev0 + 1 == this[0] and prev1 == this[1]:
                chr1_list.append(self.s1[this[0]])
                chr2_list.append('-')
            elif prev0 == this[0] and prev1 + 1 == this[1]:
                chr1_list.append('-')
                chr2_list.append(self.s2[this[1]])
            bt = bt_map[bt_key]
            this = (prev0, prev1)
        return [''.join(reversed(chr1_list)), ''.join(reversed(chr2_list))]

    def GetAlignment(self):
        if not any((self.lower, self.middle, self.upper)):
            self.GetScores()
        aligned_strings = self.GetAlignedStrings()
        return self.middle[len(self.s1)-1,len(self.s2)-1], aligned_strings

    def GetLowerScores(self, i, j):
        prev = (i-1, j)
        self.lower[i,j], self.bt_lower[i,j] = max(
            (self.lower.get(prev, 0) + self.EXT_PENALTY,
             (self.FROM_LOWER, i-1, j)),
            (self.middle.get(prev, 0) + self.OPEN_PENALTY,
             (self.FROM_MIDDLE, i-1, j)))

    def GetMiddleScores(self, i, j):
        ci, cj = self.s1[i], self.s2[j]
        self.middle[i,j], self.bt_middle[i,j] = max(
            (self.lower[i,j], (self.FROM_LOWER, i, j)),
            (self.middle.get((i-1, j-1), 0) + BLOSUM[ci][cj],
             (self.FROM_MIDDLE, i-1, j-1)),
            (self.upper[i,j], (self.FROM_UPPER, i, j)))

    def GetScores(self):
        for i in xrange(len(self.s1)):
            for j in xrange(len(self.s2)):
                self.GetLowerScores(i, j)
                self.GetUpperScores(i, j)
                self.GetMiddleScores(i, j)

    def GetUpperScores(self, i, j):
        prev = (i, j-1)
        self.upper[i,j], self.bt_upper[i,j] = max(
            (self.upper.get(prev, 0) + self.EXT_PENALTY,
             (self.FROM_UPPER, i, j-1)),
            (self.middle.get(prev, 0) + self.OPEN_PENALTY,
             (self.FROM_MIDDLE, i, j-1)))


class LinearSpace(object):
    INDEL = -5
    INVALID = -100

    def __init__(self, strings):
        self.s1, self.s2 = strings
        self.middle = (len(self.s2) - 1) / 2

    def GetScores(self, s1, tail):
        max_row = len(s1) - 1
        this_scores, prev_scores = self.InitScores(s1, tail)
        for col in xrange(-2, -1-len(tail), -1):
            this_scores, prev_scores = prev_scores, this_scores
            cc = tail[col]
            for row in xrange(max_row, -1, -1):
                this_scores[row] = prev_scores[row] + self.INDEL
                # No need to compare other elements if this is the last row.
                if row == max_row: continue
                cr = s1[row]
                this_scores[row] = max(
                    this_scores[row],
                    this_scores[row+1] + self.INDEL,
                    prev_scores[row+1] + BLOSUM[cr][cc])
        return this_scores

    def GetMiddleEdge(self):
        rev_s1 = list(reversed(self.s1))
        rev_s2 = list(reversed(self.s2[:self.middle+1]))
        left_scores = self.GetScores(rev_s1, rev_s2)
        left_scores = list(reversed(left_scores))
        right_scores = self.GetScores(self.s1, self.s2[self.middle+1:])
        rows = self.GetMiddleEdgeRows(left_scores, right_scores)
        return [(rows[0], self.middle), (rows[1], self.middle+1)]

    def GetMiddleEdgeRows(self, left_scores, right_scores):
        max_row0 = len(self.s1) - 1
        max_row1 = max_row0
        max_score = left_scores[max_row1] + right_scores[max_row1] + self.INDEL
        cc = self.s2[self.middle]
        for row in xrange(max_row0-1, -1, -1):
            cr = self.s1[row]
            score, row1 = max(
                (left_scores[row] + right_scores[row] + self.INDEL, row),
                (left_scores[row] + right_scores[row+1] + BLOSUM[cr][cc], row+1))
            if score > max_score:
                max_score = score
                max_row0, max_row1 = row, row1
        return max_row0, max_row1

    def InitScores(self, s1, s2):
        prev_scores = [0] * len(s1)
        this_scores = [0] * len(s1)
        this_scores[-1] = BLOSUM[s1[-1]][s2[-1]]
        for row in xrange(len(s1)-2, -1, -1):
            this_scores[row] += this_scores[row+1] + self.INDEL
        return this_scores, prev_scores


class Manhattan(object):
    INVALID = -100

    def __init__(self, n, m, down_map, right_map):
        self.n = n
        self.m = m
        self.down_map = down_map
        self.right_map = right_map
        self.distances = {(0,0): 0}

    def GetLPLength(self):
        for i in xrange(self.n+1):
            for j in xrange(self.m+1):
                if (i, j) in self.distances: continue
                down_source = (i-1, j)
                down_source_distance = self.distances.get(
                    down_source, self.INVALID) 
                down_weight = self.down_map.get(
                    down_source, self.INVALID)
                right_source = (i, j-1)
                right_source_distance = self.distances.get(
                    right_source, self.INVALID)
                right_weight = self.right_map.get(
                    right_source, self.INVALID)
                self.distances[i,j] = max(
                    down_source_distance + down_weight,
                    right_source_distance + right_weight)
                if self.distances[i,j] < 0:
                    raise Error('Invalid distance: {} for ({}, {}).'.format(
                            self.distances[i,j], i, j))
        return self.distances[self.n, self.m]


class MultiAlignment(object):
    INDEL1 = -1
    INDEL2 = -10
    INVALID = -100
    KEY = '{},{},{}'

    def __init__(self, strings):
        self.s1, self.s2, self.s3 = strings
        self.scores = {'0,0,0': 0}
        self.backtrack = {}

    def GetAlignedStrings(self):
        ch1_list, ch2_list, ch3_list = [], [], []
        this = self.KEY.format(len(self.s1), len(self.s2), len(self.s3))
        while this != '0,0,0':
            prev = self.backtrack[this]
            for th, pr, cl, st in itertools.izip(
                this.split(','),
                prev.split(','),
                [ch1_list, ch2_list, ch3_list],
                [self.s1, self.s2, self.s3]):
                ch = '-' if th == pr else st[int(pr)]
                cl.append(ch)
            this = prev
        return [''.join(reversed(cl))
                for cl in (ch1_list, ch2_list, ch3_list)]

    def GetAlignment(self):
        if len(self.scores) == 1:
            self.GetScores()
        key = self.KEY.format(len(self.s1), len(self.s2), len(self.s3))
        return self.scores[key], self.GetAlignedStrings()

    def GetScores(self):
        for i in xrange(len(self.s1)+1):
            for j in xrange(len(self.s2)+1):
                for k in xrange(len(self.s3)+1):
                    key = self.KEY.format(i, j, k)
                    if key in self.scores: continue
                    # Prefer no gap as highest priority.
                    score = self.Score(i, j, k)
                    prev6 = self.KEY.format(i-1, j-1, k-1)
                    self.scores[key] = self.scores.get(prev6, self.INVALID) + score
                    self.backtrack[key] = prev6
                    # Next, consider 1 gap.
                    prev3 = self.KEY.format(i-1, j-1, k)
                    prev4 = self.KEY.format(i-1, j, k-1)
                    prev5 = self.KEY.format(i, j-1, k-1)
                    score, prev = max(
                        (self.scores.get(prev3, self.INVALID), prev3),
                        (self.scores.get(prev4, self.INVALID), prev4),
                        (self.scores.get(prev5, self.INVALID), prev5))
                    if score > self.scores[key]:
                        self.scores[key] = score
                        self.backtrack[key] = prev
                    # Finally, consider 2 gaps.
                    prev0 = self.KEY.format(i-1, j, k)
                    prev1 = self.KEY.format(i, j-1, k)
                    prev2 = self.KEY.format(i, j, k-1)
                    score, prev = max(
                        (self.scores.get(prev0, self.INVALID), prev0),
                        (self.scores.get(prev1, self.INVALID), prev1),
                        (self.scores.get(prev2, self.INVALID), prev2))
                    if score > self.scores[key]:
                        self.scores[key] = score
                        self.backtrack[key] = prev
                    
    def Score(self, i, j, k):
        if all([i, j, k]):
            cset = set([self.s1[i-1], self.s2[j-1], self.s3[k-1]])
            return int(len(cset) == 1)
        return 0
