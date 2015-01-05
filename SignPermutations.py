from ReverseComplement import ReverseComplement
import collections
import itertools


class Error(Exception):
    """Base error class."""


def BreakPoints(permutation):
    perm_list = [0] + [int(e) for e in permutation] + [len(permutation) + 1]
    return sum(element + 1 != perm_list[index+1]
               for index, element in enumerate(perm_list[:-1]))


def Chrom2Cycle(chrom):
    return tuple(itertools.chain.from_iterable(ChromTuple(c) for c in chrom))


def ChromElement(val1, val2):
    val = max(val1, val2) / 2
    return '+{}'.format(val) if val2 > val1 else '-{}'.format(val)


def ChromTuple(chrom):
    val1 = int(chrom)
    val2 = 2 * abs(val1)
    tu = (val2 -1, val2)
    return tu if val1 > 0 else (tu[1], tu[0])


def ColorEdges(genome):
    return list(itertools.chain.from_iterable(
            EdgeHelper(chrom) for chrom in genome))
    

def Cycle2Chrom(cycle):
    return tuple(ChromElement(cycle[i], cycle[i+1])
                 for i in xrange(0, len(cycle), 2))


def EdgeHelper(chrom):
    cycle = list(Chrom2Cycle(chrom))
    cycle = cycle[1:] + [cycle[0]]
    return [(cycle[i], cycle[i+1]) for i in xrange(0, len(cycle), 2)]
        

def GetComplementPositions(kmer, complement_map):
    return itertools.chain.from_iterable(
        complement_map[k]
        for k in (kmer, ReverseComplement(kmer))
        if k in complement_map)


def Graph2Genome(edges):
    graph = {e[0]: e[1]
             for edge in edges
             for e in [edge, (edge[1], edge[0])]}
    this_edge = edges[0]
    cycle, results = [], []
    edge_set, visited = set(edges), set()
    while len(visited) < len(edges):
        cycle.extend(this_edge)
        candidate_edges = set([this_edge, (this_edge[1], this_edge[0])])
        real_edge_set = candidate_edges & edge_set
        if len(real_edge_set) != 1:
            raise Error('Invalid real_edge_set: {}.'.format(real_edge_set))
        visited.add(real_edge_set.pop())
        next_source = this_edge[1] + 1 if this_edge[1] % 2 else this_edge[1] - 1 
        if next_source == cycle[0]:
            cycle = [cycle[-1]] + cycle[:-1]
            results.append(Cycle2Chrom(cycle))
            diff_set = edge_set - visited
            if not diff_set: break
            this_edge = diff_set.pop()
            cycle = []
            continue
        this_edge = (next_source, graph[next_source])
    return results
    

def KmerMap(string, k):    
    kmap = collections.defaultdict(list)
    for i in xrange(len(string)-k+1):
        kmer = string[i:i+k]
        kmap[kmer].append(i)
    return kmap


def SharedKmers(k, string1, string2):
    complement_map = KmerMap(string2, k)
    return [(i, j)
            for i in xrange(len(string1)-k+1)
            for j in GetComplementPositions(string1[i:i+k], complement_map)]
        

class GreedySort(object):

    def __init__(self, permutation):
        self.permutation = permutation
        self.perm_map = {abs(int(e)): i for i, e in enumerate(permutation)}
        self.results = []

    def ChangeSign(self, element):
        value = -int(element)
        return '+{}'.format(value) if value > 0 else str(value)

    def Order(self, start):
        # Reverse elements and change signs.
        end = self.perm_map[start+1]
        if start != end:
            for i, element in enumerate(reversed(self.permutation[start:end+1])):
                index = start + i
                self.permutation[index] = self.ChangeSign(element)
                self.perm_map[abs(int(element))] = index
            self.results.append(self.PermutationString())
        # Check if the sign of this element needs to be reversed.
        element = self.permutation[start]
        if int(element) < 0:
            self.permutation[start] = self.ChangeSign(element)
            self.results.append(self.PermutationString())
        
    def PermutationString(self):
        return '({})'.format(' '.join(self.permutation))

    def Sort(self):
        for i in xrange(len(self.permutation)):
            self.Order(i)
        return self.results
 

class TwoBreakDistance(object):
 
    def __init__(self, p, q):
        self.p, self.q = p, q
        self.p_edges = ColorEdges(self.p)
        self.q_edges = ColorEdges(self.q)
        self.p_graph = self.GetGraph(self.p_edges)
        self.q_graph = self.GetGraph(self.q_edges)

    def GetCycle(self, p_candidates, q_candidates):
        p_edge = first_p = p_candidates.pop()
        cycle = [first_p]
        while p_candidates or q_candidates:
            q_edge = self.GetEdgeFromCandidates(p_edge, self.q_graph, q_candidates)
            cycle.append(q_edge)
            if q_edge[1] in first_p: return cycle
            p_edge = self.GetEdgeFromCandidates(q_edge, self.p_graph, p_candidates)
            cycle.append(p_edge)

    def GetDistance(self):
        p_candidates, q_candidates = set(self.p_edges), set(self.q_edges)
        cycles = []
        while p_candidates and q_candidates:
            cycles.append(self.GetCycle(p_candidates, q_candidates))
        return len(self.p_edges) - len(cycles)

    def GetEdgeFromCandidates(self, edge, graph, candidates):
        source = edge[1]
        dest = graph[source]
        edges = set([(source, dest), (dest, source)])
        edge_set = edges & candidates
        candidates -= edge_set
        return (source, dest)

    def GetGraph(self, edges):
        return dict(e
                    for edge in edges
                    for e in [edge, (edge[1], edge[0])])
            

class TwoBreakSort(TwoBreakDistance):

    def __init__(self, p, q):
        super(TwoBreakSort, self).__init__(p, q)
        self.p_edge_map = self.GetEdgeMap(self.p_edges)
        self.results = [p[:]]

    def Break(self, edge):
        pos0 = self.p_edge_map[edge[0]]
        pos1 = self.p_edge_map[edge[1]]
        if pos0 == pos1:
            raise Error('Invalid positions: {}, {} for edge: {}.'.format(
                    pos0, pos1, edge))
        # Perform break.
        other_edge = self.GetOtherEdge(pos0, pos1, edge)
        self.p_edges[pos0] = edge
        self.p_edges[pos1] = other_edge
        # Update new positions.
        self.p_edge_map[edge[0]] = self.p_edge_map[edge[1]] = pos0
        self.p_edge_map[other_edge[0]] = self.p_edge_map[other_edge[1]] = pos1
        # Store result.
        # import pdb; pdb.set_trace()
        self.results.append(Graph2Genome(self.p_edges))

    def GetEdgeMap(self, edges):
        return {node: index
                for index, edge in enumerate(edges)
                for node in edge}

    def GetOtherEdge(self, pos0, pos1, break_edge):
        edge_map = {e[0]: e[1]
                    for index, edge in enumerate(self.p_edges)
                    for e in [edge, (edge[1], edge[0])]
                    if index not in (pos0, pos1)}
        p_edge0 = self.p_edges[pos0]
        p_edge1 = self.p_edges[pos1]
        nodes = set(p_edge0 + p_edge1)
        other_edge_set = nodes - set(break_edge)
        this_node = break_edge[1]
        for _ in xrange(2*len(self.p_edges)):
            dest = this_node + 1 if this_node % 2 else this_node - 1
            if dest in other_edge_set:
                # dest should be the first element in the tuple.
                diff_set = other_edge_set - set([dest])
                return (dest, diff_set.pop())
            if dest == break_edge[0]:
                # A cycle is found.
                return tuple(sorted(other_edge_set, reverse=True))
            this_node = edge_map[dest]
        raise Error('Unable to return other tuple for {}.'.format(break_edge))

    def SkipTrivial(self):
        p_candidates, q_candidates = set(self.p_edges), set(self.q_edges)
        for edge in self.q_edges:
            edge_set = set([edge, (edge[1], edge[0])])
            p_candidates -= edge_set
            if len(p_candidates) < len(q_candidates):
                q_candidates -= edge_set
        return q_candidates.pop()

    def Sort(self):
        distance = self.GetDistance()
        for _ in xrange(self.GetDistance()):
            edge = self.SkipTrivial()
            self.Break(edge)
        return self.results
