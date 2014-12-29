import collections
import itertools
import ReverseComplement


NUCLEOTIDES = set('ACGT')


class Error(Exception):
    """Base error."""


def HammingDistance(input1, input2, max_d=None):
    max_d = max_d or float('inf')
    input1_length, input2_length = len(input1), len(input2)
    min_length = min(input1_length, input2_length)
    result = abs(input1_length - input2_length)
    for i1, i2 in itertools.izip(input1, input2):
        if i1 != i2:
            result += 1
        if result > max_d:
            return result
    return result


def AppPatternMatch(pattern, text, d):
    pattern_length = len(pattern)
    return [str(i)
            for i in xrange(len(text)-pattern_length+1)
            if HammingDistance(pattern,
                               text[i:i+pattern_length],
                               max_d=d)
            <= d]


def AppPatternCount(pattern, text, d):
    return len(AppPatternMatch(pattern, text, d))


def FreqKmerMismatch(text, k, d):
    counter = collections.Counter()
    freq = 0
    for kmer in GetKmers(text, k):
        counter[kmer] = AppPatternCount(kmer, text, d)
        freq = max(freq, counter[kmer])
    # Get k-mers that have max frequency and can be found in text.
    kmers = {kmer for kmer, count in counter.iteritems()
             if count == freq}
    # Get other possible k-mers that cannot be found in text.
    possible_kmers = {kmer2 for kmer1 in kmers for kmer2 in GetPossibleKmers(kmer1, d)
                      if (kmer2 not in counter and
                          AppPatternCount(kmer2, text, d) == freq)}
    return sorted(kmers.union(possible_kmers))


def FreqKmerMismatchReverse(text, k, d):
    counter = collections.Counter()
    freq = 0
    for kmer in GetMismatchKmers(text, k, d):
        reverse_kmer = ReverseComplement.ReverseComplement(kmer)
        counter[kmer] = GetKmersReversePatternCount(kmer, text, d)
        freq = max(freq, counter[kmer])
    kmers = {kmer for kmer, count in counter.iteritems()
             if count == freq}
    return sorted(kmers)


def GetKmers(text, k):
    kmers = {text[i:i+k] for i in xrange(len(text)-k+1)}
    return kmers


def GetMismatchKmers(text, k, d):
    kmers = {text[i:i+k] for i in xrange(len(text)-k+1)}
    possible_kmers = {kmer2
                      for kmer1 in kmers
                      for kmer2 in GetPossibleKmers(kmer1, d)}
    return kmers.union(possible_kmers)


def GetKmersReversePatternCount(kmer, text, d):
    reverse_kmer = ReverseComplement.ReverseComplement(kmer)
    return AppPatternCount(kmer, text, d) + AppPatternCount(reverse_kmer, text, d)


def GetPossibleKmers(k, d):
    for index_tuple in itertools.combinations(xrange(len(k)), d):
        # Find all possible bases that can be changed in k-mer's positions.
        for base_tuple in itertools.product(NUCLEOTIDES, repeat=d): 
            if not OverlapBases(k, index_tuple, base_tuple):
                k_list = list(k)
                for index, base in itertools.izip(index_tuple, base_tuple):
                    k_list[index] = base
                yield ''.join(k_list)


def OverlapBases(k, index_tuple, base_tuple):
    for index, base in itertools.izip(index_tuple, base_tuple):
        if k[index] == base:
            return True
    return False


def NeighborsOld(pattern, d):
    if not d:
        return pattern
    if len(pattern) == 1:
        return NUCLEOTIDES
    first, suffix = pattern[0], pattern[1:]
    neighbors = set()
    suffix_neighbors = Neighbors(suffix, d)
    for n in suffix_neighbors:
        if HammingDistance(suffix, n) < d:
            neighbors.update(s + n for s in NUCLEOTIDES)
        else:
            neighbors.add(first + n)
    return neighbors


def Neighbors(pattern, d):
    if not d:
        return set([pattern])
    if len(pattern) == 1:
        return NUCLEOTIDES
    first, suffix = pattern[0], pattern[1:]
    neighbors = set()
    for n in Neighbors(suffix, d):
        neighbors.add(first + n)
    suffix_neighbors = Neighbors(suffix, d-1)
    for s in NUCLEOTIDES:
        neighbors.update(s + n for n in suffix_neighbors)
    return neighbors
