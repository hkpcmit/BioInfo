import collections
import itertools
import random
from HammingDistance import GetPossibleKmers, HammingDistance


random.seed(7971)


def CheckKmersInDna(dna, kmer, d):
    k = len(kmer)
    for i in xrange(len(dna)-k+1):
        candidate = dna[i:i+k]
        if HammingDistance(kmer, candidate, max_d=d) <= d:
            return True
    return False


def DistancePatternDna(pattern, dna):
    k = len(pattern)
    return min(HammingDistance(pattern, dna[i:i+k])
               for i in xrange(len(dna)-k+1))


def DistancePatternStrings(pattern, dna_list):
    return sum(DistancePatternDna(pattern, dna) for dna in dna_list)


def GetMismatchKmers(kmer, k, d):
    kmers = set([kmer])
    kmers.update(ki
                 for i in range(1, d+1)
                 for ki in GetPossibleKmers(kmer, i))
    return kmers


def GetMotifKmers(dna, k, d):
    kmers = {dna[i:i+k] for i in xrange(len(dna)-k+1)}
    possible_kmers = {k2
                      for k1 in kmers
                      for i in range(1, d+1)
                      for k2 in GetPossibleKmers(k1, i)}
    return kmers.union(possible_kmers)


def GetProfileColumnProbability(col_counter):
    total = float(sum(col_counter.itervalues()))
    return [col_counter.get(c, 0) / total
            for c in 'ACGT']


def GetProfileKmerProbability(kmer, profile_matrix):
    return reduce(lambda x, y: x*y,
                  (profile_matrix[c][i] for i, c in enumerate(kmer)))


def GetProfileList(motifs, laplace=False):
    motif_list = [[c for c in m]
                  for m in motifs]

    def _GetCounter(col):
        if laplace:
            col = list(col) + ['A', 'C', 'G', 'T']
        return collections.Counter(col)

    col_counter_list = [_GetCounter(col)
                        for col in itertools.izip(*motif_list)]
    prob_list = [GetProfileColumnProbability(ct)
                 for ct in col_counter_list]
    return zip(*prob_list)


def GetRandomMotifs(dna_list, k):

    def _GetRandomKmer(dna):
        i = random.randint(0, len(dna)-k)
        return dna[i:i+k]

    return [_GetRandomKmer(dna) for dna in dna_list]


def GibbsRandomProfileKmer(text, k, profile_list):
    profile_matrix = ProfileMatrix(profile_list)
    prob_list = [GetProfileKmerProbability(text[i:i+k], profile_matrix)
                 for i in xrange(len(text)-k+1)]
    sum_prob = sum(prob_list)
    cum_prob_list = [prob_list[0] / sum_prob]
    for p in prob_list[1:]:
        cum_prob_list.append(cum_prob_list[-1] + p/sum_prob)
    rand_num = random.random()
    cum_prob_length = len(cum_prob_list)
    for i, cp in enumerate(cum_prob_list):
        if rand_num < cp:
            return text[i:i+k]


def GibbsSampler(dna_list, k, t, N):
    motifs = GetRandomMotifs(dna_list, k)
    best_motifs = motifs[:]
    for _ in xrange(N):
        i = random.randint(0, t-1)
        motifs2 = [m for l, m in enumerate(motifs) if l != i]
        profile = GetProfileList(motifs2, laplace=True)
        motifs[i] = GibbsRandomProfileKmer(dna_list[i], k, profile)
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs[:]
    return best_motifs


def GreedyMotifSearch(k, t, dna_list, laplace=False):
    best_motifs = [dna[:k] for dna in dna_list]
    for j in xrange(len(dna_list[0])-k+1):
        kmer = dna_list[0][j:j+k]
        motifs = [kmer]
        for i in xrange(1, t):
            profile = GetProfileList(motifs[:i], laplace=laplace)
            motifs.append(ProfileMostProbableKmer(dna_list[i], k, profile))
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
    return best_motifs


def MedianString(k, dna_list):
    k_list = sorted((DistancePatternStrings(pattern, dna_list), pattern)
                    for tu in itertools.product('ACGT', repeat=k)
                    for pattern in [''.join(tu)])
    min_score = k_list[0][0]
    results = [k_list[0][1]] 
    for tu in k_list:
        if min_score < tu[0]:
            return results
        results.append(tu[1])


def MotifEnumeration(dna_list, k, d):
    kmers = GetMotifKmers(dna_list[0], k, d)
    return {kmer
            for kmer in kmers
            if all(CheckKmersInDna(dna, kmer, d)
                   for dna in dna_list[1:])}


def ProfileMatrix(profile_list):
    return {c: profile_list[i] for i, c in enumerate('ACGT')}


def ProfileMostProbableKmer(text, k, profile_list):
    profile_matrix = ProfileMatrix(profile_list)

    def _kmer_prob(i):
        kmer = text[i:i+k]
        return GetProfileKmerProbability(kmer, profile_matrix), kmer

    # _, max_kmer = max(_kmer_prob(i) for i in xrange(len(text)-k+1))
    max_prob, max_kmer = -1, ''
    for i in xrange(len(text)-k+1):
        prob, kmer = _kmer_prob(i)
        if prob > max_prob:
            max_prob, max_kmer = prob, kmer
    return max_kmer


def ProfileProbability(kmer, profile_list):
    profile_matrix = ProfileMatrix(profile_list)
    return reduce(lambda x, y: x*y,
                  (profile_matrix[c][i] for i, c in enumerate(kmer)))


def Score(motifs):
    motif_c_list = [[c for c in m] for m in motifs]

    def _GetColumnScore(col):
        ct = collections.Counter(col)
        _, most_freq_item = max((freq, item) for item, freq in ct.iteritems())
        return len(col) - ct[most_freq_item]

    return sum(_GetColumnScore(col)
               for col in itertools.izip(*motif_c_list))


def RandomizedMotifSearch(k, t, dna_list):
    motifs = GetRandomMotifs(dna_list, k)
    best_motifs = motifs
    while True:
        profile = GetProfileList(motifs, laplace=True)
        motifs = [ProfileMostProbableKmer(dna, k, profile)
                  for dna in dna_list]
        if Score(motifs) >= Score(best_motifs):
            return best_motifs
        best_motifs = motifs


def RepeatRandomizedMotifSearch(k, t, dna_list, trials):            
    best_motifs = None
    for _ in xrange(trials):
        motifs = RandomizedMotifSearch(k, t, dna_list)
        if (not best_motifs) or (Score(motifs) < Score(best_motifs)):
            best_motifs = motifs
    return best_motifs
