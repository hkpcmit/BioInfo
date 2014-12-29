COMPLEMENT_MAP = {'A': 'T',
                  'C': 'G',
                  'G': 'C',
                  'T': 'A'}


def ReverseComplement(input):
    return ''.join(COMPLEMENT_MAP[c] for c in reversed(input))


def StartPositionPatternMatch(pattern, genome):
    pattern_length = len(pattern)
    return [i for i in xrange(len(genome)-pattern_length)
            if genome[i:i+pattern_length] == pattern]

