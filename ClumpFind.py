import collections


def ClumpFind(genome, k, L, t):
    # Get list of starting positions for each k-mer.
    position_map = collections.defaultdict(list)
    genome = genome.upper()
    for i in xrange(len(genome)-k+1):
        position_map[genome[i:i+k]].append(i)
    return [k_mer for k_mer, position_list in position_map.iteritems()
            if IsClump(position_list, L, t)]

def IsClump(position_list, L, t):
    list_length = len(position_list)
    if list_length < t:
        return False
    for i, position in enumerate(position_list):
        t_index = i + t - 1
        if t_index < list_length:
           t_position = position_list[t_index]
           if (t_position - position) < L:
               return True
    return False
