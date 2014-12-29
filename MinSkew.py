INC_MAP = {'C': -1, 'G': 1}


def MinSkew(genome):
    skews = [0]
    min_skew = float('inf')
    for c in genome:
        skews.append(skews[-1] + INC_MAP.get(c, 0))
        min_skew = min(min_skew, skews[-1])
    return ' '.join(str(i) for i, skew in enumerate(skews)
                    if skew == min_skew)


def MaxSkew(genome):
    skews = [0]
    max_skew = -1 * float('inf')
    for c in genome:
        skews.append(skews[-1] + INC_MAP.get(c, 0))
        max_skew = max(max_skew, skews[-1])
    return ' '.join(str(i) for i, skew in enumerate(skews)
                    if skew == max_skew)
