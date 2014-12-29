import collections


def FrequentWords(text, k):
    counter = collections.Counter()
    for i in xrange(len(text) - k):
        counter[text[i:i+k]] += 1
    max_value = max(counter.values())
    return {key for key, value in counter.iteritems()
            if value == max_value}


def PatternCount(text, pattern):
    count = 0
    pat_length = len(pattern)
    for i in xrange(len(text) - pat_length):
        if text[i:i+pat_length] == pattern:
            count += 1
    return count
        
