sequence = "ATTGTAGTTG$"
from itertools import accumulate, chain
from operator import itemgetter
from functools import partial
from array import array
from dataclasses import dataclass
alphabet = "ACGT"
@dataclass
class FMIndex:
    occ: dict
    C: dict
    sa: array

def get_suffixes(sequence):
    N = len(sequence)
    return ((2*sequence)[i:i+N] for i in range(N))

def get_bwt(sequence):
    suffixes = get_suffixes(sequence)
    return [suffix[-1] for suffix in sorted(suffixes)]

def argsort(iterable):
    return map(itemgetter(0), sorted(enumerate(iterable), key=itemgetter(1)))

def get_sa(sequence):
    suffixes = get_suffixes(sequence)
    return array("i", argsort(suffixes))

def get_occ(bwt):
    count = lambda c: array("i", chain([0], accumulate(int(e==c) for e in bwt)))
    return {c: count(c) for c in alphabet}

def get_C(sequence):
    offsets = accumulate(sequence.count(c) for c in sorted("$"+alphabet))
    return dict(zip(sorted(alphabet), counts))

def LF_map(occ, C, interval, char):
    return  (C[char] + occ[char][interval.start],
             C[char] + occ[char][interval.start+interval.length])

@dataclass
class SAInterval:
    start: int
    length: int

@dataclass
class ExactMatch:
    start: int
    end: int
    sa_int: SAInterval

def longest_subsequence(LF, sa, sequence):
    sa_intervals = accumulate(reversed(sequence), LF, initial=SAInterval(0, len(sa)))
    *_, longest = takewhile(lambda i: i.length, sa_intervals)
    return sa[longest:longest.start+last.length]

def get_fm_index(sequence):
    return FMIndex(get_occ(get_bwt(sequence)),
                   get_C(sequence),
                   get_sa(sequence))

if __name__ == "__main__":
    bwt = get_bwt(sequence)
    print("".join(bwt))
    print("Occ: ", get_occ(bwt))
    print("C: ", get_C(sequence))
    sa = get_sa(sequence)
    
    LF = partial(LF_map, get_occ(bwt), get_C(sequence))
    
    # interval = lr((0, len(sequence)), "T")
    # interval = lr(interval, "G")
    subseq = partial(longest_subsequence, LF, sa)
    print(subseq("TTG"))
