sequence = "ATTGTAGTTG$"
from itertools import accumulate, chain
from functools import partial

def get_suffixes(sequence):
    N = len(sequence)
    return ((2*sequence)[i:i+N] for i in range(N))

def get_bwt(sequence):
    suffixes = get_suffixes(sequence)
    return [suffix[-1] for suffix in sorted(suffixes)]

def argsort(iterable):
    return (i for v, i in sorted((v, i) for i, v in enumerate(iterable)))

def get_sa(sequence):
    suffixes = get_suffixes(sequence)
    return list(argsort(suffixes))

def get_occ(bwt):
    count = lambda c: [0]+list(accumulate(int(e==c) for e in bwt))
    return {c: count(c) for c in "AGCT"}

def get_C(sequence):
    counts = accumulate(sum(e==c for e in sequence) for c in sorted("$AGCT"))
    return dict(zip(sorted("AGCT"), counts))

def lr_map(occ, C, interval, char):
    return  (C[char] + occ[char][interval[0]],
             C[char] + occ[char][interval[1]])

def longest_subsequence(lr, sa, sequence):
    interval = (0, len(sa))
    for c in reversed(sequence):
        old_interval, interval = interval, lr(interval, c)
        if interval[0]==interval[1]:
            break
    else:
        old_interval=interval
    return sa[slice(*old_interval)]


if __name__ == "__main__":
    bwt = get_bwt(sequence)
    print("".join(bwt))
    print("Occ: ", get_occ(bwt))
    print("C: ", get_C(sequence))
    sa = get_sa(sequence)
    
    lr = partial(lr_map, get_occ(bwt), get_C(sequence))
    
    # interval = lr((0, len(sequence)), "T")
    # interval = lr(interval, "G")
    subseq = partial(longest_subsequence, lr, sa)
    print(subseq("TTG"))
