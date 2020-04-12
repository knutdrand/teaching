from fmindex import * # get_bwt, get_sa, get_occ, get_C
from dataclasses import dataclass
from itertools import tee, takewhile, groupby, accumulate
from operator import itemgetter
import matplotlib.pyplot as plt
"""
The adaption introuced in bwa mem, is to also include the reverse complement of the sequence in the FM index. 
And to find mathces for the query sequence and it's reverse compliment at the same time. 

This brings with it the property that there are exactly the same number of matches for a query sequence Q as for its reverse compliment. 
"""
@dataclass
class BiInterval:
    start: int
    r_start: int
    length: int

class EM:
    start: int
    end: int
    bint: BiInterval

compliments = str.maketrans("AGCT", "TCGA")
reference = "ATCGTTGTGC"

def reverse_compliment(sequence):
    return sequence.translate(compliments)[::-1]

def get_C(sequence):
    counts = accumulate(sum(e==c for e in sequence) for c in sorted("AGCT"))
    return dict(zip(sorted("AGCT"), chain([2], (c+2 for c in counts))))

whole_sequence = reference+"$"+ reverse_compliment(reference)+"$"

bwt = get_bwt(whole_sequence)
occ = get_occ(bwt)
C = get_C(whole_sequence)

def lr_map(occ, C, bint, char):
    bint = BiInterval(0, 0, len(occ["A"])-1) if bint is None else bint
    new_start = C[char] + occ[char][bint.start]
    new_length = occ[char][bint.start+bint.length]-occ[char][bint.start]
    counts=[occ[c][bint.start+bint.length]-occ[c][bint.start] for c in "ACGT".translate(compliments)]
    new_r_start = bint.r_start+bint.length-sum(counts["ACGT".find(char.translate(compliments)):])
    return BiInterval(new_start, new_r_start, new_length)

def plot_bint(lines, bint, L):
    print(bint)
    lines = list(list(line) for line in lines)
    M, N = len(lines), len(lines[0])
    for i, line in enumerate(reversed(lines)):
        for j, c in enumerate(line):
            plt.text(j+1, i+1, c, ha="center", va="center")
    start = M-bint.start-bint.length+0.5
    r_start = M-bint.r_start-bint.length+0.6
    plt.hlines((start, start+bint.length), 0, N+1, color="green")
    plt.hlines((r_start, r_start+bint.length), 0, N+1, color="red")
    plt.vlines(L+0.5, 0, M+1)
    plt.xlim(0, M+1)
    plt.ylim(0, N+1)
    plt.show()

def unique_lastseen(iterable, key=None):
    return (last for _, (*_, last) in groupby(iterable, key))

def find_smem_candidates(backward_extend, query, position):
    print(position)
    init_bint = backward_extend(None, query[position])
    bints = accumulate(reversed(query[:position]), backward_extend, initial=init_bint)
    bints = takewhile(lambda bint: bint.length, bints)
    candidates = unique_lastseen(enumerate(bints), lambda cand: cand[1].length)
    return list(candidates)

def finalize_smem_candidates(forward_extend, query, position, candidates):
    def filter_candidates(candidates, c):
        extended = ((j, forward_extend(bint, c)) for (j, bint) in candidates)
        return [(j, bint) for j, bint in extended if bint.length]
    candidates = list(reversed(candidates))
    remaining_candidates = accumulate(query[position+1:], filter_candidates, initial=candidates)
    remaining_candidates = takewhile(lambda x: x, remaining_candidates)
    longest_candidates = ((*candidates[0], i) for i, candidates in enumerate(remaining_candidates))
    smems = unique_lastseen(longest_candidates, key=itemgetter(0))
    return [((position-j, position+i+1), bint) for j, bint, i in smems]

def forward(lr, bint, c):
    r_bint = BiInterval(bint.r_start, bint.start, bint.length)
    new_bint = lr(r_bint, c.translate(compliments))
    return BiInterval(new_bint.r_start, new_bint.start, new_bint.length)

def get_smems(lr, forward, query, position):
    candidates = find_smem_candidates(lr, query, position)
    return finalize_smem_candidates(forward, query, position, candidates)

def find_all_smems(backward, forward, query):
    i = 0
    all_smems = []
    while i<len(query):
        smems = get_smems(backward, forward, query, i)
        i = max(end for (start, end), bint in smems)
        all_smems.extend(smems)
    return all_smems
lr=partial(lr_map, occ, C)
forward_ext = partial(forward, lr)

smems = find_all_smems(lr, forward_ext, "ATCGTGGTGC")
for (start, end), bint in smems:
    plot_bint(list(sorted(get_suffixes(whole_sequence))), bint, end-start)
