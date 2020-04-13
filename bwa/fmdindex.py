from fmindex import * # get_bwt, get_sa, get_occ, get_C
from dataclasses import dataclass, replace
from itertools import tee, takewhile, groupby, accumulate
from operator import itemgetter
import matplotlib.pyplot as plt
"""
The adaption introuced in bwa mem, is to also include the reverse complement of the sequence in the FM index. 
And to find mathces for the query sequence and it's reverse compliment at the same time. 

This brings with it the property that there are exactly the same number of matches for a query sequence Q as for its reverse compliment. 
"""
compliments = str.maketrans("AGCT", "TCGA")
reverse_compliment = lambda seq: seq.translate(compliments)[::-1]

@dataclass
class BiInterval:
    start: int
    r_start: int
    length: int

@dataclass
class EM:
    start: int
    end: int
    bint: BiInterval

def lr_map(occ, C, bint, char):
    bint = BiInterval(0, 0, len(occ["A"])-1) if bint is None else bint
    new_start = C[char] + occ[char][bint.start]
    new_length = occ[char][bint.start+bint.length]-occ[char][bint.start]
    counts=[occ[c][bint.start+bint.length]-occ[c][bint.start] for c in "ACGT".translate(compliments)]
    new_r_start = bint.r_start+bint.length-sum(counts["ACGT".find(char.translate(compliments)):])
    return replace(bint, start=new_start, r_start=new_r_start, length=new_length)

def forward(lr, bint, c):
    r_bint = BiInterval(bint.r_start, bint.start, bint.length)
    new_bint = lr(r_bint, c.translate(compliments))
    return BiInterval(new_bint.r_start, new_bint.start, new_bint.length)

def _backward_search(lr, em, c):
    return replace(em, start=em.start-1, bint=lr(em.bint, c))

def _forward_search(forward_extend, em, c):
    return replace(em, end=em.end+1, bint=forward_extend(em.bint, c))

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

def _find_smem_candidates(backward_extend, query, position):
    init_bint = backward_extend(None, query[position])
    bints = accumulate(reversed(query[:position]), backward_extend, initial=init_bint)
    bints = takewhile(lambda bint: bint.length, bints)
    candidates = unique_lastseen(enumerate(bints), lambda cand: cand[1].length)
    return list(candidates)

def find_smem_candidates(backward_search, query, position):
    init_em = backward_search(EM(position+1, position+1, None), query[position])
    ems = accumulate(reversed(query[:position]), backward_search, initial=init_em)
    ems = takewhile(lambda em: em.bint.length, ems)
    return list(unique_lastseen(ems, lambda em: em.bint.length))

def _finalize_smem_candidates(forward_extend, query, position, candidates):
    def filter_candidates(candidates, c):
        print(candidates)
        extended = ((j, forward_extend(bint, c)) for (j, bint) in candidates)
        return [(j, bint) for j, bint in extended if bint.length]
    candidates = list(reversed(candidates))
    remaining_candidates = accumulate(query[position+1:], filter_candidates, initial=candidates)
    remaining_candidates = takewhile(lambda x: x, remaining_candidates)
    longest_candidates = ((*candidates[0], i) for i, candidates in enumerate(remaining_candidates))
    smems = unique_lastseen(longest_candidates, key=itemgetter(0))
    return [((position-j, position+i+1), bint) for j, bint, i in smems]

def finalize_smem_candidates(forward_search, query, position, candidates):
    def filter_candidates(candidates, c):
        extended = (forward_search(em, c) for em in candidates)
        return [em for em in extended if em.bint.length]
    candidates = list(reversed(candidates))
    print("#", candidates)
    remaining_candidates = accumulate(query[position+1:], filter_candidates, initial=candidates)
    remaining_candidates = takewhile(lambda x: x, remaining_candidates)
    longest_candidates = list(map(itemgetter(0), remaining_candidates))
    return unique_lastseen(longest_candidates, key=lambda em: em.start)
# return [((position-j, position+i+1), bint) for j, bint, i in smems]

def get_smems(lr, forward, query, position):
    candidates = list(find_smem_candidates(lr, query, position))
    print(candidates)
    return finalize_smem_candidates(forward, query, position, candidates)

def find_all_smems(backward, forward, query):
    find_smems = partial(get_smems, backward, forward,  query)
    i = 0
    all_smems = []
    while i<len(query):
        smems = list(find_smems(i))
        i = max(em.end for em in smems)
        all_smems.extend(smems)
    return all_smems

reference = "ATCGTTGTGC"
whole_sequence = reference+"$"+ reverse_compliment(reference)+"$"

bwt = get_bwt(whole_sequence)
occ = get_occ(bwt)
C = get_C(whole_sequence)

lr=partial(lr_map, occ, C)
forward_ext = partial(forward, lr)
backward_search=partial(_backward_search, lr)
forward_search=partial(_forward_search, forward_ext)

#smems = find_all_smems(lr, forward_ext, "ATCGTGGTGC")
#print(list(get_smems(backward_search, forward_search, "ATCGTGGTGC", 6)))
smems = find_all_smems(backward_search, forward_search, "ATCGTGGTGC")
suffixes = list(sorted(get_suffixes(whole_sequence)))
def plot_em(em):
    plot_bint(suffixes, em.bint, em.end-em.start)
for em in smems:
    plot_em(em)
