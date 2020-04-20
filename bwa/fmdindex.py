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
class ExactMatch:
    start: int
    end: int
    bint: BiInterval

def count_occurances(occ, char, length, start=0):
    return occ[char][start+length]-occ[char][start]

def get_LF_map(fm):
    N = len(fm.occ["A"])-1
    def count_occ(char, length, start=0):
        return fm.occ[char][start+length]-fm.occ[char][start]

    def LF_map(bint, char):
        bint = BiInterval(0, 0, N) if bint is None else bint
        new_start = fm.C[char] + count_occ(char, bint.start)
        new_length = count_occ(char, bint.length, bint.start)
        counts=[count_occ(c, bint.length, bint.start) for c in "ACGT".translate(compliments)]
        new_r_start = bint.r_start+bint.length-sum(counts["ACGT".find(char.translate(compliments)):])
        return replace(bint, start=new_start, r_start=new_r_start, length=new_length)

    return LF_map

def forward_extend(backward_extend, bint, c):
    r_bint = replace(bint, start=bint.r_start, r_start=bint.start)
    new_bint = backward_extend(r_bint, c.translate(compliments))
    return replace(new_bint, start=new_bint.r_start, r_start=new_bint.start)

def _backward_search(LF, em, c):
    return replace(em, start=em.start-1, bint=LF(em.bint, c))

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
    return (last for *_, last in map(itemgetter(1), groupby(iterable, key)))

def find_smem_candidates(backward_search, query, position):
    init_em = backward_search(ExactMatch(position+1, position+1, None), query[position])
    ems = accumulate(reversed(query[:position]), backward_search, initial=init_em)
    ems = takewhile(lambda em: em.bint.length, ems)
    return list(unique_lastseen(ems, lambda em: em.bint.length))

def finalize_smem_candidates(forward_search, query, position, candidates):
    def filter_candidates(candidates, c):
        extended=map(partial(forward_search, c=c), candidates)
        return [em for em in extended if em.bint.length]
    remaining_candidates = accumulate(query[position+1:], filter_candidates, initial=candidates)
    remaining_candidates = takewhile(lambda x: x, remaining_candidates)
    longest_candidates = map(itemgetter(-1), remaining_candidates)
    return unique_lastseen(longest_candidates, key=lambda em: em.start)

def get_smems(LF, forward, query, position):
    candidates = list(find_smem_candidates(LF, query, position))
    return finalize_smem_candidates(forward, query, position, candidates)

def find_all_smems(backward, forward, query):
    find_smems = partial(get_smems, backward, forward,  query)
    i = 0
    all_smems = []
    while i < len(query):
        smems = list(find_smems(i))
        i = smems[-1].end  # max(em.end for em in smems)
        all_smems.extend(smems)
    return all_smems

def get_smem_finder(LF_map):
    backward_search=partial(_backward_search, LF)
    forward_ext = partial(forward_extend, LF)
    forward_search=partial(_forward_search, forward_ext)
    return partial(find_all_smems, backward_search, forward_search)

def translate_bint(sa, bint):
    return sa[bint.start:bint.start+bint.length]

def plot_em_matches(ref, query, sa, em):
    ref_starts = translate_bint(sa, em.bint)
    for j, c in enumerate(ref):
        plt.text(j, 1, c, ha="center", va="center")
    for j, c in enumerate(query):
        plt.text(j, -1, c, ha="center", va="center")
    L = em.end-em.start
    plt.hlines(range(2, 2+len(ref_starts)), ref_starts, [s+L for s in ref_starts])
    plt.hlines(-2, em.start, em.end)
    plt.xlim(0, len(ref)+1)
    plt.show()

reference = "ATCGTTGTGC"
whole_sequence = reference+"$"+ reverse_compliment(reference)+"$"
fm_index = get_fm_index(whole_sequence)
LF=get_LF_map(fm_index) # partial(LF_map, fm_index.occ, C)
smem_finder = get_smem_finder(LF_map)
#
#forward_ext = partial(forward, LF)
#backward_search=partial(_backward_search, LF)
#forward_search=partial(_forward_search, forward_ext)
#
##smems = find_all_smems(LF, forward_ext, "ATCGTGGTGC")
##print(list(get_smems(backward_search, forward_search, "ATCGTGGTGC", 6)))
query = "ATCGTGGTGC"
smems = smem_finder(query)
suffixes = list(sorted(get_suffixes(whole_sequence)))
def plot_em(em):
    plot_bint(suffixes, em.bint, em.end-em.start)
for em in smems:
    plot_em_matches(whole_sequence, query, fm_index.sa, em)
    # plot_em(em)
