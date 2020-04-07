def parse_bed_line(line):
    chrom, start, end, _, _, strand = line.split(None, 6)[:6]
    return (chrom, int(start), int(end), strand)
reads = (parse_bed_line(line) for line in open("data/test.bed"))

def extend_read(read, fragment_length):
    chrom, start, end, strand = read
    start = start if strand=="+" else end-fragment_length
    end = end if strand=="-" else start+fragment_length
    return chrom, start, end

fragments = [extend_read(read, 296) for read in reads]
import numpy as np
import matplotlib.pyplot as plt

from itertools import chain
def get_sparse_pileup(fragments):
    index_diff_pairs = chain(*([(start, 1), (end, -1)] for _, start, end in fragments))
    indexes, diffs = zip(*sorted(index_diff_pairs))
    return indexes, np.cumsum(diffs)


treat_pileup = get_sparse_pileup(fragments)

def plot_sparse_pileup(indexes, values, chrom_size=500):
    plt.xlim(0, chrom_size);
    plt.step(indexes, values, where="post")
    plt.show()
plot_sparse_pileup(*treat_pileup)

def expand_to_window_size(read, window_size):
    chrom, start, end, strand = read
    midpoint = start if strand=="+" else end
    return chrom, midpoint-window_size//2, midpoint+window_size//2
input_reads = (parse_bed_line(line) for line in open("data/test.bed"))
window_size = 500
covered_areas = (expand_to_window_size(read, window_size) for read in input_reads)

indexes, counts = get_sparse_pileup(covered_areas)
averages = counts/window_size

plot_sparse_pileup(indexes, averages)

from collections import namedtuple
SparsePileup = namedtuple("SparsePileup", ["Indexes", "Values"])
def get_local_average(reads, window_size):
    covered_areas = (expand_to_window_size(read, window_size) for read in reads)
    indexes, counts = get_sparse_pileup(covered_areas)
    return SparsePileup(indexes, counts/window_size)

import heapq
from itertools import groupby
from operator import setitem, itemgetter
def apply_func(sparse_pileups, func): 
     cur_values = [0 for _ in sparse_pileups] 
     updated_values = lambda vps: [setitem(cur_values, p, v) for _, v, p in vps] and cur_values 
     key_values = lambda p, sp: ((i, v, p) for (i, v) in zip(*sp)) 
     all_key_values = [key_values(*pair) for pair  in enumerate(sparse_pileups)] 
     index_newvalue_pairs = [(i, func(*updated_values(vps))) for i, vps in groupby(heapq.merge(*all_key_values), key=itemgetter(0))] 
     return SparsePileup(*zip(*index_newvalue_pairs))
reads = [parse_bed_line(line) for line in open("data/test.bed")]
sparse_pileups = [get_local_average(reads, w) for w in (100, 200)]
[plt.step(*sp, where="post") for sp in sparse_pileups]; plt.xlim(0, 500);plt.show()
control_lambda = apply_func(sparse_pileups, max)
plt.xlim(0, 500);plt.step(*control_lambda, where="post")

from scipy.stats import poisson
get_p_score = lambda count, l: -poisson.logsf(count, l)/np.log(10)
log_p = apply_func((treat_pileup, control_lambda), get_p_score)
plt.xlim(0, 500);plt.step(*log_p, where="post")

def get_peaks(scores, threshold=2):
    value_grouped = groupby(zip(*scores), lambda pair: pair[1]>-np.log10(threshold))
    indexes, values = zip(*((next(pairs)[0],v) for v, pairs in value_grouped))
    return indexes, values
peaks = get_peaks(log_p, 0.0000001)
plt.xlim(0, 500);plt.step(*peaks, where="post")

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def merge_peaks(peaks, max_gap):
    return SparsePileup(*zip(*((i, v) for (i, v), (j, _) in pairwise(zip(*peaks)) if not(j-i<=max_gap and v))))

def remove_small_peaks(pekas, min_size):
    return SparsePileup(*zip(*((i, v) for (i, v), (j, _) in pairwise(zip(*peaks)) if not(j-i<=max_gap and v))))
