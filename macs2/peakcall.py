import math
import matplotlib.pyplot as plt
import heapq
from itertools import groupby, accumulate, chain
from operator import setitem, itemgetter
from scipy.stats import poisson

def parse_bed_line(line):
    _, start, end, _, _, strand = line.split(None, 6)[:6]
    return (int(start), int(end), strand)

def extend_read(read, fragment_length):
    start, end, strand = read
    start = start if strand=="+" else end-fragment_length
    end = end if strand=="-" else start+fragment_length
    return start, end

def get_sparse_pileup(fragments):
    steps = chain(*(((start, 1), (end, -1)) for start, end in fragments))
    return accumulate(sorted(steps), lambda x, y: (y[0], x[1]+y[1]))

def plot_step_function(steps, chrom_size=20000):
    plt.xlim(0, chrom_size)
    plt.step(*zip(*steps),  where="post")
    # plt.show()

def expand_to_window_size(read, window_size):
    start, end, strand = read
    midpoint = start if strand=="+" else end
    return midpoint-window_size//2, midpoint+window_size//2

def get_local_average(reads, window_size):
    covered_areas = (expand_to_window_size(read, window_size) for read in reads)
    steps = get_sparse_pileup(covered_areas)
    return apply_func(steps, lambda v: v/window_size)

def uniq_indexes(pileup):
    return (last_value for index, (*_, last_value) in groupby(pileup, itemgetter(0)))

def uniq_values(pileup):
    return (next(group) for _, group in groupby(pileup, itemgetter(1)))


def sync_indexes(steps_a, steps_b):
    update = lambda v, tv: (v[0], tv[1]) if tv[0] else (tv[1], v[1])
    all_steps = heapq.merge(

def sync_indexes(pileups):
    cur_values = [0 for _ in pileups]
    updated_values = lambda j, v: setitem(cur_values, j, v) or cur_values
    all_steps = ([(step, j) forv step in steps]  for j, steps in enumerate(pileups))
    steps = ((i, updated_values(j, v)) for (i, v), j in heapq.merge(*all_steps))
    return uniq_indexes(steps)

def apply_func(steps, func):
    return [(i, func(v)) for i, v in steps]

from itertools import dropwhile
def get_peaks(scores, threshold=2):
    thresholded = apply_func(scores, lambda x: x>-math.log(threshold, 10))
    cleaned = dropwhile(lambda iv: not iv[1], uniq_values(thresholded))
    return [i for i, v in cleaned]

def paired(iterator):
    return zip(*([iter(iterator)]*2))

def merge_peaks(peaks, max_gap=10):
    holes = paired(peaks[1:-1])
    filtered_holes = chain(*((end, start) for end, start in holes if start-end>=max_gap))
    return list(chain(peaks[:1], filtered_holes, peaks[-1:]))

def remove_small_peaks(peaks, min_size=100):
    return [(start, end) for start, end in paired(peaks) if end-start>=min_size]

def plot_peaks(peaks, y=-1):
    for start, end in peaks:
        plt.hlines(y, start, end)
###
treat_reads = [parse_bed_line(line) for line in open("data/treat.bed")]
fragments = [extend_read(read, 100) for read in treat_reads]
treat_pileup = list(uniq_indexes(get_sparse_pileup(fragments)))
plot_step_function(treat_pileup)
input_reads = [parse_bed_line(line) for line in open("data/input.bed")]
sparse_pileups = [get_local_average(input_reads, w) for w in (200, 400)]
control_lambda = apply_func(sync_indexes(sparse_pileups), max)
scale_factor = len(treat_reads)/len(input_reads)*100
scaled_lambda = apply_func(control_lambda, lambda v: v*scale_factor)
plot_step_function(scaled_lambda)
get_p_score = lambda pair: -poisson.logsf(pair[0]-1, pair[1])/math.log(10)
log_p = apply_func(sync_indexes((treat_pileup, scaled_lambda)), get_p_score)
peaks = get_peaks(log_p, 0.001)
plot_peaks(paired(peaks))
merged_peaks = merge_peaks(peaks)
plot_peaks(paired(merged_peaks), -2)
final = remove_small_peaks(merged_peaks)
plot_peaks(final, -3)
plt.show()
