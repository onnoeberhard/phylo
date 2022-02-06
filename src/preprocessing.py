"""
Turn the raw DNA sequences into distance matrices. Two matrices are created and
saved: 'dat/K_len.npy' and 'dat/K_lev.npy', where 'len' is the "length"-distance 
and 'lev' the Levenshtein distance. In 'K_len', many some distances are NaN.

Onno Eberhard, Jan 2022
"""

import numpy as np
import pandas as pd
from Levenshtein import distance as levenshtein
from tqdm.auto import trange


def len_distance(a, b):
    return abs(len(a) - len(b))

# Load data
species = pd.read_csv('dat/data.csv')

# Build matrix of length-distances
K_len = np.full((len(species), len(species)), np.nan)
for i in trange(len(species)):
    a = species.dna[i]
    for j in range(i, len(species)):
        b = species.dna[j]
        K_len[i, j] = K_len[j, i] = len_distance(a, b)

# Build partial matrix of Leveshtein distance, where length-distance < 2000
# This step took 24h on a 2015 2.2 GHz Quad-Core Intel Core i7 MacBook Pro
K_lev = np.full((len(species), len(species)), np.nan)
for i in trange(len(species)):
    a = species.dna[i]
    K_lev[i, i] = 0
    for j in trange(i + 1, len(species), leave=False):
        b = species.dna[j]
        if K_len[i, j] < 2000 and (len(a) + len(b)) // 2 < 40_000:
            K_lev[i, j] = K_lev[j, i] = levenshtein(a, b)

# Save matrices
np.save('dat/K_len.npy', K_len)
np.save('dat/K_lev.npy', K_lev)
