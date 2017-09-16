import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rnamake import sqlite_library
import rnamake.motif_type
import rnamake.motif_tree
import rnamake.cluster
import rnamake.motif


mlib = sqlite_library.MotifSqliteLibrary("nway")
mlib.load_all()
motifs = [m for m in mlib.all() if len(m.ends)==3]

# save motif sequences
junction_seqs = {}
for i, m in enumerate(motifs):
    seq_list = m.secondary_structure.sequence().split('&')
    if len(seq_list) == 3:
        junction_seqs[i] = pd.Series(seq_list + [m.name], index=['j1', 'j2', 'j3', 'tw_name'])
    else:
        print m.name
junction_seqs = pd.concat(junction_seqs).unstack()

# reduce duplicates
junction_seqs_red = {}
for name, group in junction_seqs.groupby('tw_name'):
    # take the one with the longest j1 sequence
    idx = group.j1.str.len().sort_values(ascending=False).index[0]
    junction_seqs_red[idx] = group.loc[idx]
junction_seqs_red = pd.concat(junction_seqs_red).unstack()

junction_seqs_red.to_csv('~/JunctionLibrary/seq_params/three_way_junctions.dat', sep='\t', index=False)
	
