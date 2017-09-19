import os
import numpy as np
import pandas as pd
import itertools
from fittinglibs import fileio
from tectolibs import tectplots
from hjh.junction import Junction
from hjh.helix import Helix

def find_seq_and_data(junction_seq, helix_sub):
    """Assemble the seqs"""
    junction_len = 0 # don't replace anything in the flaning helix
    total_length = len(helix_sub.side1) # total length of flaning helix
    h1_length = int(np.floor(total_length/2.)) # put junction 1 bp into flanking helix

    helix = Helix(helix_sub, junction_len)
    helix_df = helix.formatHelix(helix_sub, total_length, h1_length)
    
    # put bulge seq just on one side
    seq = pd.Series({'side1':helix_df.before.side1 + junction_seq.side1 + helix_df.after.side1,
                     'side2':helix_df.before.side2 + junction_seq.side2 + helix_df.after.side2})
    name = '_'.join(junction_seq) + ':' + ';'.join([''.join(bases) for bases in zip(helix_sub.side1, helix_sub.side2[::-1])])
    seq_data = pd.Series({'j_name':name})
    return pd.concat([seq, seq_data])    



# make a bunch of AA mismatches (with different flanks)
flanking_bps = pd.concat([Junction(('C', 'M', 'C', 'G')).sequences,
                          Junction(('C', 'C', 'M', 'G')).sequences,
                          Junction(('W', 'W', 'W', 'W')).sequences])

# first do single AA mismatch
junction_seqs = pd.concat([pd.Series({'side1':'A', 'side2':'A'}),
                           pd.Series({'side1':'A', 'side2':'U'})], axis=1).transpose()
AA_seqs = []
for (idx1, junction_seq), (idx, helix_sub) in itertools.product(junction_seqs.iterrows(), flanking_bps.iterrows()):
    AA_seqs.append(find_seq_and_data(junction_seq, helix_sub))

# save
AA_seqs = pd.concat(AA_seqs, axis=1).transpose()
AA_seqs.to_csv('~/JunctionLibrary/seq_params/AA_junctions.dat', sep='\t', index=False)