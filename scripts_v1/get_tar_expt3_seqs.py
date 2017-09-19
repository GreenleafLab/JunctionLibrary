import os
import pandas as pd
import itertools
from hjh.junction import Junction
from hjh.helix import Helix

def find_seq_and_data(bulge_seq, helix_sub):
    """Assemble the seqs"""
    junction_len = 0 # don't replace anything in the flaning helix
    h1_length = 1 # put junction 1 bp into flanking helix
    total_length = 3 # total length of flaning helix
    helix = Helix(helix_sub, junction_len)
    helix_df = helix.formatHelix(helix_sub, total_length, h1_length)
    
    # put bulge seq just on one side
    seq = pd.Series({'side1':helix_df.before.side1 + bulge_seq + helix_df.after.side1,
                     'side2':helix_df.before.side2 + helix_df.after.side2})
    name = bulge_seq + ':' + ';'.join([''.join(bases) for bases in zip(helix_sub.side1, helix_sub.side2[::-1])])
    seq_data = pd.Series({'j_name':name, 'bulge_seq':bulge_seq, 'flank':'_'.join(helix_sub)})
    return pd.concat([seq, seq_data])    

# dinfe the bulge sequences
starting_bulge_seqs = (['UCU'] + # wt seq
    ['U'*i for i in range(1, 5)] + #
    ['C'*i for i in range(1, 5)] +
    [''])
flanking_bps = Junction(('W+', 'W+', 'A')).sequences
tar_seqs = []
for bulge_seq, (idx, helix_sub) in itertools.product(starting_bulge_seqs, flanking_bps.iterrows()):
    tar_seqs.append(find_seq_and_data(bulge_seq, helix_sub))

# do the same but for WT sequences also vary the top bp
bulge_seq = ['UCU']
flanking_bps_top = Junction(('W', 'W', 'W+')).sequences
for bulge_seq, (idx, helix_sub) in itertools.product([starting_bulge_seqs[0]], flanking_bps_top.iterrows()):
    tar_seqs.append(find_seq_and_data(bulge_seq, helix_sub))
    
# save
tar_seqs = pd.concat(tar_seqs, axis=1).transpose()
tar_seqs_red = tar_seqs.set_index(['side1', 'side2']).groupby(level=[0,1]).first().reset_index()