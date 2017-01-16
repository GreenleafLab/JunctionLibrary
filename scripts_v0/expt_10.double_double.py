import sys
import math
import random
import os
import hjh.tecto_assemble
import pandas as pd
import numpy as np
  
saveDir = '150311_library_v2/double_double'
if not os.path.exists(saveDir): os.mkdir(saveDir)

seq_params_file = os.path.expanduser('~/JunctionLibrary/seq_params/seq_params.txt')
seq_params    = pd.read_table(seq_params_file , index_col=[0,1])

seqs = pd.read_table('150311_library_v2/double_double_sequences.txt', names=['seq'], header=0, sep=' ')

params = pd.Series({'helix':np.nan,
'junction':'double_double',
'receptor':       '11nt',
'loop':           'GGAA',
'base':         'normal',
'adapters':     'truseq',
'length':           10,
'offset':            0,
'side':             'up'})
tectoSeq = hjh.tecto_assemble.TectoSeq(seq_params, params)
allSeqs = pd.DataFrame(index=np.arange(len(seqs)), columns=tectoSeq.params.index.tolist() + ['ss'])

for i in seqs.index:
    
    allSeqs.loc[i, params.index] = params
    allSeqs.loc[i, 'junction'] += '_%d'%i
    allSeqs.loc[i, 'tecto_sequence'] = seqs.loc[i, 'seq'].replace('T', 'U')
    seq = seq_params.loc[('adapters', params.loc['adapters'])]
    allSeqs.loc[i, 'sequence'] = seq.loc['side1'] +  seqs.loc[i, 'seq'].replace('U', 'T') + seq.loc['side2']

out_file = os.path.join(saveDir, 'all.sequences.library')
   
neworder = [ 'junction', 'length', 'offset', 'helix', 'receptor', 'loop',  'side',
            'flank', 'no_flank', 'n_flank', 'junction_seq',  'junction_SS',  'ss_correct',
            'adapters', 'base', 'helix_one_length', 'effective_length', 'helix_seq',
            'tecto_sequence', 'sequence', 'ss']
allSeqs.loc[:, neworder].to_csv(out_file+'.txt', sep='\t')



