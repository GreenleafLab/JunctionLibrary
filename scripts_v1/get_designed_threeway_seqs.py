import os
import pandas as pd
import itertools
import functools
from fittinglibs import seqfun
from hjh.junction import Junction
from hjh.helix import Helix

# bulge seq has two As
rc = functools.partial(seqfun.rc, rna=True)
bases = ['A', 'C', 'G', 'U']
bulge_seqs = ['AA', 'UU']
junction_seqs = {}
for bulge_seq in bulge_seqs:
    for (base1, base2, base3) in itertools.product(*([bases]*3)):
        junction = pd.Series({'j1': base1 + base2,
                              'j2': rc(base2) + bulge_seq + base3,
                              'j3': rc(base3) + rc(base1) })
        junction_seqs[(bulge_seq,  base1+base2+ base3)] = junction       
junction_seqs = pd.concat(junction_seqs, names=['bulge_seq', 'flank_seq']).unstack().reset_index()
junction_seqs.to_csv('~/JunctionLibrary/seq_params/three_way_designed.dat', sep='\t', index=False)


# cut out the 'extra'
junction_seq_controls = {}
for idx, row in junction_seqs.iterrows():
    for loop_context in ['L1','L2']:
        if loop_context == 'L1':
            two_way_seq = pd.Series({'side1':row.j1,
                                     'side2':row.j2[0] + row.j3[-1]})

        elif loop_context == 'L2':
            two_way_seq = pd.Series({'side1':row.j1[0] + row.j2[-1],
                                     'side2':row.j3})
        junction_seq_controls[(idx, loop_context)] = pd.concat([
                two_way_seq, row.drop(['j1', 'j2', 'j3'])])

        
junction_seq_controls = pd.concat(junction_seq_controls, names=['index', 'loop_context']).unstack().swaplevel(0,1).sort_index()
for loop_context in ['L1', 'L2']:
    junction_seq_controls.loc[loop_context].to_csv('~/JunctionLibrary/seq_params/three_way_designed_controls_%s.dat'%loop_context, sep='\t', index=False)
    
    
    
    