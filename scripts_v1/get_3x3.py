import os
import numpy as np
import pandas as pd
import itertools
from fittinglibs import fileio
from tectolibs import tectplots
from hjh.junction import Junction
from hjh.helix import Helix


# make a bunch of AA mismatches (with different flanks)
junction_seqs = pd.concat([Junction(('U', 'C', 'M', 'M', 'M', 'G', 'U')).sequences,
                           Junction(('U', 'C', 'M', 'M', 'W', 'G', 'U')).sequences,
                          Junction(('U', 'C', 'W', 'W', 'W', 'G', 'U')).sequences,],
    ignore_index=True)
junction_seq_data = {}
for idx, junction_seq in junction_seqs.iterrows():
    junction_seq_abbrev = pd.Series({key:junction_seq.loc[key][2:-2] for key in ['side1', 'side2']}) 
    no_flank = '_'.join(junction_seq_abbrev)
    name = ';'.join([''.join(bases) for bases in zip(junction_seq_abbrev.side1, junction_seq_abbrev.side2[::-1])])
    junction_seq_data[idx] = pd.Series({'no_flank':no_flank, 'name':name})
junction_seq_data = pd.concat(junction_seq_data).unstack()

pd.concat([junction_seqs, junction_seq_data], axis=1).to_csv('~/JunctionLibrary/seq_params/3x3_junctions.dat', sep='\t', index=False)