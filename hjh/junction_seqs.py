import hjh.mutations
import pandas as pd
import numpy as np
import pickle

def getSingleMutationIres():
    pdb_junctions = pd.read_table('150311_library_v2/ires_junctions.txt', index_col=0)
    otherSide = {'side1':'side2', 'side2':'side1'}
    # make single mutants of all
    pdb_junctions_mut = {}
    for ires in pdb_junctions.index:
        key = ires+'_1m'
        pdb_junctions_mut[key] = {}
        for side in ['side1', 'side2']:
            s = ''
            
            # find seqs and strip flanks
            side_x = pdb_junctions.loc[ires, side]
            side_y = pdb_junctions.loc[ires, otherSide[side]]
            flank = [side_x[i] for i in [0,-1]]
            no_flank_x = ''.join([side_x[i] for i in np.arange(1, len(side_x)-1)])
            no_flank_y = ''.join([side_y[i] for i in np.arange(1, len(side_y)-1)])
            
            # initialize data frame to store all mutations
            pdb_junctions_mut[key][side] = {}
            num = len(no_flank_x)*3 + 1
            pdb_junctions_mut[key][side][s] = pd.DataFrame(index=np.arange(num),
                                                         columns = ['side1', 'side2', 'n_flank', 'seq'])
            pdb_junctions_mut[key][side][s].loc[:, 'n_flank'] = 1
            
            # find single mutations of side X
            pdb_junctions_mut[key][side][s].loc[:, side]            = hjh.mutations.singles(no_flank_x) 
            pdb_junctions_mut[key][side][s].loc[:, otherSide[side]] =  hjh.mutations.complement(flank[-1]) + no_flank_y + hjh.mutations.complement(flank[0])
            
            # add flankers back in
            pdb_junctions_mut[key][side][s].loc[:, side] =  flank[0] + pdb_junctions_mut[key][side][s].loc[:, side] + flank[-1]
            #pdb_junctions_mut[ires][side][s].loc[:, otherSide[side]] = hjh.mutations.complement(flank[-1]) + pdb_junctions_mut[ires][side][s].loc[:, otherSide[side]] + hjh.mutations.complement(flank[0])
            
            # concatenate into data fram
            pdb_junctions_mut[key][side] = pd.concat(pdb_junctions_mut[key][side], names=['single_ind', 'junction_ind'])
        pdb_junctions_mut[key] = pd.concat(pdb_junctions_mut[key], names=['mut_side'])
        
    pdb_junctions_mut = pd.concat(pdb_junctions_mut, names=['junction'])
    pdb_junctions_mut.loc[:, 'seq'] = ['_'.join(pdb_junctions_mut.loc[i, ['side1', 'side2']] )for i in pdb_junctions_mut.index]
    pdb_junctions_mut.drop_duplicates(subset=['seq'], inplace=True)
    return pdb_junctions_mut

def loadAllseqs(filenames, load_pickle=None):
    if load_pickle is None: load_pickle is False
    allSeqs = []
    # load allSeqs
    for filename in filenames:
        allSeqSub = pd.read_table(filename+'.txt', index_col=0)
        allSeqSub.loc[:, 'tecto_object'] = np.nan
        if load_pickle:
            try:
                with open(filename+'.pkl', 'rb') as input:
                    for loc in allSeqSub.index:
                        allSeqSub.loc[loc, 'tecto_object'] = pickle.load(input)
            except IOError: pass
        allSeqs.append(allSeqSub)
    allSeqs = pd.concat(allSeqs, ignore_index=True)
    return allSeqs