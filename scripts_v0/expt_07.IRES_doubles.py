from hjh.junction import Junction
import itertools
import os
import pandas as pd
import hjh.mutations
from hjh.junction import Junction
import numpy as np
import sys
import hjh.junction_seqs

# functions
def getAllJunctionSeqs():

    otherSide = {'side1':'side2', 'side2':'side1'}
    # make single mutants of all
    pdb_junctions_mut = hjh.junction_seqs.getSingleMutationIres()
    
    pdb_junctions_2mut = {}
    for ires in ['HCV', 'SVV']:
        key = ires+'_2m'
        pdb_junctions_2mut[key] = {}
        for side in ['side1', 'side2']:
            pdb_junctions_2mut[key][side] = {}


            for i in pdb_junctions_mut.loc[ires+'_1m'].index:
                s = '%d:%s'%(i[-1], i[0])
                
                # find seqs and strip flanks
                side_x = pdb_junctions_mut.loc[ires+'_1m'].loc[i, side]
                side_y = pdb_junctions_mut.loc[ires+'_1m'].loc[i,  otherSide[side]]
                flank = [side_x[i] for i in [0,-1]]
                no_flank_x = ''.join([side_x[i] for i in np.arange(1, len(side_x)-1)])
                no_flank_y = ''.join([side_y[i] for i in np.arange(1, len(side_y)-1)])
                
                # initialize data structure
                num = len(no_flank_x)*3 + 1
                pdb_junctions_2mut[key][side][s] = pd.DataFrame(index=np.arange(num),
                                                           columns = ['side1', 'side2', 'n_flank', 'seq'])
                pdb_junctions_2mut[key][side][s].loc[:, 'n_flank'] = 1
                
                # save all single mutations
                pdb_junctions_2mut[key][side][s].loc[:, side]            = hjh.mutations.singles(no_flank_x)
                pdb_junctions_2mut[key][side][s].loc[:, otherSide[side]] = hjh.mutations.complement(flank[-1]) + no_flank_y +  hjh.mutations.complement(flank[0])
                
                # add back flanker
                pdb_junctions_2mut[key][side][s].loc[:, side] = flank[0] + pdb_junctions_2mut[key][side][s].loc[:, side] + flank[-1]
                
            pdb_junctions_2mut[key][side] = pd.concat(pdb_junctions_2mut[key][side], names=['single_ind', 'junction_ind'])
            pdb_junctions_2mut[key][side].loc[:, 'seq'] = ['_'.join(pdb_junctions_2mut[key][side].loc[i, ['side1', 'side2']] )for i in pdb_junctions_2mut[key][side].index]
        pdb_junctions_2mut[key] = pd.concat(pdb_junctions_2mut[key], names=['mut_side'])
    pdb_junctions_2mut = pd.concat(pdb_junctions_2mut, names=['junction'])        
    
    # check uniqueness
    pdb_junctions_2mut.drop_duplicates(subset=['seq'], inplace=True)
    
    return pdb_junctions_2mut

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/ires_doubles/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()

    # downsample by half
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files
    lengths = [ 9, 10, 11]
    sides = ['up', 'down']

    num_lengths = len(lengths)
    expt_map_constant = pd.Series(index = ['junction','receptor', 'helix','loop', 'base',   'adapters', 'offset'],
                                  data  = ['defunct', '11nt',  'wc',  'GGAA', 'normal', 'truseq', 0], dtype=str)
    cols = expt_map_constant.index.tolist()


    expt_map = pd.DataFrame(index = np.arange(num_lengths),
                        columns = cols + ['length'])
    expt_map.loc[0, cols] = expt_map_constant
    expt_map.loc[np.arange(len(sides)), 'side'] = sides
    expt_map.loc[np.arange(len(lengths)), 'length'] = lengths
    filename = os.path.join(saveDir, 'expt.length.map')
    expt_map.to_csv(filename, index=False, sep='\t')

    print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s -ss"%(filename,
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     )
    
    sys.exit()  # and run make_library commands
    
    filenames = [os.path.join(saveDir, 'all.junctions_to_compare.length_%d'%length) for length in lengths]
    allSeqs = []
    # load allSeqs
    for filename in filenames:
        allSeqSub = pd.read_table(filename+'.txt', index_col=0)
        allSeqSub.loc[:, 'tecto_object'] = np.nan
        with open(filename+'.pkl', 'rb') as input:
            for loc in allSeqSub.index:
                allSeqSub.loc[loc, 'tecto_object'] = pickle.load(input)
        allSeqs.append(allSeqSub)
    allSeqs = pd.concat(allSeqs, ignore_index=True)
    
    # print varna sequences fro moving A bulge and think about it
    subset = allSeqs.loc[(allSeqs.loc[:, 'no_flank'] == 'A_')&(allSeqs.loc[:, 'flank'] == 'GCGC')&(allSeqs.loc[:, 'side'] == 'up')]
    for loc in subset.index:
        cmnd = subset.loc[loc, 'tecto_object'].printVarnaCommand()[0]
        os.system(cmnd.replace('test.png', 'junction_%d.png'%loc))