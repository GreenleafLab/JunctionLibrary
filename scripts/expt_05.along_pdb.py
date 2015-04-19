from hjh.junction import Junction
import itertools
import os
import pandas as pd
import hjh.mutations
from hjh.junction import Junction
import numpy as np
import sys

# functions
def getAllJunctionSeqs():
    # save junctions : first B1 junctions
    # also get PDB sequences
    pdb_junctions = pd.read_table('150311_library_v2/all_pdb_junctions.txt', index_col=0, names=['side1', 'side2', 'type', 'ss'], header=None)
    pdb_junctions.loc[:, 'ss_side1'] = [pdb_junctions.loc[i, 'ss'].split('&')[0] for i in pdb_junctions.index]
    pdb_junctions.loc[:, 'ss_side2'] = [pdb_junctions.loc[i, 'ss'].split('&')[1] for i in pdb_junctions.index]
    pdb_junctions.loc[:, 'length_diff'] = [len(pdb_junctions.loc[i, 'side1']) - len(pdb_junctions.loc[i, 'side2'])  for i in pdb_junctions.index]

    # flip everything that's longer on one side than the other
    index = pdb_junctions.loc[:, 'length_diff'] < 0
    side1 = pdb_junctions.loc[index, 'side2']
    side2 = pdb_junctions.loc[index, 'side1']
    ss_side1 = [pdb_junctions.loc[i, 'ss_side2'].replace(')', '(') for i in index.loc[index].index]
    ss_side2 = [pdb_junctions.loc[i, 'ss_side1'].replace('(', ')') for i in index.loc[index].index]
    pdb_junctions.loc[index, ['side1']] = side1
    pdb_junctions.loc[index, ['side2']] = side2
    pdb_junctions.loc[index, ['ss']] = ['&'.join(s) for s in itertools.izip(ss_side1, ss_side2)]
    
    # filter out junctions that aren't properly formatted
    index = [pdb_junctions.loc[i, 'ss_side1'][0] == '.' or pdb_junctions.loc[i, 'ss_side1'][-1] == '.'
             or pdb_junctions.loc[i, 'ss_side2'][0] == '.' or pdb_junctions.loc[i, 'ss_side2'][-1] == '.' for i in pdb_junctions.index]
    pdb_junctions = pdb_junctions.loc[np.logical_not(index)]
    
    # uniqueify
    pdb_junctions.loc[:, 'junction_seq'] = ['_'.join(pdb_junctions.loc[i, ['side1', 'side2']]) for i in pdb_junctions.index]
    pdb_junctions.drop_duplicates(subset='junction_seq', inplace=True)
    
    # add flanking base pair
    for i in pdb_junctions.index:
        pdb_junctions.loc[i, 'side1'] = 'C' + pdb_junctions.loc[i, 'side1'] + 'G'
        pdb_junctions.loc[i, 'side2'] = 'C' + pdb_junctions.loc[i, 'side2'] + 'G'
        s = pdb_junctions.loc[i, 'ss'].split('&')
        pdb_junctions.loc[i, 'ss'] = '('+ s[0] + '(' + '&' + ')' + s[1] + ')'
        
    # find flank
    pdb_junctions.loc[:, 'flank'] = ''
    pdb_junctions.loc[:, 'n_flank'] = 2
    for i in pdb_junctions.index:
        ind = np.in1d(list(pdb_junctions.loc[i, 'ss'].split('&')[0]), '(')
        pdb_junctions.loc[i, 'flank'] = ''.join(np.array(list(pdb_junctions.loc[i, 'side1']))[ind])
    
    pdb_junctions.loc[:, 'junction'] = ''
    for i in pdb_junctions.index:
        s = pdb_junctions.loc[i, 'ss'].split('&')
        pdb_junctions.loc[i, 'junction'] = 'pdb_' + 'x'.join([str(np.sum(np.in1d(list(x), '.'))) for x in s])

    return pdb_junctions.loc[:, ['junction', 'flank', 'side1', 'side2', 'n_flank']].sort('junction')

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/junction_conformations_pdb/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)
    
    # make map files
    lengths = [8, 9, 10, 11, 12]
    offsetsPerLength = {8:[0], 9:[-1,0,1], 10:[-1,0,1], 11:[0], 12:[0] }
    filenames = {}
    num_lengths = len(lengths)
    expt_map_constant = pd.Series(index = ['helix','junction','receptor','loop', 'base',   'adapters'],
                                  data  = ['wc',   'defunct', '11nt',    'GGAA', 'normal', 'truseq'], dtype=str)
    cols = expt_map_constant.index.tolist()
    sides = ['up', 'down']
    for length in lengths:
        offsets = offsetsPerLength[length]
        num = max(len(offsets), len(sides))
        expt_map = pd.DataFrame(index = np.arange(num),
                            columns = cols + ['length','offset','side'])
        expt_map.loc[0, cols] = expt_map_constant
        expt_map.loc[0, 'length'] = length
        expt_map.loc[np.arange(len(offsets)), 'offset'] = offsets
        expt_map.loc[np.arange(len(sides)), 'side'] = sides
        filenames[length] = os.path.join(saveDir, 'expt.length_%d.map'%length)
        expt_map.to_csv(filenames[length], index=False, sep='\t')
    
    filenames = {}
    for length in lengths:
        filenames[length] = os.path.join(saveDir, 'expt.length_%d.map'%length)  
        print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s -out %s"%(filenames[length],
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.length_%d'%length))
    
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