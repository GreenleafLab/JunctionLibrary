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
    junctionMotifs = []
    junctionSeqs = {}
    flanks = [['G', 'C', 'G', 'C']]

    maxNumSeqs = 12
    for motif in [ '_', 'B1', 'B1,B1', 'B1,B1,B1', 'M', 'M,M', 'M,M,M', 'M,B1', 'M,M,B1', 'M,B1,B1']:
        junctionSeqs[motif] = {}
        
        for flank in flanks:

            baseNum = len(flank)/2
            junctionMotif = ','.join(flank[:baseNum] + [motif] + flank[baseNum:])
    
            junctionSeq = Junction(tuple(junctionMotif.split(','))).sequences
            junctionSeq.loc[:, 'n_flank'] = baseNum
            numSeqs = len(junctionSeq)
            
            # reduce total number of sequences 
            if numSeqs > maxNumSeqs:
                index = np.linspace(0, numSeqs - 1, maxNumSeqs).astype(int)
                
                junctionSeq = junctionSeq.loc[index]
                
            junctionSeqs[motif][''.join(flank)] = junctionSeq
        junctionSeqs[motif] = pd.concat(junctionSeqs[motif], names=['flank', 'junction_num'])
    
    return pd.concat(junctionSeqs, names=['junction'])

## for junction conformation expt, do two flanking base pairs, 7 different positions
if __name__ == '__main__':
 
    saveDir = '150311_library_v2/tertiary_contacts/'
    if not os.path.exists(saveDir): os.mkdir(saveDir)
    
    # get sequences
    junctionSeqs = getAllJunctionSeqs()
    junctionSeqs.to_csv(os.path.join(saveDir, 'all.junctions_to_compare.junctions'), sep='\t', index=True)

    # load receptor/loop combos
    filenames = {}
    for s in ['loops', 'receptors', 'receptors_abbrev']:
        filenames[s] = os.path.join(saveDir, 'expt_map.%s.txt'%s)
    
    for key, filename in filenames.items():
        print "%%run ~/JunctionLibrary/hjh/make_library.py -map %s -jun %s -o %s -con"%(filename,
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.junctions'),
                                                                                     os.path.join(saveDir, 'all.junctions_to_compare.%s.junctions'%key))
    
    sys.exit()  # and run make_library commands
    
    
    filenames = [os.path.join(saveDir, 'all.junctions_to_compare.%s.junctions'%s) for s in ['loops', 'receptors','receptors_abbrev']]
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
    subset = allSeqs.loc[(allSeqs.loc[:, 'junction'] == '_')&(allSeqs.loc[:, 'length']==10)]
    for ind in subset.index:
        name = os.path.join(saveDir, '_'.join(subset.loc[ind, ['receptor', 'loop']])+'.png')
        cmnd = subset.loc[ind, 'tecto_object'].printVarnaCommand(name=name)
        os.system(cmnd)