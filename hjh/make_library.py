#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

# This function calls dependent functions to make the subset of
# final library that is all three by three junctions and less


##### IMPORT #####
import numpy as np
import pandas as pd
import os
import sys
import argparse
import itertools
from operator import mul
import subprocess
import cPickle as pickle
from matplotlib import gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import functools
import multiprocessing

# load custom libraries
import hjh.tecto_assemble
from hjh.helix import Helix
from hjh.junction import Junction
from hjh.tecto_assemble import TectoSeq


#set up command line argument parser
parser = argparse.ArgumentParser(description="script for making library")
parser.add_argument('-map','--map_file', help='file that contains parameters to vary', required=True)
parser.add_argument('-max','--max_number', help='maximum number of seuqences per juncito motif', type=int)
parser.add_argument('-par','--seq_params_file', help='location of seq_param file. Default is ~/JunctionLibrary/seq_params/seq_params.txt')
parser.add_argument('-jun','--junction_sequence_file', help='overwrite the junction motifs given in map with particular sequenecs in this file')

parser.add_argument('-out','--out_file', help='file to save output')

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()
    
def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def findInitialDistance(allSeqs, indx):
    seq, ss = findSecondaryStructure(allSeqs.loc[indx])
    script = '~/RNAMake/rnamake/lib/RNAMake/simulate_tectos/test_chip_seq.py'
    call = 'python %s -cseq "%s" -css "%s"'%(script, seq,ss )
    distance=np.nan
    try:
        result = subprocess.check_output(call, shell=True).strip().split('\n')
        if result[-2].find('WARNING') != 0:
            distance = float(result[-1])
    except:
        pass
    print '%d\t%s\t%s\t%4.2f'%(allSeqs.loc[indx, 'offset'], allSeqs.loc[indx, 'junction_seq'], call, distance)
    return pd.Series([allSeqs.loc[indx, 'offset'], allSeqs.loc[indx, 'junction_seq'],distance, seq], index=['offset', 'junction_seq', 'distance', 'tecto_sequence'])



#parse command line arguments
args = parser.parse_args()
if args.out_file is None:
    if args.junction_sequence_file is None:
        args.out_file = (os.path.splitext(args.map_file)[0] + '.library').replace(',', '')
    else:
        args.out_file = (os.path.splitext(args.junction_sequence_file)[0] + '.library').replace(',', '')
if args.seq_params_file is None:
    args.seq_params_file = os.path.expanduser('~/JunctionLibrary/seq_params/seq_params.txt')

# load file setting parameters like loop and receptor sequences
expt_params = pd.read_table(args.map_file)
seq_params    = pd.read_table(args.seq_params_file , index_col=[0,1])
numCores = 20

# if seq file is defined, replace junction motifs with this
if args.junction_sequence_file is not None:
    expt_params.loc[:, 'junction'] = ['defunct'] +[np.nan]*(len(expt_params)-1)
    junctionSeqs = pd.read_table(args.junction_sequence_file, index_col=0)
else:
    junctionSeqs = pd.DataFrame(columns=['side1', 'side2'])
    for motif in expt_params.loc[:, 'junction']:
        junctionSeqs = junctionSeqs.append(Junction(tuple(motif.split(','))).sequences, ignore_index=True)
    expt_params.loc[:, 'junction'] = ['defunct'] +[np.nan]*(len(expt_params)-1)


# iterate through all
cols = TectoSeq(seq_params, expt_params.loc[0]).params.index
numToLog =  np.product((~expt_params.isnull()).sum().values)
allSeqs = pd.DataFrame(columns=cols)
logSeqs = pd.DataFrame(index=np.arange(numToLog), columns=expt_params.columns.tolist()+['number'])

# initiate structure to store, per motif, and many contexts, whether it could form
with open(args.out_file + '.pkl', 'wb') as output:
    count = 0
    for i, params in enumerate(itertools.product(*[expt_params.loc[:,name].dropna() for name in expt_params])):
        pass
        params = pd.Series(params, index=expt_params.columns.tolist())
      
        # if up, junction remains the same. else switch sides
        junction = Junction(sequences=junctionSeqs.loc[:, ['side1', 'side2']])
        if params.side == 'up':
            pass
        elif params.side == 'down':
            junction.sequences = pd.DataFrame(junction.sequences.loc[:, ['side2', 'side1']].values,
                                                         columns=['side1', 'side2'],
                                                         index=[-name for name in junction.sequences.index])
        
        # initialize saving
        allSeqSub = pd.DataFrame(index=junction.sequences.index, columns=cols)
        logSeqs.loc[i] = params
        logSeqs.loc[i, 'number'] = len(junction.sequences)
        print '%4.1f%% complete'%(100*i/float(numToLog))
        # cycle through junction sequences
        f = functools.partial(hjh.tecto_assemble.findTecto, params, junction, seq_params)
        workerPool = multiprocessing.Pool(processes=numCores)
        allSeqSub = workerPool.map(f, np.array_split(junction.sequences.index.tolist(), numCores))
        workerPool.close(); workerPool.join()
        
        allSeqs = allSeqs.append(pd.concat(allSeqSub), ignore_index=True)

# check if successful secondary structure        
allSeqs.loc[:, 'ss'] = hjh.tecto_assemble.getAllSecondaryStructures(allSeqs.loc[:, 'tecto_sequence'])

#check if any module was successful
numberFlankingBasepairs = 2
allSeqs.loc[:, 'numberFlanking'] = numberFlankingBasepairs
allSeqs.loc[:, 'no_flank'] = ['_'.join([allSeqs.loc[loc, 'tecto_object'].junction['side1'][numberFlankingBasepairs:-numberFlankingBasepairs],
                                       allSeqs.loc[loc, 'tecto_object'].junction['side2'][numberFlankingBasepairs:-numberFlankingBasepairs]]) for loc in allSeqs.index] 
allSeqs.loc[:, 'flank'] = ''
for loc in allSeqs.index:
    if allSeqs.loc[loc, 'side'] == 'up':
        side = 'side1'
    else:
        side = 'side2'
    allSeqs.loc[loc, 'flank'] = allSeqs.loc[loc, 'tecto_object'].junction[side][:numberFlankingBasepairs] + allSeqs.loc[loc, 'tecto_object'].junction[side][-numberFlankingBasepairs:]

workerPool = multiprocessing.Pool(processes=numCores)
indices = np.array_split(allSeqs.index.tolist(), numCores)
allSeqSub = workerPool.map(hjh.tecto_assemble.getSecondaryStructureMultiprocess,
                           [allSeqs.loc[index] for index in indices])
workerPool.close(); workerPool.join()

allSeqs = pd.concat(allSeqSub) 
allSeqs.drop(['tecto_object'], axis=1).to_csv(args.out_file+'.txt', sep='\t')

flanks = np.unique(allSeqs.loc[:, 'flank'])
no_flanks = [f for f in np.unique(allSeqs.loc[:, 'no_flank']) if f[0] != '_']
iterables = [no_flanks, expt_params.length.dropna(), expt_params.offset.dropna().values, expt_params.side.dropna().values]
index = pd.MultiIndex.from_product(iterables, names=['no_flank', 'length', 'offset', 'side'])
allsuccess = pd.DataFrame(index=flanks, columns=index)

for loc in allSeqs.index:
    #if loc%10 == 0: print loc
    index = allSeqs.loc[loc, 'flank']
    no_flank = allSeqs.loc[loc, 'no_flank']
    if allSeqs.loc[loc, 'side'] == 'down':
        no_flank = no_flank[::-1]
    if loc%100==0:
        print loc, no_flank
    
    if no_flank in no_flanks:
        columns = (no_flank,
                   allSeqs.loc[loc, 'length'],
                   allSeqs.loc[loc, 'offset'],
                   allSeqs.loc[loc, 'side'])
        allsuccess.loc[index, columns] = allSeqs.loc[loc, 'ss_correct']
allsuccess.to_csv(args.out_file+'.success.mat', sep='\t')

allsuccessTrimmed = allsuccessAll.dropna()
worksInAll = pd.Series(index=allsuccessTrimmed.index, data=[allsuccessTrimmed.loc[loc].values.astype(bool).sum() for loc in allsuccessTrimmed.index])
worksInAll.sort()

no_flanks = ['A_', 'G_', 'C_', 'U_'] + [''.join(x)+'_' for x in itertools.product('ACGU', 'ACGU')]
worksInAll = pd.DataFrame(index=allsuccessTrimmed.index, columns = no_flanks)
for seq in no_flanks:
    worksInAll.loc[:, seq]=[allsuccessTrimmed.loc[loc, seq].values.astype(bool).sum() for loc in allsuccessTrimmed.index]

plt.figure(figsize=(20,20))
a = worksInAll.loc[(worksInAll >=12).all(axis=1), 'A_'].copy()
a.sort()
ax = sns.heatmap(worksInAll.loc[a.index], square=True,
                    linewidth=0.1, vmax=24, vmin=0,
                    annot=True, fmt='d',
                    cbar=False)
ax.set_yticklabels(a.index.tolist(), rotation=0, fontsize=12)

# plot those that work in two sets
b = worksInAll.loc[a.index].loc[(worksInAll.loc[a.index, ['A_', 'U_']] == 24).all(axis=1)].sum(axis=1).copy()
b.sort()
c = worksInAll.loc[a.index].loc[(worksInAll.loc[a.index, ['G_', 'C_']] == 24).all(axis=1)].sum(axis=1).copy()
c.sort()
plt.figure(figsize=(10,10))
ax = sns.heatmap(worksInAll.loc[pd.concat([b,c]).index], square=True,
                    linewidth=0.1, vmax=24, vmin=0,
                    annot=True, fmt='d',
                    cbar=False)
ax.set_yticklabels(pd.concat([b,c]).index.tolist()[::-1], rotation=0, fontsize=12)
ax.set_xticklabels(no_flanks, rotation=90, fontsize=12)
plt.savefig(args.out_file + '.work_in_AU_or_CG.pdf')
# plot those that work in other two sets
b = worksInAll.loc[a.index].loc[(worksInAll.loc[a.index, ['A_', 'C_']] == 24).all(axis=1)].sum(axis=1).copy()
b.sort()
c = worksInAll.loc[a.index].loc[(worksInAll.loc[a.index, ['G_', 'U_']] == 24).all(axis=1)].sum(axis=1).copy()
c.sort()
plt.figure(figsize=(10,10))
ax = sns.heatmap(worksInAll.loc[pd.concat([b,c]).index], square=True,
                    linewidth=0.1, vmax=24, vmin=0,
                    annot=True, fmt='d',
                    cbar=False)
ax.set_yticklabels(pd.concat([b,c]).index.tolist()[::-1], rotation=0, fontsize=12)
ax.set_xticklabels(no_flanks, rotation=90, fontsize=12)
plt.savefig(args.out_file + '.work_in_AC_or_GU.pdf')



