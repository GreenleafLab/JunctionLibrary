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
if args.out_file is None: args.out_file = (os.path.splitext(args.map_file)[0] + '.library').replace(',', '')
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
junction = Junction(sequences=junctionSeqs.loc[:, ['side1', 'side2']])

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

workerPool = multiprocessing.Pool(processes=numCores)
indices = np.array_split(allSeqs.index.tolist(), numCores)
allSeqSub = workerPool.map(hjh.tecto_assemble.getSecondaryStructureMultiprocess,
                           [allSeqs.loc[index] for index in indices])
workerPool.close(); workerPool.join()

allSeqs = pd.concat(allSeqSub) 
allSeqs.drop(['tecto_object'], axis=1).to_csv(args.out_file+'.txt', sep='\t')



#check if any module was successful
allSeqs.loc[:, 'no_flank'] = ['_'.join([allSeqs.loc[loc, 'tecto_object'].junction['side1'][2:-2],
                                       allSeqs.loc[loc, 'tecto_object'].junction['side2'][2:-2]]) for loc in allSeqs.index] 
flanks = np.unique([allSeqs.loc[loc, 'tecto_object'].junction['side1'][:2] + allSeqs.loc[loc, 'tecto_object'].junction['side1'][-2:] for loc in allSeqs.index])

iterables = [['A_', 'G_', 'C_', 'U_', '_'], expt_params.length.dropna(), expt_params.offset.dropna().values, expt_params.side.dropna().values]
index = pd.MultiIndex.from_product(iterables, names=['no_flank', 'length', 'offset', 'side'])
allsuccess = pd.DataFrame(index=flanks, columns=index)

for loc in allSeqs.index:
    if loc%10 == 0: print loc
    index = allSeqs.loc[loc, 'tecto_object'].junction['side1'][:2] + allSeqs.loc[loc, 'tecto_object'].junction['side1'][-2:]
    no_flank = allSeqs.loc[loc, 'no_flank']
    if allSeqs.loc[loc, 'side'] == 'down':
        no_flank = no_flank[::-1]
    columns = (no_flank,
               allSeqs.loc[loc, 'length'],
               allSeqs.loc[loc, 'offset'],
               allSeqs.loc[loc, 'side'])
    allsuccess.loc[index, columns] = allSeqs.loc[loc, 'ss_correct']
sys.exit()

worksInAll = pd.DataFrame(index = ['_'.join(x) for x in itertools.product(['GC', 'CG', 'UA', 'AU'], ['GC', 'CG', 'UA', 'AU'])],
                          columns=[allsuccess.index.tolist()])
for flank1, flank2 in itertools.product(['GC', 'CG', 'UA', 'AU'], ['GC', 'CG', 'UA', 'AU']):
    for ind in allsuccess.index.tolist():
        worksInAll.loc['%s_%s'%(flank1, flank2), ind] = allsuccess.loc[ind,(flank1, flank2)].dropna().values.astype(bool).sum()

numPlots = np.max([1, np.around(worksInAll.shape[1]/50.)])
inds = np.array_split(np.arange(worksInAll.shape[1]), numPlots)
for i, ind in enumerate(inds):
    plt.figure(figsize=(0.16*(len(ind)-4)+4, 5))
    sns.heatmap(worksInAll.iloc[:,ind].astype(int), square=True,
                        annot=True, fmt="d",
                        cbar=False,)
    plt.tight_layout()
    plt.savefig('%s.successes.%d.png'%(args.out_file, i))

sys.exit()

# plot somehow
plt.figure(figsize=(6,10))
sns.heatmap(allsuccess.astype(int), square=True,
                    
                    cbar=False,)



for flank1, flank2 in itertools.product(['GC', 'CG', 'UA', 'AU'], ['GC', 'CG', 'UA', 'AU']):
    for i, length in enumerate([8.0, 9., 10., 11.]):
        
        print '%s\t%s'%(flank1, flank2,)
        print allsuccess.loc[:,(flank1, flank2, '%3.1f'%length)].sum(axis=1)
    
    offsets = [-2., -1., 0., 1., 2.]
    fig = plt.figure(figsize=(5,7))
    
    gs = gridspec.GridSpec(4, 2, wspace=0.1, hspace=0.1)
    for i, length in enumerate([8.0, 9., 10., 11.]):
        
        print '%s\t%s'%(flank1, flank2,)
        print allsuccess.loc[:,(flank1, flank2, '%3.1f'%length)].sum(axis=1)
        
        ax1 = fig.add_subplot(gs[i,0])
        ax2 = fig.add_subplot(gs[i,1])
        
        if i == 0:
            ax1.set_title('side1')
            ax2.set_title('side2')
            ax1.annotate('%s %s'%(flank1, flank2), xy=(1, 1.1),
                xytext=(0.8, 0.95), textcoords='figure fraction',
                horizontalalignment='center', verticalalignment='top',
                )
            
        if i == 3:
            xticklabels=offsets
            ax1.set_xlabel('offsets')
            ax2.set_xlabel('offsets')
        else:
            xticklabels=['']
            ax1.set_xlabel('')
            ax2.set_xlabel('')
        
        sns.heatmap(allsuccess.loc[:,(flank1, flank2, '%3.1f'%length,)].iloc[:, np.arange(0, len(offsets)*2, 2)],
                    ax=ax1,
                    square=True,
                    annot=True, fmt="d",
                    cbar=False,
                    yticklabels=allsuccess.index.tolist(),
                    xticklabels=xticklabels)
        
        sns.heatmap(allsuccess.loc[:,(flank1, flank2, '%3.1f'%length,)].iloc[:, np.arange(1, len(offsets)*2, 2)],
                    ax=ax2,
                    square=True,
                    annot=True, fmt="d",
                    cbar=False,
                    xticklabels=xticklabels)


# module success
iterables = [['GC', 'CG', 'UA', 'AU'], ['GC', 'CG', 'UA', 'AU'], np.array([8.0, 9., 10., 11.,]).astype(str),
    np.array([-2., -1., 0., 1., 2.]).astype(str), ['none', 'swap']]
index = pd.MultiIndex.from_product(iterables, names=['flank1', 'flank2', 'length', 'offset', 'side'])
allsuccess = pd.DataFrame(index = np.unique(allSeqs.loc[:, 'junction_seq_noflank']), columns=index)


allsuccess.to_csv(filename, sep='\t', index=True, header=True)

# laod again
num = logSeqs.number.sum()
tectoSeqs = ['']*num
with open(args.out_file, 'rb') as input:
    for j in range(num):

        tectoSeqs[j] = pickle.load(input)
        seq, ss = tectoSeqs[j].findSecondaryStructure(tectoSeq=tectoSeqs[j].params.tecto_sequence)
        print '%d\t%s\t%s\t%s'%(j, tectoSeqs[j].params.junction_seq, ss, tectoSeqs[j].printVarnaCommand())

# print number that worked per flanking sequence
numSuccesful = pd.DataFrame(index=['GC', 'CG', 'UA', 'AU'], columns=['GC', 'CG', 'UA', 'AU'])
for row, col in itertools.product(numSuccesful.index, numSuccesful.columns):
    numSuccesful.loc[row, col] = np.sum([tectoSeq.isSecondaryStructureCorrect()[0] for tectoSeq in tectoSeqs[row][col]])


numCores = 20
f = functools.partial(findInitialDistance, allSeqs)
workerPool = multiprocessing.Pool(processes=numCores) #create a multiprocessing pool that uses at most the specified number of cores
x = workerPool.map(f, allSeqs.index)
workerPool.close()
workerPool.join()

results = pd.concat([pd.DataFrame(x[i], columns=[i]).transpose() for i in  allSeqs.index])

topology = 'B1'
sns.lmplot('offset' , 'distance', results, hue='junction_seq', fit_reg=False, palette='Set1')
plt.savefig(os.path.join(figDirectory, 'initial_distance.%s.pdf'%topology))
sns.lmplot('offset' , 'distance', results,  palette='Set1', x_estimator=np.mean, fit_reg=False)
plt.savefig(os.path.join(figDirectory, 'initial_distance.%s.mean.pdf'%topology))

# colormap

