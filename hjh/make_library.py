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
parser.add_argument('-jun','--junction_sequences', help='overwrite the junction motifs given in map with particular sequenecs in this file')

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

# if seq file is defined, replace junction motifs with this
if args.junction_sequences is not None:
    expt_params.loc[:, 'junction'] = ['user_defined'] +[np.nan]*(len(expt_params)-1)

# either load all sequences r

# iterate through all
cols = expt_params.columns.tolist()+['junction_seq', 'junction_seq_noflank', 'helix_seq','tecto_sequence', 'sequence', 'tecto_object']
numToLog =  np.product((~expt_params.isnull()).sum().values)
allSeqs = pd.DataFrame(columns=cols)
logSeqs = pd.DataFrame(index=np.arange(numToLog), columns=expt_params.columns.tolist()+['number'])

# initiate structure to store, per motif, and many contexts, whether it could form

with open(args.out_file + '.pkl', 'wb') as output:
    count = 0
    for i, params in enumerate(itertools.product(*[expt_params.loc[:,name].dropna() for name in expt_params])):
        
        params = pd.Series(params, index=expt_params.columns.tolist())
        # if a particular junction file was specified, replace junctions in map file with those
        if params.junction == 'user_defined':
            sequences = pd.read_table(args.junction_sequences, index_col=0).loc[:,['side1', 'side2']]
            junction = Junction(sequences=sequences)
        
        # if it wasn't specified, we can assume just one helix split per junction and can just do this once, here
        else:                
            junction = Junction(tuple(params.junction.split(',')))
            helix = Helix(seq_params.loc[('helix', params.helix)], junction.length, params.offset, int(params.length))
                
        # if up, junction remains the same. else switch sides
        if params.side == 'up':
            pass
        elif params.side == 'down':
            junction.sequences = pd.DataFrame(junction.sequences.loc[:, ['side2', 'side1']].values,
                                                         columns=['side1', 'side2'],
                                                         index=['%d.s'%name for name in junction.sequences.index])
        
        # initialize saving
        allSeqSub = pd.DataFrame(index=junction.sequences.index, columns=cols)
        logSeqs.loc[i] = params
        logSeqs.loc[i, 'number'] = len(junction.sequences)
        print '%4.1f%% complete'%(100*i/float(numToLog))
        # cycle through junction sequences
        for loc in junction.sequences.index:
            
            # user defined junctions can have any length- so helix splitting depends on sequence, not motif
            if params.junction == 'user_defined':
                helix = Helix(seq_params.loc[('helix', x.helix)],
                              junction.findEffectiveJunctionLength(sequence=junction.sequences.loc[loc]),
                              params.offset, int(params.length))
            
            tectoSeq = TectoSeq(seq_params, params, helix, junction.sequences.loc[loc])
            pickle.dump(tectoSeq, output, pickle.HIGHEST_PROTOCOL)
            allSeqSub.loc[loc] = tectoSeq.params

            allSeqSub.loc[loc, 'tecto_object'] = tectoSeq
            #tectoSeqs[count] = tectoSeq; count+=1
            #tectoSeqs[junction.motif[0]][junction.motif[-1]].append(tectoSeq)
        allSeqs = allSeqs.append(allSeqSub, ignore_index=True)
        
# check if successful secondary structure        
allSeqs.loc[:, 'ss'] = hjh.tecto_assemble.getAllSecondaryStructures(allSeqs.loc[:, 'tecto_sequence'])
allSeqs.loc[:, 'ss_correct'] = False
allSeqs.loc[:, 'junction_SS'] = ''
for loc in allSeqs.index:
    ss_correct, junctionSS, ss = allSeqs.loc[loc, 'tecto_object'].isSecondaryStructureCorrect(ss=allSeqs.loc[loc, 'ss'])
    allSeqs.loc[loc, 'ss_correct'] = ss_correct
    allSeqs.loc[loc, 'junction_SS'] = junctionSS
    
allSeqs.drop(['tecto_object'], axis=1).to_csv(args.out_file+'.txt', sep='\t')

#check if any module was successful
iterables = [['GC', 'CG', 'UA', 'AU'], ['GC', 'CG', 'UA', 'AU'], expt_params.length.dropna(), expt_params.offset.dropna().values, expt_params.side.dropna().values]
index = pd.MultiIndex.from_product(iterables, names=['flank1', 'flank2', 'length', 'offset', 'side'])
allsuccess = pd.DataFrame(index = np.unique(allSeqs.loc[allSeqs.side=='up', 'junction_seq_noflank']), columns=index)

for loc in allSeqs.index:
    index = allSeqs.loc[loc, 'junction_seq_noflank']
    index = allSeqs.loc[loc, 'junction_seq_noflank']
    if allSeqs.loc[loc, 'side'] == 'down':
        index = index[::-1]
    columns = (allSeqs.loc[loc, 'junction'].split(',')[0],
               allSeqs.loc[loc, 'junction'].split(',')[-1],
               allSeqs.loc[loc, 'length'],
               allSeqs.loc[loc, 'offset'],
               allSeqs.loc[loc, 'side'])
    allsuccess.loc[index, columns] = allSeqs.loc[loc, 'ss_correct']
    
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

