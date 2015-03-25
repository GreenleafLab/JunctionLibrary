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

# load custom libraries
from hjh.helix import Helix
from hjh.junction import Junction


#set up command line argument parser
parser = argparse.ArgumentParser(description="script for making library")
parser.add_argument('-map','--map_file', help='file that contains parameters to vary', required=True)
parser.add_argument('-sub','--subset_number', help='maximum number of sequences per junction motif. Default is all possible ')

parser.add_argument('-par','--seq_params', help='location of seq_param file. Default is ~/JunctionLibrary/seq_params/seq_params.txt')
parser.add_argument('-jun','--junction_sequences', help='overwrite the junction motifs given in map with particular sequenecs in this file')
parser.add_argument('-ter','--tertiary_contacts', help='overwrite the tertiary contacts given in map with particular sequenecs in this file')
parser.add_argument('-seq','--helix_sequences', help='overwrite the helices given in map with particular sequenecs in this file')
parser.add_argument('-out','--out_file', help='file to save output')

if not len(sys.argv) > 1:
    parser.print_help()
    sys.exit()
    
def findSecondaryStructure(x):
    seq = x.loc['tecto_sequence']
    seq, ss, energy = subprocess.check_output("echo %s | RNAfold"%seq, shell=True).split()
    return seq, ss


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

def threadTogether(inside, outside):
    """
    Given a sequence 'inside' and a tuple 'outside' ('side1', 'side2'),
    add  the sequences from 'outside' to either side of 'inside'.
    
    outside(side 1) + inside + outside(side 2)
    """
    
    seq = (outside[0] + inside + outside[1])
    return seq

def findSequence(loop, helix, junction, receptor, base):
    sequence = loop
    sequence = threadTogether(sequence, [helix.split.loc['side1', 'after'], helix.split.loc['side2', 'before']])
    sequence = threadTogether(sequence, junction)
    sequence = threadTogether(sequence, [helix.split.loc['side1', 'before'], helix.split.loc['side2', 'after']])
    sequence = threadTogether(sequence, receptor)
    sequence = threadTogether(sequence, base)
    
    colormap = [4]*len(loop)
    colormap = threadTogether(colormap, [[0.5]*len(helix.split.loc['side1', 'after']), [0.5]*len(helix.split.loc['side2', 'before'])])
    colormap = threadTogether(colormap, [[2]*len(junction.loc['side1']), [2]*len(junction.loc['side2'])])
    colormap = threadTogether(colormap, [[0.5]*len(helix.split.loc['side1', 'before']), [0.5]*len(helix.split.loc['side2', 'after'])])
    colormap = threadTogether(colormap, [[3.5]*len(receptor.loc['side1']), [3.5]*len(receptor.loc['side2'])])
    colormap = threadTogether(colormap, [[0]*len(base.loc['side1']), [0]*len(base.loc['side2'])])
        
    return sequence, colormap
    

#parse command line arguments
args = parser.parse_args()
if args.out_file is None: args.out_file = os.path.splitext(args.map_file)[0] + '.library'
if args.seq_params is None: args.seq_params = os.path.expanduser('~/JunctionLibrary/seq_params/seq_params.txt')

varnaScript = '~/VARNAv3-92.jar'
# load the parameters
seq_params = pd.read_table(args.seq_params, index_col=[0,1])

# read map file
expt_params = pd.read_table(args.map_file)

# if seq file is defined, replace junction motifs with this
if args.junction_sequences is not None:
    expt_params.loc[:, 'junction'] = ['user_defined'] +[np.nan]*(len(expt_params)-1)

# iterate through all
cols = expt_params.columns.tolist()+['junction_seq', 'helix_seq','tecto_sequence', 'sequence']
numToLog =  np.product((~expt_params.isnull()).sum().values)
allSeqs = pd.DataFrame(columns=cols)
logSeqs = pd.DataFrame(index=np.arange(numToLog), columns=expt_params.columns.tolist()+['number'])

for i, x in enumerate(itertools.product(*[expt_params.loc[:,name].dropna() for name in expt_params])):
    
    x = pd.Series(x, index=expt_params.columns.tolist())
    if x.junction == 'user_defined':
        junction = Junction(sequences=pd.read_table(args.junction_sequences, index_col=0).loc[:,['side1', 'side2']])
    else:                
        junction = Junction(tuple(x.junction.split(',')))
        helix = Helix(seq_params.loc[('helix', x.helix)], junction.length   , x.offset, int(x.length))
    
    allSeqSub = pd.DataFrame(index=junction.sequences.index, columns=cols)
    logSeqs.loc[i] = x
    logSeqs.loc[i, 'number'] = len(junction.sequences)

    
    for loc in junction.sequences.index:
        if x.junction == 'user_defined':    # user defined junctions can have any length- so it depends on sequence, not motif
            helix = Helix(seq_params.loc[('helix', x.helix)],
                          junction.findEffectiveJunctionLength(sequence=junction.sequences.loc[loc]),
                          x.offset, int(x.length))        
        sequence, colormap = findSequence(seq_params.loc[('loop', x.loop), 'loop'],
                                          helix,
                                          junction.sequences.loc[loc],
                                          seq_params.loc[('receptor', x.loc['receptor']), ['side1', 'side2']],
                                          seq_params.loc[('base', x.loc['base']), ['side1', 'side2']])
        x.loc['junction_seq'] = '_'.join(junction.sequences.loc[loc])
        x.loc['helix_seq'] = '&'.join(['_'.join(helix.split.loc['side1']), '_'.join(helix.split.loc['side2'])])
        x.loc['tecto_sequence'] = sequence
        x.loc['sequence']  = threadTogether(sequence, seq_params.loc[('adapters', x.adapters), ['side1', 'side2']]).replace('T', 'U')
        allSeqSub.loc[loc] = x
        allSeqSub.loc[loc, 'junction'] = junction.sequences.loc[loc].name
        seq, ss = findSecondaryStructure(x)
        print 'java -cp %s fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "%s" -structureDBN "%s" -colorMap "%s" -colorMapStyle blue -o %s'%(varnaScript, seq, ss, ';'.join(np.array(colormap, dtype=str)), 'test.png' )
    allSeqs = allSeqs.append(allSeqSub, ignore_index=True)
    
sys.exit()

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

