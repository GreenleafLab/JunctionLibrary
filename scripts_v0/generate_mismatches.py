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
import hjh.mutations

# find all the paths between base pair states
junctionSeqStarts = Junction(tuple(['W','W'])).sequences
junctionSeqStarts = pd.DataFrame(data=junctionSeqStarts.values, index=junctionSeqStarts.loc[:, 'side1'], columns=['side1', 'side2'] )



allNumPaths = pd.DataFrame(data=0, index=junctionSeqStarts.loc[:, 'side1'],columns=junctionSeqEnds.loc[:, 'side1'])
allPaths = {}
for i, startind in enumerate(junctionSeqStarts.index):
    allPaths[startind] = {}
    for j, endind in enumerate(junctionSeqStarts.index):

        if (junctionSeqStarts.loc[startind] == junctionSeqStarts.loc[endind]).all():
            pass
        else:
            paths = hjh.mutations.getPaths(junctionSeqStarts.loc[startind], junctionSeqStarts.loc[endind])
            numPaths = len(paths.index.levels[0])
            allNumPaths.loc[startind,endind] = numPaths
            if numPaths == 24:
                allPaths[startind][endind] = paths

    if allPaths[startind]:
        allPaths[startind] = pd.concat(allPaths[startind])
    else:
        allPaths.pop(startind)
allPaths = pd.concat(allPaths)


# see how many unique sequences per starting point

for startind in allPaths.index.levels[0]:
    seqs = [''.join(s) for s in allPaths.loc[['AA', 'GG']].values]
    print '%s\t%d'%(startind, len(np.unique(seqs)))
    allPaths.loc[startind]
    
# find paths with GU intermediates
numGUs = pd.DataFrame(index=junctionSeqStarts.index, columns=junctionSeqStarts.index)
for startind in allPaths.index.levels[0]:
    endinds = allPaths.loc[startind].index.levels[0]
    for endind in endinds:
        numGuPath = 0
        for pathind in allPaths.loc[startind].loc[endind].index.levels[0]:
            seqs = [''.join(s) for s in allPaths.loc[startind].loc[endind].loc[pathind].values]
            anyGu = False
            for seq in seqs:
                if seq[0] == 'G' and seq[-1] == 'U': anyGu=True
                if seq[-1] == 'G' and seq[0] == 'U': anyGu=True
                if seq[1] == 'G' and seq[2] == 'U': anyGu=True
                if seq[2] == 'G' and seq[1] == 'U': anyGu=True
            numGuPath += int(anyGu)
        numGUs.loc[startind, endind] = numGuPath

plt.figure(figsize=(5,4))
ax = sns.heatmap(allNumPaths.astype(int), cbar=True, square=True, annot=False, fmt='d')
ax.set_yticklabels(junctionSeqStarts.index[::-1], rotation=0, fontsize=10)
ax.set_xticklabels(junctionSeqStarts.index, rotation=90, fontsize=10)
plt.savefig('150311_library_v2/junction_conformations/length_of_path.MM.pdf')

plt.figure(figsize=(5,4))
ax = sns.heatmap(numGUs.astype(int), cbar=True, square=True, annot=False, fmt='d')
ax.set_yticklabels(junctionSeqStarts.index[::-1], rotation=0, fontsize=10)
ax.set_xticklabels(junctionSeqStarts.index, rotation=90, fontsize=10)
plt.savefig('150311_library_v2/junction_conformations/numberGu.MM.pdf')


for startind in allPaths.index.levels[0]:
    seqs = [''.join(s) for s in allPaths.loc[['AA', 'GG']].values]
    print '%s\t%d'%(startind, len(np.unique(seqs)))
    allPaths.loc[startind]