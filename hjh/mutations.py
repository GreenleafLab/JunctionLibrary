#!/usr/bin/env python

# Author: Sarah Denny, Stanford University 

# Provides modular tools for making a helix-junction
# helix DNA library

##### IMPORT #####
import numpy as np
import pandas as pd
import nwalign as nw

def make_mutation(base):
    if base == 'T': base = 'U'
    bases = ['G', 'C', 'U', 'A']
    bases.remove(base)
    return bases

def complement(base):
    seqDict = {'A':'U', 'G':'C', 'C':'G', 'U':'A'}
    return seqDict[base]
    
def singlesDF(sequenceDF):  
    combinatorics1 = pd.DataFrame(columns=sides)
    combinatorics1.loc[:, 'side1'] = singles(sequenceDF.loc['side1'])
    combinatorics1.loc[:, 'side2'] = sequenceDF.loc['side2']
    
    combinatorics2 = pd.DataFrame(columns=sides)
    combinatorics2.loc[:, 'side2'] = singles(sequenceDF.loc['side2'])
    combinatorics2.loc[:, 'side1'] = sequenceDF.loc['side1']
    
    return pd.concat([combinatorics1, combinatorics2], ignore_index=True).drop_duplicates()

def singles(sequence):
    total_number_bases = len(sequence)
    new_sequences = [sequence]
    for i, base in enumerate(sequence):
        other_bases = make_mutation(base)
        for new_base in other_bases:
            newseq = list(sequence)
            del newseq[i]
            newseq.insert(i,new_base )
            new_sequences += [''.join(newseq)]
    return new_sequences

def singleStep(vecStart, vecEnd):
    vecSteps = []
    for i, base in enumerate(vecEnd):
        vecStep = ''.join(list(vecStart)[:i] + [base] + list(vecStart)[(i+1):])
        if vecStep != vecStart:
            vecSteps.append(vecStep)
        if vecStep == vecEnd:
            pathFinished = True
        else:
            pathFinished = False
    return vecSteps
       
def getPathsFromFormatted(junctionSeqStart, junctionSeqEnd):
    
    # if length of either side is different, align them
    for side in ['side1', 'side2']:
        if len(junctionSeqStart.loc[side]) != len(junctionSeqEnd.loc[side]):
            a, b = nw.global_align(junctionSeqStart.loc[side], junctionSeqEnd.loc[side], gap_open=-10)
            print side, a, b
            junctionSeqStart.loc[side] = a.replace('-', '_')
            junctionSeqEnd.loc[side] = b.replace('-', '_')
         
    
    vecInit = junctionSeqStart.loc['side1'] + junctionSeqStart.loc['side2']
    vecEnd = junctionSeqEnd.loc['side1'] + junctionSeqEnd.loc['side2']
    numSteps = max(len(vecInit), len(vecEnd))
    
    # add bases to vecInit
    while len(vecInit) < len(vecEnd):
        vecInit += '_'
    
    while len(vecInit) > len(vecEnd):
        vecEnd += '_'
        
    steps = getSteps(vecInit, vecEnd, numSteps)
    graph = getGraphFromSteps(steps)
    paths = np.array(list(dfs_paths(graph, vecInit, vecEnd)))
    paths_parsed = formatPaths(paths, numSteps/2)
    return steps, graph, paths, paths_parsed

def getSteps(vecInit, vecEnd, numSteps):
    
    # make sure these are different sequences
    if vecInit == vecEnd:
        print 'start is the same as end'
        return
    
    vecStepsAll = {}
    vecStepsAll[0] = {vecInit:singleStep(vecInit, vecEnd)}
    
    for i in range(1,numSteps):
        vecStepsAll[i] = {}
        for vecStart in np.ravel(vecStepsAll[i-1].values()):
            vecStepsAll[i][vecStart] = singleStep(vecStart, vecEnd)
    return vecStepsAll

def getGraphFromSteps(steps):
    # to make a graph, you have to flatten it
    flattened = {}
    for step in steps.keys():
        for key, values in steps[step].items():
            flattened[key] = set(values)
    return flattened

def formatPaths(paths, numBases):
    # find paths
    iterables = [np.arange(paths.shape[0]), np.arange(paths.shape[1])]
    index = pd.MultiIndex.from_product(iterables, names=['path', 'step'])
    paths_parsed = pd.DataFrame(index=index, columns=['side1', 'side2'])
    for i, path in enumerate(paths):
        for j, seq in enumerate(path):
            paths_parsed.loc[i].loc[j, 'side1'] = seq[:numBases]
            paths_parsed.loc[i].loc[j, 'side2'] = seq[-numBases:]
    return paths_parsed

def dfs_paths(graph, start, goal):
    stack = [(start, [start])]
    while stack:
        (vertex, path) = stack.pop()
        for next in graph[vertex] - set(path):
            if next == goal:
                yield path + [next]
            else:
                stack.append((next, path + [next]))