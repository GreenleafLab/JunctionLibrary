import luigi
import sciluigi
import os
import sys
import subprocess
import logging
import itertools
import pandas as pd
import numpy as np
import ipdb
from hjh import processing
from rnamake import sqlite_library
log = logging.getLogger('sciluigi-interface')

class MyWorkflow(sciluigi.WorkflowTask):
    # only required parameter
    outdir = luigi.Parameter()
    seq_params_dir = os.path.expanduser('~/JunctionLibrary/seq_params/')
    # tl/tlr subsets
    receptors_filename = luigi.Parameter(default=seq_params_dir+'receptors_11ntR.dat')

    # 3way junctions
    threeway_junctions_natural = luigi.Parameter(default=seq_params_dir+'three_way_junctions.dat')
    #threeway_junctions_AA = luigi.Parameter(default=seq_params_dir+'')
    
    # helix_seq subsets
    helix_contexts_filename = luigi.Parameter(default=seq_params_dir+'threeway_helices_10.dat')

    
    # base subsets
    base_1_filename = luigi.Parameter(default=seq_params_dir+'base_1.dat')
    
    # others
    adapters_filename = luigi.Parameter(default=seq_params_dir+'adapters_1.dat')
    #loops_filename = luigi.Parameter(default=seq_params_dir+'loop_1.dat')
    
    def workflow(self):
        
        ## convert some of the input parameters to task outputs
        make_receptors_all = self.new_task('make_receptors_all', ParamToOut, filename=self.tlr_1_filename)

        ##### SUBLIBRARY 1 ##### 
        # all natural junctions
        getjunctions = self.new_task('get_junctions', GetNaturalJunctions, outdir = self.outdir)
        
        # also get control seqs
        getjunctionscont = self.new_task('getjunctionscont', Get3wayJunctionControls, outdir = self.outdir)
        getjunctionscont.in_sequences = getjunctions.out_sequences
        
        # combine each to make librar
        threewaytoseqs = self.new_task('threewaytoseqs', ThreeWayToSeqs, outdir=self.outdir, )
        threewaytoseqs.in_sequences = getjunctions.out_sequences
        
        make_sublib1 = self.new_task('make_sublib1', MakeLibrary, lib_name='sublib1', outdir=self.outdir,
                                     loops=self.loops_filename, adapters=self.adapters_filename)
        make_sublib1.in_receptors = make_receptors_all.out_file
        make_sublib1.in_helices = make_helix50.out_file
        make_sublib1.in_bases = make_base1.out_file
        
        ##### SUBLIBRARY 2 ##### 

        return make_sublib1

class GetNaturalJunctions(sciluigi.Task):
    """Make a motif using HOMER and consensus seq."""
    outdir = luigi.Parameter()
    # input

    # outputs
    def out_sequences(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'threeway_junctions_natural.dat'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        # load
        mlib = sqlite_library.MotifSqliteLibrary("nway")
        mlib.load_all()
        motifs = [m for m in mlib.all() if len(m.ends)==3]
        
        # find motif sequences
        junction_seqs = {}
        for i, m in enumerate(motifs):
            seq_list = m.secondary_structure.sequence().split('&')
            if len(seq_list) == 3:
                junction_seqs[i] = pd.Series(seq_list + [m.name], index=['j1', 'j2', 'j3', 'tw_name'])
            else:
                print m.name
        junction_seqs = pd.concat(junction_seqs).unstack()
        
        # reduce duplicates
        junction_seqs_red = {}
        for name, group in junction_seqs.groupby('tw_name'):
            # take the one with the longest j1 sequence
            idx = group.j1.str.len().sort_values(ascending=False).index[0]
            junction_seqs_red[idx] = group.loc[idx]
        junction_seqs_red = pd.concat(junction_seqs_red).unstack()
        
        # make permutations
        junction_seqs_perm = {}
        for idx, row in junction_seqs_red.iterrows():
            for perm in [0,1,2]:
                seqs = processing.convert_junction_to_seq(row, perm)
            
                junction_seqs_perm[(idx, perm)] = pd.concat([
                    pd.Series(seqs, index=processing.return_junction_fields()),
                    pd.Series({'perm':perm}),
                    row.drop(processing.return_junction_fields())])
        junction_seqs_perm = pd.concat(junction_seqs_perm).unstack()
        
        junction_seqs_perm.to_csv(self.out_sequences().path, sep='\t', index=False)

class Get3wayJunctionControls(sciluigi.Task):
    """Make a motif using HOMER and consensus seq."""
    outdir = luigi.Parameter()
    # input
    in_sequences = None
    
    # outputs
    def out_sequences(self):
        return sciluigi.TargetInfo(self,
                                   os.path.join(self.outdir, get_basename(self.in_sequences().path) + '_controls.dat'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # load sequences
        junction_seqs = processing.load_file(self.in_sequences().path)
        
        # cut out the 'extra'
        junction_seq_controls = {}
        for idx, row in junction_seqs.iterrows():
            for loop_context in ['L1','L2']:
                if loop_context == 'L1':
                    two_way_seq = pd.Series({'side1':row.j1,
                                             'side2':row.j2[0] + row.j3[-1]})
                elif loop_context == 'L2':
                    two_way_seq = pd.Series({'side1':row.j1[0] + row.j2[-1],
                                             'side2':row.j3})
                junction_seq_controls[(idx, loop_context)] = pd.concat([
                    two_way_seq, pd.Series({'loop_context':loop_context})])             
        junction_seq_controls = pd.concat(junction_seq_controls).unstack()    
        junction_seq_controls.to_csv(self.out_sequences().path, sep='\t', index=False)

class Make3Ways(sciluigi.Task):
    """Make a motif using HOMER and consensus seq."""
    outdir = luigi.Parameter()
    receptor_filename = luigi.Parameter()
    base_filename = luigi.Parameter()
    adapter_filename = luigi.Parameter()
    helices_filename = luigi.Parameter()
    
    #inputs
    in_sequences = None

    # outputs
    def out_sequences(self):
        return sciluigi.TargetInfo(self,
                                   os.path.join(self.outdir, 'designed_threeways.dat'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        task_call1 = "python -m assemble_threeways -c ~/JunctionLibrary/seq_params/three_way_helices.dat -r designed_threeways.dat -out designed_threeways_assembled.dat"
        task_call2 = "python -m assemble_tecto_seqs -a designed_threeways_assembled.dat -b ~/JunctionLibrary/seq_params/receptors_11nt.dat ~/JunctionLibrary/seq_params/base_1.dat -out designed_threeways_assembled_receptor.dat"
        
        

class ThreeWayToSeqs(sciluigi.Task):
    """Make a motif using HOMER and consensus seq."""
    outdir = luigi.Parameter()

    # input
    in_sequences = None
    in_helices = None
    
    # outputs
    def out_sequences(self):
        return sciluigi.TargetInfo(self,
                                   os.path.join(self.outdir, get_basename(self.in_sequences().path) + '_assembled.dat'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        task_call = 'python -m assemble_threeways -c %s -r %s -out %s'%(
            self.in_helices().path, self.in_sequences().path, self.out_sequences().path)
        

def get_basename(filename):
    return os.path.basename(os.path.splitext(filename)[0])

  
if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow)