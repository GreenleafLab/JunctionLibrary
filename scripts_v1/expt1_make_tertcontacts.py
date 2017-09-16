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
log = logging.getLogger('sciluigi-interface')

class MyWorkflow(sciluigi.WorkflowTask):
    # only required parameter
    outdir = luigi.Parameter()
    seq_params_dir = os.path.expanduser('~/JunctionLibrary/seq_params/')
    # tl/tlr subsets
    tlr_1_filename = luigi.Parameter(default=seq_params_dir+'receptors_all.dat')
    tlr_2_filename = luigi.Parameter(default=seq_params_dir+'receptors_expt2.dat')
    tlr_3_filename = luigi.Parameter(default=seq_params_dir+'receptors_expt3.dat')
    tlr_4_filename = luigi.Parameter(default=seq_params_dir+'receptors_expt4.dat')
    tlr_5_filename = luigi.Parameter(default=seq_params_dir+'receptors_2.dat')
    tlr_6_filename = luigi.Parameter(default=seq_params_dir+'receptors_evol_paths.dat')
    
    # helix_seq subsets
    helix_50_filename = luigi.Parameter(default=seq_params_dir+'helix_sides_50.dat')
    helix_30_filename = luigi.Parameter(default=seq_params_dir+'helix_sides_30.dat')
    helix_8_filename = luigi.Parameter(default=seq_params_dir+'helix_sides_8.dat')
    helix_5_filename = luigi.Parameter(default=seq_params_dir+'helix_sides_5.dat')
    
    # base subsets
    base_1_filename = luigi.Parameter(default=seq_params_dir+'base_1.dat')
    base_6expt_filename = luigi.Parameter(default=seq_params_dir+'bases_6.dat')
    
    # others
    adapters_filename = luigi.Parameter(default=seq_params_dir+'adapters_1.dat')
    loops_filename = luigi.Parameter(default=seq_params_dir+'loop_1.dat')
    
    def workflow(self):
        
        ## convert some of the input parameters to task outputs
        make_receptors_all = self.new_task('make_receptors_all', processing.ParamToOut, filename=self.tlr_1_filename)
        make_receptors_e2 = self.new_task('make_receptors_e2', processing.ParamToOut, filename=self.tlr_2_filename)
        make_receptors_e3 = self.new_task('make_receptors_e3', processing.ParamToOut, filename=self.tlr_3_filename)
        make_receptors_e4 = self.new_task('make_receptors_e4', processing.ParamToOut, filename=self.tlr_4_filename)
        make_receptors_e5 = self.new_task('make_receptors_e5', processing.ParamToOut, filename=self.tlr_5_filename)
        make_receptors_evol = self.new_task('make_receptors_evol', processing.ParamToOut, filename=self.tlr_6_filename)
        make_helix50 = self.new_task('make_helix_50', processing.ParamToOut, filename=self.helix_50_filename)
        make_helix30 = self.new_task('make_helix_30', processing.ParamToOut, filename=self.helix_30_filename)
        make_helix5 = self.new_task('make_helix_5', processing.ParamToOut, filename=self.helix_5_filename)
        make_base1 = self.new_task('make_base1', processing.ParamToOut, filename=self.base_1_filename)
        make_base6 = self.new_task('make_base6', processing.ParamToOut, filename=self.base_6expt_filename)
        
        ##### SUBLIBRARY 1 ##### 
        # all receptors with 50 helical contexts
        make_sublib1 = self.new_task('make_sublib1', processing.MakeLibrary, lib_name='sublib1', outdir=self.outdir,
                                     loops=self.loops_filename, adapters=self.adapters_filename)
        make_sublib1.in_receptors = make_receptors_all.out_file
        make_sublib1.in_helices = make_helix50.out_file
        make_sublib1.in_bases = make_base1.out_file
        
        ##### SUBLIBRARY 2 ##### 
        # single muts of a couple positions of the receptor
        mutate_receptors_e2 = self.new_task('mutate_receptors_e2', processing.MakeMuts, outdir=self.outdir, positions='1,2;2,3,4')
        mutate_receptors_e2.in_filename = make_receptors_e2.out_file
        
        # assemble
        make_sublib2 = self.new_task('make_sublib2', processing.MakeLibrary, lib_name='sublib2', outdir=self.outdir,
                                     loops=self.loops_filename, adapters=self.adapters_filename)
        make_sublib2.in_receptors = mutate_receptors_e2.out_file
        make_sublib2.in_helices   = make_helix5.out_file
        make_sublib2.in_bases     = make_base1.out_file        
        
        ##### SUBLIBRARY 3 ##### 
        # single muts of fewer receptors and more contexts
        mutate_receptors_e3 = self.new_task('mutate_receptors_e3', processing.MakeMuts, outdir=self.outdir, positions='0,1,2;2,3,4,5')
        mutate_receptors_e3.in_filename = make_receptors_e3.out_file
        
        # assemble
        make_sublib3 = self.new_task('make_sublib3', processing.MakeLibrary, lib_name='sublib3', outdir=self.outdir,
                                     loops=self.loops_filename, adapters=self.adapters_filename)
        make_sublib3.in_receptors = mutate_receptors_e3.out_file
        make_sublib3.in_helices   = make_helix50.out_file
        make_sublib3.in_bases     = make_base1.out_file
        
        ##### SUBLIBRARY 4 #####
        # change the base helix below the receptor
        make_sublib4 = self.new_task('make_sublib4', processing.MakeLibrary, lib_name='sublib4', outdir=self.outdir,
                                     loops=self.loops_filename, adapters=self.adapters_filename)
        make_sublib4.in_receptors = make_receptors_e4.out_file
        make_sublib4.in_helices   = make_helix30.out_file
        make_sublib4.in_bases     = make_base6.out_file
        
        ##### SUBLIBRARY 5 #####
        # make double mutants
        mutate_receptors2 = self.new_task('mutate_receptors2', processing.MakeMuts, outdir=self.outdir,)
        mutate_receptors2.in_filename = make_receptors_e5.out_file
        mutate_receptors_again = self.new_task('mutate_receptors_again', processing.MakeMuts, outdir=self.outdir,)
        mutate_receptors_again.in_filename = mutate_receptors2.out_file
        # assemble
        make_sublib5 = self.new_task('make_sublib5', processing.MakeLibrary, lib_name='sublib5', outdir=self.outdir,
                                     loops=self.loops_filename, adapters=self.adapters_filename)
        make_sublib5.in_receptors = mutate_receptors_again.out_file
        make_sublib5.in_helices   = make_helix5.out_file
        make_sublib5.in_bases     = make_base1.out_file

        ##### SUBLIBRARY 6 #####
        # evol intermediates of certain terts
        make_sublib6 = self.new_task('make_sublib6', processing.MakeLibrary, lib_name='sublib6', outdir=self.outdir,
                                     loops=self.loops_filename, adapters=self.adapters_filename)
        make_sublib6.in_receptors = make_receptors_evol.out_file
        make_sublib6.in_helices   = make_helix30.out_file
        make_sublib6.in_bases     = make_base1.out_file
        
        ##### COMBINE ALL #####
        combine_all = self.new_task('combine_all', CombineAll, outdir=self.outdir)
        combine_all.in_tables = [make_sublib1.out_sequences,
                                 make_sublib2.out_sequences,
                                 make_sublib3.out_sequences,
                                 make_sublib4.out_sequences,
                                 make_sublib5.out_sequences,
                                 make_sublib6.out_sequences]
        return combine_all
        #return make_sublib6


class CombineAll(sciluigi.Task):
    """Dummy class for taking a variable and making it into a ciluigi Target."""
    outdir = luigi.Parameter()
    in_tables = None
    
    # outputs
    def out_table(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'library_all.dat'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            
        # load all targets
        final_table = pd.concat({'tert_contacts_%d'%(i+1):processing.load_file(target().path) for i, target in enumerate(self.in_tables)}, names=['sublibrary'])
        final_table.reset_index(level=0, inplace=True)
        
        final_table.to_csv(self.out_table().path, sep='\t', index=False)
        

if __name__ == '__main__':
    luigi.run(local_scheduler=True, main_task_cls=MyWorkflow)