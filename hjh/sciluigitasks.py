import os
import pandas as pd
import subprocess
import luigi
import sciluigi
import logging
#### SCILUIGI TASKS ####

class ParamToOut(sciluigi.Task):
    """Dummy class for taking a variable and making it into a ciluigi Target."""
    filename = luigi.Parameter()
    
    def out_file(self):
        return sciluigi.TargetInfo(self, self.filename)
    
    def run(self):
        if not os.path.exists(self.filename):
            log.error('File %s does not exist'%self.filename)
            sys.exit()
            


class MakeLibrary(sciluigi.Task):
    """Make a motif using HOMER and consensus seq."""
    lib_name = luigi.Parameter()
    outdir = luigi.Parameter()
    adapters = luigi.Parameter()
    loops = luigi.Parameter()
    #receptors = luigi.Parameter()
    #helices = luigi.Parameter()
    #bases = luigi.Parameter()
    # input
    in_receptors = None
    in_helices = None
    in_bases = None
    
    # outputs
    def out_sequences(self):
        return sciluigi.TargetInfo(self, os.path.join(self.outdir, 'library_%s'%self.lib_name + '.dat'))

    # run
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        # call peak summits with macs2
        # call peak summits with macs2
        tmp_outfile = os.path.splitext(self.out_sequences().path)[0] + '_noadapters.dat'
        task_call1 = ('python -m assemble_tecto_seqs -a {loop} '
        '-b  {helices} {receptor} {base} -out {out}').format(loop=self.loops,
            helices = self.in_helices().path,
            base = self.in_bases().path,
            receptor = self.in_receptors().path,
            out = tmp_outfile)
        log.info(task_call1)
        subprocess.call(task_call1, shell=True)                                                                         

        # run again to add adapters
        task_call2 = ('python -m assemble_tecto_seqs -a {old_out} '
        '-b  {adapters} -out {out}').format(
            adapters=self.adapters,
            old_out=tmp_outfile,
            out=self.out_sequences().path)                                                               
        
        # print and run task
        log.info(task_call2)
        subprocess.call(task_call2, shell=True)
        
class MakeMuts(sciluigi.Task):
    """Dummy class for taking a variable and making it into a ciluigi Target."""
    outdir = luigi.Parameter()
    positions = luigi.Parameter(default=None)
    
    in_filename = None
    
    def out_file(self):
        return sciluigi.TargetInfo(self,  os.path.join(self.outdir, os.path.splitext(os.path.basename(self.in_filename().path))[0] + '_muts.dat'))
    
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        if self.positions:
            add_str = ' -p "%s"'%self.positions
        else:
            add_str = ''
        task_call = ('python ~/JunctionLibrary/mutate_seqs.py -s {in_file} -o {out_file}' + add_str).format(
            in_file=self.in_filename().path, out_file=self.out_file().path, pos=self.positions)
            
        # print and run task
        log.info(task_call)
        subprocess.call(task_call, shell=True)

class MakeHelixSeq(sciluigi.Task):
    """Dummy class for taking a variable and making it into a ciluigi Target."""
    outdir = luigi.Parameter()

    in_junctions = None
    in_helices = None
    in_lenght_positions = None
    
    def out_file(self):
        return sciluigi.TargetInfo(self,  os.path.join(self.outdir, os.path.splitext(os.path.basename(self.in_filename().path))[0] + '_muts.dat'))
    
    def run(self):
        # make out directory if it doesn't exist
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        
        if self.positions:
            add_str = ' -p "%s"'%self.positions
        else:
            add_str = ''
        task_call = ('python ~/JunctionLibrary/mutate_seqs.py -s {in_file} -o {out_file}' + add_str).format(
            in_file=self.in_filename().path, out_file=self.out_file().path, pos=self.positions)
            
        # print and run task
        log.info(task_call)
        subprocess.call(task_call, shell=True)