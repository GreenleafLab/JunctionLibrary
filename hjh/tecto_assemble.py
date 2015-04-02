import pandas as pd
import os
import subprocess
import numpy as np

from hjh.helix import Helix

class TectoSeq():
    def __init__(self, seq_params,  params=None, helix=None, junction=None, sequence=None, reverse_junction=None):

        self.loop     = seq_params.loc[('loop', params.loop), 'loop']
        self.receptor = seq_params.loc[('receptor', params.loc['receptor']),  ['side1', 'side2']]
        self.helix    = helix
        self.base     = seq_params.loc[('base', params.loc['base']),  ['side1', 'side2']]
        self.adapters = seq_params.loc[('adapters', params.loc['adapters']),  ['side1', 'side2']]
        self.sequence = sequence
        
        # decide whether to flip junction or not
        if reverse_junction is None: reverse_junction = False
        if reverse_junction:
            junction.index = junction.index[::-1]
        self.junction = junction
        
        # save params in a series
        self.params = params.append(self.returnInfo())
        self.name = '_'.join(params.loc[['junction', 'length', 'offset']].astype(str))
        
        # if junction has a name, change params name
        if junction is not None:
            if 'name' in junction.index:
                self.params.loc['junction'] = junction.loc['name']
        
    def threadTogether(self, inside, outside):
        return (outside[0] + inside + outside[1])
    
    def makeTecto(self):    
        seq = self.loop
        seq = self.threadTogether(seq, [self.helix.split.loc['side1', 'after'], self.helix.split.loc['side2', 'before']])
        seq = self.threadTogether(seq, [self.junction.loc['side1'], self.junction.loc['side2']])
        seq = self.threadTogether(seq, [self.helix.split.loc['side1', 'before'], self.helix.split.loc['side2', 'after']])
        seq = self.threadTogether(seq, self.receptor)
        seq = self.threadTogether(seq, self.base)
        return seq
    
    def makeTectoWAdapters(self):
        self.threadTogether(self.makeTecto(), self.adapters)
        
    def makeColormap(self, tectoSeq=None):
        
        # if tectoSeq is defined, make everything the same color
        if tectoSeq is not None:
            colormap = np.ones(len(tectoSeq))
        else:
            regions = self.getRegions()
            colormap = regions.copy()
            colormap.loc[regions == 'loop'] = 0.5
            colormap.loc[regions == 'helix'] = 0
            colormap.loc[regions == 'junction'] = 2
            colormap.loc[regions == 'receptor'] = 3.5
            colormap.loc[regions == 'base'] = 0
        return colormap
    
    def getRegions(self):
        regions = ['loop']*len(self.loop)
        regions = self.threadTogether(regions, [['helix']*len(self.helix.split.loc['side1', 'after']), ['helix']*len(self.helix.split.loc['side2', 'before'])])
        regions = self.threadTogether(regions, [['junction']*len(self.junction.loc['side1']), ['junction']*len(self.junction.loc['side2'])])
        regions = self.threadTogether(regions, [['helix']*len(self.helix.split.loc['side1', 'before']), ['helix']*len(self.helix.split.loc['side2', 'after'])])
        regions = self.threadTogether(regions, [['receptor']*len(self.receptor.loc['side1']), ['receptor']*len(self.receptor.loc['side2'])])
        regions = self.threadTogether(regions, [['base']*len(self.base.loc['side1']), ['base']*len(self.base.loc['side2'])])     
        return pd.Series(regions)
    
    def findSecondaryStructure(self, tectoSeq=None):
        if tectoSeq is None:
            tectoSeq = self.makeTecto()
        seq, ss, energy = subprocess.check_output("echo %s | RNAfold"%tectoSeq, shell=True).split()
        return seq, ss
    
    def printVarnaCommand(self, tectoSeq=None, name=None):
        if name is None: name = 'test.png'
        varnaScript = '~/VARNAv3-92.jar'
        seq, ss = self.findSecondaryStructure(tectoSeq)
        colormap = self.makeColormap(tectoSeq)
        commandString = 'java -cp %s fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN "%s" -structureDBN "%s" -colorMap "%s" -colorMapStyle blue -spaceBetweenBases "0.75" -o %s'%(varnaScript, seq, ss, ';'.join(colormap.astype(str)), name)
        return commandString, seq, ss
    
    def returnInfo(self):
        
        params = pd.Series(index = ['junction_seq', 'n_flank', 'no_flank', 'flank', 'helix_one_length', 'helix_seq', 'tecto_sequence', 'sequence'] )
        try:
            
            params.loc['junction_seq'] = '_'.join(self.junction.loc[['side1', 'side2']])
            params.loc['helix_seq'] = '&'.join(['_'.join(self.helix.split.loc['side1']), '_'.join(self.helix.split.loc['side2'])])
            params.loc['tecto_sequence'] = self.makeTecto()
            params.loc['sequence']  = self.threadTogether(params.loc['tecto_sequence'], self.adapters.loc[['side1', 'side2']]).replace('U', 'T')
            params.loc['helix_one_length'] = len(self.helix.split.loc['side1', 'before'])
            # things that require junction had other settings (i.e. name and n_flank)
            params.loc['n_flank'] = int(self.junction.loc['n_flank'])
            params.loc['no_flank'] = '_'.join([self.junction.loc[side][params.loc['n_flank']:-params.loc['n_flank']] for side in ['side1', 'side2']])
            params.loc['flank'] = ''.join([self.junction.loc['side1'][:params.loc['n_flank']], self.junction.loc['side1'][-params.loc['n_flank']:]])
            #params.loc['junction'] = self.junction.loc['name']
            
        except AttributeError: pass
        
        return params
    
    def isSecondaryStructureCorrect(self, ss=None, numberFlankingBases=None):
        regions = self.getRegions()
        if numberFlankingBases is None:
            numberFlankingBases = 1
        if ss is None:
            seq, ss = self.findSecondaryStructure()
        ss = pd.Series(list(ss), dtype=str)
        
        expectedSecondaryStructure = {'base':'(((((())))))',
                                      'receptor': '(..(())...)',
                                      'loop': '....',
                                      }
        helixSS = ''.join(ss.loc[regions == 'helix'])
        helixLength = len(helixSS)/2
        expectedSecondaryStructure['helix'] = ''.join(['(']*helixLength + [')']*helixLength)
        
        isCorrect = True
        # base should be all base paired
        if ''.join(ss.loc[regions == 'base']) == expectedSecondaryStructure['base']:
            pass
        else:
            isCorrect = False
            print '\tissue with base region. Expected/actual: %s\t%s'%(expectedSecondaryStructure['base'], ''.join(ss.loc[regions == 'base']) )
        
        # receptor: assume tetraloop receptor
        if ''.join(ss.loc[regions == 'receptor']) == expectedSecondaryStructure['receptor']:
            pass
        else:
            isCorrect = False
            print '\tissue with receptor region. Is it 11nt receptor? Expected/actual: %s\t%s'%(expectedSecondaryStructure['receptor'], ''.join(ss.loc[regions == 'receptor']))
            
        # receptor: assume tetraloop receptor
        if helixSS == expectedSecondaryStructure['helix']:
            pass
        else:
            isCorrect = False
            print '\tissue with helix region. Expected:/actual %s\t%s'%(expectedSecondaryStructure['helix'], helixSS)
        
        # print junction secondary structure
        junctionSS = ''.join(ss.loc[regions == 'junction'])
        if (junctionSS[:numberFlankingBases]  == ''.join(['(']*numberFlankingBases) and
            junctionSS[-numberFlankingBases:] == ''.join([')']*numberFlankingBases) and
            junctionSS.find(''.join(['(']*numberFlankingBases) + ''.join([')']*numberFlankingBases)) != -1):
            pass
        else:
            isCorrect = False
            print '\tissue with junction region. Expected:/actual %s*()*%s\t%s'%(''.join(['(']*numberFlankingBases), ''.join([')']*numberFlankingBases), junctionSS)
            
        return isCorrect, junctionSS, ''.join(ss)

def findTecto(params, junction, seq_params, locs):
    cols = TectoSeq(seq_params, params).params.index
    allSeqSub = pd.DataFrame(index=locs, columns=cols)
    for loc in locs:
        helix = Helix(seq_params.loc[('helix', params.helix)],
                      junction.findEffectiveJunctionLength(sequence=junction.sequences.loc[loc]),
                      params.offset, int(params.length))
                
        tectoSeq = TectoSeq(seq_params, params, helix, junction.sequences.loc[loc])
        allSeqSub.loc[loc] = tectoSeq.params
        allSeqSub.loc[loc, 'tecto_object'] = tectoSeq
    return allSeqSub

def makeFasta(tectoSeqs, outputFile):
    f = open(outputFile, 'w')
    for line in tectoSeqs.index:
        f.write('>%d\n'%line)
        f.write('%s\n'%tectoSeqs.loc[line])
    f.close()
    return

def getAllSecondaryStructures(tectoSeqs):
    outputFile = 'test.fasta'
    makeFasta(tectoSeqs, outputFile)
    outputString = subprocess.check_output('cat %s | RNAfold --noPS'%outputFile, shell=True).lstrip('>').split('>')
    secondaryStructures = ['']*len(outputString)
    indices = np.zeros(len(outputString))
    for i, line in enumerate(outputString):
        indices[i] = line.split('\n')[0]
        secondaryStructures[i] = line.split()[2]
    return pd.Series(secondaryStructures, index=indices )

def getSecondaryStructureMultiprocess(allSeqSub):

    allSeqSub.loc[:, 'ss_correct'] = False
    allSeqSub.loc[:, 'junction_SS'] = ''
    for loc in allSeqSub.index:
        ss_correct, junctionSS, ss = allSeqSub.loc[loc, 'tecto_object'].isSecondaryStructureCorrect(ss=allSeqSub.loc[loc, 'ss'], numberFlankingBases=allSeqSub.loc[loc, 'n_flank'],)
        allSeqSub.loc[loc, 'ss_correct'] = ss_correct
        allSeqSub.loc[loc, 'junction_SS'] = junctionSS
    
    return allSeqSub