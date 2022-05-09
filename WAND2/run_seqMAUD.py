## modify MAUD file - parameter set
# to refine or not

# nate peterson - UVa
# dec 2020

import os,re,sys,contextlib
import pandas as pd
import numpy as np
from math import isnan
from distutils.util import strtobool
from subprocess import Popen
from subprocess import CREATE_NEW_CONSOLE
from tqdm import tqdm

class DummyFile(object):
    
    file = None
    def __init__(self, file):
        self.file = file
    
    def write(self, x):
        # Avoid print() second call (useless \n)
        if len(x.rstrip()) > 0:
            tqdm.write(x, file=self.file)

@contextlib.contextmanager
def nostdout():
    
    save_stdout = sys.stdout
    sys.stdout = DummyFile(sys.stdout)
    yield
    sys.stdout = save_stdout
    
def tqdm_write( s ):
    
    """
    writes string with tqdm
    
    """
    
    assert isinstance(s,str), 'must print string...'

    with nostdout():
        print(s)

class ParFile(object):

    """
    par-file class
    """

    def __init__( self, file, refine_step=0 ):

        ref_reke = "(\(\d+\.\S+\))" #gets any number inside ()
        quo_reke = "'\s*(.*?)\s*'" #gets anything inside a ''      

        if os.path.exists(file):
            
            self.fname = file
            self.lines = {}
            self.lines_refinable = []
            self.phases = []
            
            tqdm_write('   loading in par file..')
            
            with open(self.fname, encoding="utf8") as f:
                
                self.alllines = f.readlines()
                # set this for now.. can be changed ahead
                phase_name = None                
                
                for li,line in enumerate(self.alllines):
                    
                    # skip blank lines
                    if line != '\n':
                        
                        # initialize attributes
                        self.lines[li] = {}
                        self.lines[li]['line'] = line
                        self.lines[li]['key'] = ''
                                         
                        # store these for use later in try block
                        # where a list of values is
                        if 'loop_' in line: 
                            loop_start = li
                        elif '_pd_phase_name' in line:
                            phase_name = line.split('\'')[1]
                            self.phases.append(phase_name)
                            
                        # should be a normal parameter
                        keys = line.split(' ')
                        if '_' in keys[0]:

                            # has multiple items
                            if len(keys) > 1:
                                
                                self.lines[li]['key'] = keys[0].rstrip()

                                # was refined earlier
                                if re.search(ref_reke,keys[1]) is not None:
                                    self.lines[li]['val'] = keys[1].rstrip()
                                    self.lines_refinable.append(li)
                                    self.lines[li]['refined'] = True   
                                    self.lines[li]['phase']   = phase_name                                    
                                
                                # is a boolean key.. may or may not be a true variable
                                elif keys[1].rstrip().lower() == 'true' or keys[1].rstrip().lower() == 'false':
                                    
                                    self.lines[li]['val']   = keys[1].rstrip()
                                    self.lines[li]['phase'] = phase_name
                                    self.lines[li]['refined'] = 'bool'

                                # this is a line with something in '' 
                                elif re.search(quo_reke,line) is not None:
                                    
                                    match = re.search(quo_reke,line)
                                    self.lines[li]['val'] = line[match.start(0):match.end(0)]
                                    self.lines[li]['phase'] = phase_name
                                    self.lines[li]['refined'] = 'NA'

                                # random space in front of variable
                                elif keys[0] == '':
                                                                
                                    self.lines[li]['key'] = ''
                                    self.lines[li]['val'] = None
                                    self.lines[li]['refined'] = 'NA'
                                    self.lines[li]['phase']   = phase_name

                                elif line.startswith('#'):
                                    
                                    self.lines[li]['key'] = '' 
                                    self.lines[li]['val'] = None
                                    self.lines[li]['refined'] = 'NA'
                                    self.lines[li]['phase']   = phase_name                                    

                                else:
                                    # this is variable
                                    try: # is a number
                                        float(keys[1])
                                        self.lines[li]['val'] = keys[1].rstrip()
                                        self.lines[li]['refined'] = False
                                        self.lines[li]['phase']   = phase_name
                                    except: # is not a number
                                        self.lines[li]['val'] = keys[1].rstrip()
                                        self.lines[li]['refined'] = 'NA'
                                        self.lines[li]['phase']   = phase_name
                                        
                            # single item on line
                            elif len(keys) == 1:
                            
                                self.lines[li]['key'] = '' 
                                self.lines[li]['val'] = None
                                self.lines[li]['refined'] = 'NA'
                                self.lines[li]['phase']   = phase_name
                            
                        # should be values under a 'loop_'
                        elif keys[0] == '' and keys[1] == '':
                            
                            # get a loop counter
                            count   = li - (loop_start+1)
                            key_txt = self.alllines[loop_start+1].rstrip()
                            
                            self.lines[li]['key'] = '{}_{:d}'.format(key_txt,
                                                                     count)
                            self.lines[li]['val'] = keys[2].rstrip() 
                            self.lines_refinable.append(li)                                
                            self.lines[li]['phase']   = phase_name                                
                            
                            if re.search(ref_reke,keys[2]) is not None:
                                # need to take line before loop
                                self.lines[li]['refined'] = True
                            else:
                                self.lines[li]['refined'] = False

                        else:
                            try: # see if it is intensity
                                float(keys[1])
                                self.lines.pop(li)
                            except: # nope, skip
                                pass
                            
    
            self.df = pd.DataFrame.from_dict(self.lines,orient='index')   
 
            # go back and label phases
            for pi, phase in enumerate(self.phases):
                
                # we are indexing starting at 1
                key = '_pd_phase_atom_%_{:d}'.format(pi+1)
                
                self.df.loc[ self.df['key'] == key, 'phase' ] = phase
            
            self._refine_step = refine_step
            self._refineable_params = []

        else: raise FileNotFoundError('{}'.format(file))

    def _checkPar( self, key ):
        
        """
        
        search par file for keyword
        
        """
        
        pass

    def _refinePar( self, key, val, debug, search='exact' ):
        
        """

        set a parameter to refine

        """

        assert search in ['exact', 'fuzzy'], 'search (type) must be \'exact\' or \'fuzzy\'' 
        assert val.lower() in ['y', 'n'], 'val must by \'y\' or \'n\''

        # raise NotImplementedError('Still a work in progress..')

        if search == 'exact':
            param_index = self.df.loc[self.df['key'] == key ].index
        elif search == 'fuzzy':
            param_index = self.df.loc[self.df['key'].str.contains(key)].index

        if param_index.empty:
            raise KeyError('key supplied to toggle not present')

        else:
            
            tqdm_write('     For key: \'{}\' - changing to \'{}\''.format(key,val))
            
            for pi in param_index:
                
                if debug:
                    suffix = '@ line# : {}'.format(pi)
                else:
                    suffix = ''
                
                if self.df.loc[pi]['val'] is None: 
                    
                    raise KeyError('        key supplied \'{}\' doesn\'t have a value.. should it be refined?'.format(key))
                
                else:
                    
                    if self.df.loc[pi]['refined'] == 'NA':
                        
                        tqdm_write('        parameter is not refinable... will swap instead')
                        self._setPar( key, val, debug=debug, search=search )
                        
                    elif val == 'n':
                        
                        # we can skip here..
                        if self.df.loc[pi,'refined'] is False:
                            tqdm_write('        Skipped: {} {}'.format(repr(self.df.loc[pi]['val']),suffix))  
                        elif self.df.loc[pi]['refined'] == 'bool' and self.df.loc[pi]['val'] == 'false':
                            tqdm_write('        Skipped: {} {}'.format(repr(self.df.loc[pi]['val']),suffix))  
                        else:
                            # it was refined previously
                            old_val = self.df.loc[pi]['val']
                            if self.df.loc[pi]['refined'] == 'bool': new_val = 'false'
                            else: new_val = old_val.split('(')[0]
                            
                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = False
                            
                            tqdm_write('        Replaced: {} with {} {}'.format(repr(old_val),repr(new_val),suffix))
                    
                    elif val == 'y':
                        
                        # we can skip here..
                        if self.df.loc[pi,'refined'] is True:
                            tqdm_write('        Skipped: {} {}'.format(repr(self.df.loc[pi]['val']),suffix))  
                        elif self.df.loc[pi]['refined'] == 'bool' and self.df.loc[pi]['val'] == 'true':
                            tqdm_write('        Skipped: {} {}'.format(repr(self.df.loc[pi]['val']),suffix))  
                        else:
                            # it was not refined previously
                            old_val = self.df.loc[pi]['val']
                            if self.df.loc[pi]['refined'] == 'bool': new_val = 'true'
                            else: new_val = '{}(0.0)'.format(old_val)
                            
                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = True
                        
                            tqdm_write('        Replaced: {} with {} {}'.format(repr(old_val),repr(new_val),suffix))

    def _setPar( self, key, val, debug=False, search='exact' ):

        """
        set a parameter to a given value
        """

        assert search in ['exact', 'fuzzy'], 'search type must be \'exact\' or \'fuzzy\'' 

        # raise NotImplementedError('Still a work in progress..')

        if search == 'exact':
            param_index = self.df.loc[self.df['key'] == key ].index
        elif search == 'fuzzy':
            param_index = self.df.loc[self.df['key'].str.contains(key)].index

        if param_index.empty:
            raise KeyError('key supplied to toggle not present')

        else:
            tqdm_write('For key: \'{}\''.format(key))
            for pi in param_index:
                
                if debug:
                    suffix = '@ line# : {}'.format(pi)
                else:
                    suffix = ''

                if self.df.loc[pi]['val'] is None: 
                    
                    raise KeyError('key supplied \'{}\' doesn\'t have a value.. should it be refined?'.format(key))
                
                else:
                    
                    ## boolean parameter - should have NaN for refined
                    if isnan(self.df.loc[pi]['refined']):

                        if val.lower() in ['true','false']: pass
                        else: assert isinstance(val,bool), 'parameter is bool.. val should be bool'
                        # it was refined previously
                        old_val = self.df.loc[pi]['val']
                        new_val = str(val).lower()
                        
                        old_line = self.df.loc[pi]['line']
                        new_line = old_line.replace(old_val,new_val)
                        
                        self.df.loc[pi]['line'] = new_line
                        self.df.loc[pi]['refined'] = False
                        
                        tqdm_write('   Replaced: {} with {} {}'.format(old_val,new_val,suffix)) 

                    else:
                        
                        assert isinstance(val,str)

                        if self.df.loc[pi]['refined'] is True:
                            
                            # it was refined previously
                            old_val = self.df.loc[pi]['val']
                            # need to split out err
                            old_err = old_val.split('(')[1]
                            new_val = '{}({}'.format(val,old_err)

                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = False
                            
                            tqdm_write('   Replaced: {} with {} {}'.format(repr(old_val),repr(new_val),suffix))
                            
                        elif self.df.loc[pi]['refined'] is False:

                            # it was not refined previously                     
                            old_val = self.df.loc[pi]['val']
                            new_val = val
                            
                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = True
                            
                            tqdm_write('   Replaced: {} with {} {}'.format(repr(old_val),repr(new_val),suffix))

    def _togglePar( self, key, debug=False, directSearch=True, fixType=None ):
        
        """
        
        toggle a specific keyword
        
        """
        
        if directSearch:
            param_index = self.df.loc[self.df['key'] == key ].index
        else:
            param_index = self.df.loc[self.df['key'].str.contains(key)].index

        if param_index.empty:
            raise KeyError('key supplied to toggle not present')

        else:
            
            tqdm_write('For key: \'{}\''.format(key))
            
            for pi in param_index:
                
                if debug:
                    suffix = '@ line# : {}'.format(pi)
                else:
                    suffix = ''
                
                try: 
                    
                    isnan(self.df.loc[pi]['val'])
                    raise KeyError('key supplied \'{}\' doesn\'t have a value.. should it be refined?'.format(key))
                
                except:
                    
                    if isnan(self.df.loc[pi]['refined']):
                        
                        # boolean
                        if self.df.loc[pi]['val'] == 'true' and fixType != True:
                        
                            # it was refined previously
                            old_val = self.df.loc[pi]['val']
                            new_val = 'false'
                            
                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = False
                            
                            tqdm_write('Replaced: {} with {} {}'.format(old_val,new_val,suffix))                            
                    
                        elif self.df.loc[pi]['val'] == 'false' and fixType != False:
                                   
                            # it was not refined previously
                            old_val = self.df.loc[pi]['val']
                            new_val = 'true'
                            
                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = True
                            
                            tqdm_write('Replaced: {} with {} {}'.format(old_val,new_val,suffix))                            
                            
                        else:
                            pass
                    
                    elif self.df.loc[pi]['refined']:
                        
                        if fixType == True:
                            tqdm_write('already refineable.. going to skip {}'.format(suffix))
                            pass
                        else:
                            # it was refined previously
                            old_val = self.df.loc[pi]['val']
                            new_val = old_val.split('(')[0]
                            
                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = False
                            
                            tqdm_write('Replaced: {} with {} {}'.format(repr(old_val),repr(new_val),suffix))
                        
                    elif self.df.loc[pi]['refined'] is False:
                        # it was not refined previously
                        
                        if fixType == False:
                            tqdm_write('already not refinable.. going to skip {}'.format(suffix))
                            pass
                        else:                        
                            old_val = self.df.loc[pi]['val']
                            new_val = '{}(0.0)'.format(old_val)
                            
                            old_line = self.df.loc[pi]['line']
                            new_line = old_line.replace(old_val,new_val)
                            
                            self.df.loc[pi]['line'] = new_line
                            self.df.loc[pi]['refined'] = True
                            
                            tqdm_write('Replaced: {} with {} {}'.format(repr(old_val),repr(new_val),suffix))

    def _writeBat( self, batch_file, analysis_file, iter_num, new_analysis_file, result_file, wizard_index=None):

        """
        write a batch file

        eventually this should only ask for iter_num..
        """

        with open(batch_file,'w') as fout:

            fout.write('loop_\n')
            fout.write('_riet_analysis_file\n')
            fout.write('_riet_analysis_iteration_number\n')
            if wizard_index is not None:
                fout.write('_riet_analysis_wizard_index\n')
            fout.write('_riet_analysis_fileToSave\n')
            fout.write('_riet_append_simple_result_to\n\n')

            if wizard_index is not None:
                fout.write('\'{}\' {} {} \'{}\' \'{}\'\n'.format(analysis_file, iter_num, wizard_index, new_analysis_file, result_file  ))
            else:
                fout.write('\'{}\' {} \'{}\' \'{}\'\n'.format(analysis_file, iter_num, new_analysis_file, result_file ))

    def change( self, paramKey:str, val, paramHint=None, debug=False ):
        
        """

        change parameter
        simple options - need to expand on this
        
        """
        
        VALID_PARAMS = {'bkgd_pol','scale','lat_par','phase_fr','texture'}
        
        if paramKey not in VALID_PARAMS:
            raise ValueError('paramKey must be one of following: {}'.format(VALID_PARAMS))
        else:
            if paramKey == 'bkgd_pol':
                
                key = '_riet_par_background_pol'
                self._setPar(key, val, debug=debug, search='fuzzy')
                
            elif paramKey == 'scale':

                assert paramHint in ['TOF', 'CW'],"Must select the whether TOF or CW"
                
                if paramHint == 'TOF':
                    
                    raise NotImplementedError('need to add key for incident_inten_spec')
                
                elif paramHint == 'CW':
                    
                    key = '_pd_proc_intensity_incident'
                    self._setPar(key, val, debug=debug, search='exact')
            
            elif paramKey == 'phase_fr':
                
                ## not going to be easy..
                if paramHint in self.phases:
                    # we specifying a phase..
                    raise NotImplementedError('need to test this out..')                    
                
                else:
                    
                    assert isinstance(paramHint,int), 'please tell me what phase to change'
                    
                    # loop over phases.. need to keep #1 off
                    for pi,phase in enumerate(self.phases):
                      
                        if pi > 0:
                            key = '_pd_phase_atom_%_{:d}'.format(pi)
                            self._setPar(key, val, debug=debug, search='exact')
                
            elif paramKey == 'lat_par':
                
                if paramHint in self.phases:
                    # we specifying a phase..
                    raise NotImplementedError('need to test this out..')
                    
                else:
                    
                    key = 'cell_length'
                    self._setPar(key, val, debug=debug, search='fuzzy')
      
        
            elif paramKey == 'texture':
                
                if paramHint in self.phases:
                    # we specifying a phase..
                    raise NotImplementedError('need to test this out..')
                    
                else:
                    
                    key = '_rita_odf_refinable'   
                    self._setPar(key, val, debug=debug, search='fuzzy')

    def refine( self, paramKey:str, val, paramHint=None, debug=False ):

        """

        change parameter
        simple options - need to expand on this
        
        """
        
        VALID_PARAMS = {'bkgd_pol','scale','phase_fr','lat_par','biso','microstructure','texture'}

        assert paramKey in VALID_PARAMS, 'paramKey must be one of following: {}'.format(VALID_PARAMS)
        assert val.lower() in ['y', 'n'], 'val must by \'y\' or \'n\''

        if paramKey == 'bkgd_pol':
            
            key = '_riet_par_background_pol'
            self._refinePar(key=key, val=val, debug=debug, search='fuzzy')
            
        elif paramKey == 'scale':

            assert paramHint in ['TOF', 'CW'],"Must select the whether TOF or CW"
            
            if paramHint == 'TOF':
                
                raise NotImplementedError('need to add key for incident_inten_spec')
            
            elif paramHint == 'CW':
                
                key = '_pd_proc_intensity_incident'
                self._refinePar(key=key, val=val, debug=debug, search='exact')

        elif paramKey == 'phase_fr':
            
            if paramHint in self.phases:
                # we specifying a phase..
                raise NotImplementedError('need to test this out..')                    
            
            elif isinstance(paramHint,int) is True:
                ## not tested..
                key = '_pd_phase_atom_%_{:d}'.format(paramHint)
                self._refinePar(key, val, debug=debug, search='exact')                
                
            else:
                # loop over phases.. need to keep #1 off
                for pi,phase in enumerate(self.phases):
                  
                    if pi > 0:
                        key = '_pd_phase_atom_%_{:d}'.format(pi+1)
                        self._refinePar(key, val, debug=debug, search='exact')
                    else:
                        tqdm_write('     Skipping phase#1.. MAUD needs one phase off')
        
        elif paramKey == 'lat_par':
            
            if paramHint in self.phases:
                # we specifying a phase..
                raise NotImplementedError('need to test this out..')
                
            else:
                
                key = 'cell_length'
                self._refinePar(key=key, val=val, debug=debug, search='fuzzy')
    
        elif paramKey == 'texture':
            
            if paramHint in self.phases:
                # we specifying a phase..
                raise NotImplementedError('need to test this out..')
                
            else:
                
                key = '_rita_odf_refinable'   
                self._refinePar(key=key, val=val, debug=debug, search='fuzzy')

    def save( self, newFile ):

        """
        write out the par file to a new one
        or same file
        """            

        with open(newFile,'w') as fout:

            for li,line in enumerate(self.alllines):

                if li not in self.df.index:
                    fout.write(line)
                else:
                    fout.write(self.df.loc[li]['line'])

# location of main file directory
rootDir = '/'

""" 
needs columns: 
[Time, Run #, Chi, Phi, Sample ]  
"""

stdout=None
stderr=None

# output directory
outDir      = os.path.join(rootDir,'maud')
# experiment metadata
expInfoFile = os.path.join(rootDir,'expSummary.csv')
expInfo     = pd.read_csv(expInfoFile)

# setup run sets for fitting
samples    = pd.unique(expInfo['Sample'])
batchFile  = os.path.join(outDir,'maud_batch.ins')

# progress bar
pbar = tqdm(samples, file=sys.stdout)

# need to use tqdm to setup iter, using stdout as the file
for sample in pbar:

    if 'Steel' in sample:

        sampleDir   = os.path.join(outDir,'samples_texred_15deg')

        # make this play nice with unix
        sample_desc = sample.replace(' #','_')     
        # from reduction script
        initParFile = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_out.par'.format(sample_desc))
        
        # ====> First Refinement Step
        # ===> Background & scale

        tqdm_write('\n --- {} --- '.format(sample_desc))

        # load in par file
        initPar     = ParFile(initParFile)
        # set initial refinable parameters
        initPar.refine('bkgd_pol', 'y')
        initPar.refine('phase_fr', 'n', paramHint=1) #this wasn't turned off in the init files
        initPar.refine('phase_fr', 'n')
        initPar.refine('scale', 'y', paramHint='CW')
        initPar.refine('texture', 'n') #shut tex off
        initPar.refine('lat_par', 'n') #shut lat param off

        # generate new par file to be used
        newParFile       = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_Step1.par'.format(sample_desc))
        newParResultFile = os.path.join(sampleDir,sample_desc,'refineResult.txt')
        relParFile       = os.path.relpath(newParFile,outDir)
        relParResult     = os.path.relpath(newParResultFile,outDir)
        initPar.save(newParFile)
        
        # write batch file that MAUD can read in
        # those paths need to be relative to the batchFile
        initPar._writeBat(batchFile, relParFile, 8, 'HB2C_MAUD_{}_Step1.par'.format(sample_desc), 'refineResult.txt')
        
        # call the batch script in new CMD window
        pbar.set_description('{} : Step 1'.format(sample))
        pbar.refresh()
        tqdm_write('\n --- Calling MAUD for Step 1 --- ')
        proc = Popen(['cmd','/C','maud_batch.bat'],cwd=outDir,creationflags=CREATE_NEW_CONSOLE,stdout=stdout,stderr=stderr)
        proc.wait()        

        # ====> 2nd Refinement Step
        # ===> phase fraction

        tqdm_write('\n --- {} --- '.format(sample_desc))

        # load in par file
        firstPar     = ParFile(newParFile)
        # set initial refinable parameters
        firstPar.refine('phase_fr', 'y')

        # generate new par file to be used
        newParFile       = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_Step2.par'.format(sample_desc))
        newParResultFile = os.path.join(sampleDir,sample_desc,'refineResult.txt')
        relParFile       = os.path.relpath(newParFile,outDir)
        relParResult     = os.path.relpath(newParResultFile,outDir)
        firstPar.save(newParFile)
        
        # write batch file that MAUD can read in
        # those paths need to be relative to the batchFile
        firstPar._writeBat(batchFile, relParFile, 5, 'HB2C_MAUD_{}_Step2.par'.format(sample_desc), 'refineResult.txt')

        # call the batch script in new CMD window
        pbar.set_description('{} : Step 1'.format(sample))
        pbar.refresh()
        tqdm_write('\n --- Calling MAUD for Step 2 --- ')
        proc = Popen(['cmd','/C','maud_batch.bat'],cwd=outDir,creationflags=CREATE_NEW_CONSOLE,stdout=stdout,stderr=stderr)
        proc.wait()        

        # ====> 3rd Refinement Step
        # ===> Lattice

        # reload in after running refinement
        secondPar = ParFile(newParFile)
        # set initial refinable parameters
        secondPar.refine('lat_par', 'y') #shut lat param off      

        # generate new par file to be used
        newParFile       = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_Step3.par'.format(sample_desc))
        newParResultFile = os.path.join(sampleDir,sample_desc,'refineResult.txt')
        relParFile       = os.path.relpath(newParFile,outDir)
        relParResult     = os.path.relpath(newParResultFile,outDir)
        secondPar.save(newParFile)

        # write batch file that MAUD can read in
        # those paths need to be relative to the batchFile
        secondPar._writeBat(batchFile, relParFile, 5, 'HB2C_MAUD_{}_Step3.par'.format(sample_desc), 'refineResult.txt')
        
        # call the batch script in new CMD window
        pbar.set_description('{} : Step 3'.format(sample))
        pbar.refresh()
        tqdm_write('\n --- Calling MAUD for Step 3 --- ')
        proc = Popen(['cmd','/C','maud_batch.bat'],cwd=outDir,creationflags=CREATE_NEW_CONSOLE,stdout=stdout,stderr=stderr)
        proc.wait()   

        # ====> 4th Refinement Step
        # ===> Texture

        # reload in after running refinement
        thirdPar = ParFile(newParFile)
        # set initial refinable parameters
        thirdPar.refine('texture', 'y') 

        # generate new par file to be used
        newParFile       = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_Step4.par'.format(sample_desc))
        newParResultFile = os.path.join(sampleDir,sample_desc,'refineResult.txt')
        relParFile       = os.path.relpath(newParFile,outDir)
        relParResult     = os.path.relpath(newParResultFile,outDir)
        thirdPar.save(newParFile)

        # write batch file that MAUD can read in
        # those paths need to be relative to the batchFile
        thirdPar._writeBat(batchFile, relParFile, 5, 'HB2C_MAUD_{}_Step4.par'.format(sample_desc), 'refineResult.txt')
        
        # call the batch script in new CMD window
        pbar.set_description('{} : Step 4'.format(sample))
        pbar.refresh()
        tqdm_write('\n --- Calling MAUD for Step 4 --- ')
        proc = Popen(['cmd','/C','maud_batch.bat'],cwd=outDir,creationflags=CREATE_NEW_CONSOLE,stdout=stdout,stderr=stderr)
        proc.wait()           
        
        # ### ====== Move to other type ======= ###
        # sampleDir   = os.path.join(outDir,'samples_texred_5deg')

        # # from reduction script
        # initParFile = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_out.par'.format(sample_desc))

        # # ====> First Refinement Step
        # # ===> Background & scale

        # tqdm_write(' --- {} --- \n'.format(sample_desc))

        # # load in par file
        # initPar     = ParFile(initParFile)
        # # set initial refinable parameters
        # initPar.refine('bkgd_pol', 'y')  
        # initPar.refine('scale', 'y', paramHint='CW')
        # initPar.refine('texture', 'n') #shut tex off
        # initPar.refine('lat_par', 'n') #shut lat param off

        # # generate new par file to be used
        # newParFile       = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_Step1.par'.format(sample_desc))
        # newParResultFile = os.path.join(sampleDir,sample_desc,'refineResult_Step1.txt')
        # relParFile       = os.path.relpath(newParFile,outDir)
        # relParResult     = os.path.relpath(newParResultFile,outDir)
        # initPar.save(newParFile)
        
        # # write batch file that MAUD can read in
        # # those paths need to be relative to the batchFile
        # initPar._writeBat(batchFile, relParFile, 5, relParFile, relParResult)
        
        # # call the batch script in new CMD window
        # pbar.set_description('{} : Step 1'.format(sample))
        # pbar.refresh()
        # tqdm_write('\n --- Calling MAUD --- ')
        # proc = Popen(['cmd','/C','maud_batch.bat'],cwd=outDir,creationflags=CREATE_NEW_CONSOLE,stdout=stdout,stderr=stderr)
        # proc.wait()        

        # # ====> 2nd Refinement Step
        # # ===> Lattice

        # # reload in after running refinement
        # firstPar = ParFile(newParFile)
        # # set initial refinable parameters
        # initPar.refine('lat_par', 'y') #shut lat param off      

        # # generate new par file to be used
        # newParFile       = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_Step2.par'.format(sample_desc))
        # newParResultFile = os.path.join(sampleDir,sample_desc,'refineResult_Step2.txt')
        # relParFile       = os.path.relpath(newParFile,outDir)
        # relParResult     = os.path.relpath(newParResultFile,outDir)
        # firstPar.save(newParFile)

        # # write batch file that MAUD can read in
        # # those paths need to be relative to the batchFile
        # firstPar._writeBat(batchFile, relParFile, 5, relParFile, relParResult)
        
        # # call the batch script in new CMD window
        # pbar.set_description('{} : Step 2'.format(sample))
        # pbar.refresh()
        # tqdm_write('\n --- Calling MAUD --- ')
        # proc = Popen(['cmd','/C','maud_batch.bat'],cwd=outDir,creationflags=CREATE_NEW_CONSOLE,stdout=stdout,stderr=stderr)
        # proc.wait()   

        # # ====> 3rd Refinement Step
        # # ===> Texture

        # # reload in after running refinement
        # secondPar = ParFile(newParFile)
        # # set initial refinable parameters
        # secondPar.refine('texture', 'y') 

        # # generate new par file to be used
        # newParFile       = os.path.join(sampleDir,sample_desc,'HB2C_MAUD_{}_Step3.par'.format(sample_desc))
        # newParResultFile = os.path.join(sampleDir,sample_desc,'refineResult_Step3.txt')
        # relParFile       = os.path.relpath(newParFile,outDir)
        # relParResult     = os.path.relpath(newParResultFile,outDir)
        # secondPar.save(newParFile)

        # # write batch file that MAUD can read in
        # # those paths need to be relative to the batchFile
        # secondPar._writeBat(batchFile, relParFile, 5, relParFile, relParResult)
            
        # # call the batch script in new CMD window
        # pbar.set_description('{} : Step 3'.format(sample))
        # pbar.refresh()
        # tqdm_write('\n --- Calling MAUD --- ')
        # proc = Popen(['cmd','/C','maud_batch.bat'],cwd=outDir,creationflags=CREATE_NEW_CONSOLE,stdout=stdout,stderr=stderr)
        # proc.wait()              
            
            
            
            
            
            
            
            
            
        