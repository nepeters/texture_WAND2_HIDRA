# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:21:37 2019

@author: npk
"""

"""
Peak fit WAND2 xye files

export pole figs for MTEX
"""

import logging
import os,sys,itertools,time
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from lmfit.model import ModelResult
from lmfit.models import PseudoVoigtModel, LinearModel, ConstantModel
from tqdm import tqdm
from tqdm.utils import _term_move_up
from datetime import datetime
from pyFAI.calibrant import Cell

import concurrent.futures as cf
import multiprocessing as mp

from ..fit_utils import rotate_project_q, export_pfs, write_MTEX

plt.ion()
plt.rcParams.update({
    "text.usetex": False,
    "axes.prop_cycle": mpl.cycler( color=['#e41a1c',
                                          '#377eb8',
                                          '#4daf4a',
                                          '#984ea3',
                                          '#ff7f00',
                                          '#ffff33',
                                          '#a65628',
                                          '#f781bf',
                                          '#999999'] )})

tqdm_colors = itertools.cycle(['RED','GREEN','YELLOW','BLUE','MAGENTA','CYAN','WHITE'])

def get_II(func, amp, cen, sig, frac, tt_ran, i):

    calc = func(tt_ran,amplitude=amp,center=cen,sigma=sig,fraction=frac)
    comp_II = np.trapz(y = calc,
                        x = tt_ran,
                        dx = tt_ran[1]-tt_ran[0])

    return comp_II, i

def runMP( samples, expInfo, rootDir, dataDir, mode, sample_rotations, n_tasks=None ):

    """
    run fitting routine - using concurrent.futures
    
    can't pass Cell object thru concurrent.futures..

    """

    # https://alexwlchan.net/2019/10/adventures-with-concurrent-futures/

    # number of tasks to be scheduled at once
    if n_tasks is None:
        # play it safe...
        n_tasks = mp.cpu_count() - 2

    assert mode.lower() in ['auto', 'texture'], 'mode must either be: \'auto\' or \'texture\''

    samples = iter(samples)

    with cf.ProcessPoolExecutor(max_workers=n_tasks) as executor:

        # initialize futures
        futures = {}

        # loop over slice of samples - based on n_tasks
        for pbaridx, sample in enumerate(itertools.islice(samples, n_tasks)):

            # get metadata for a given sample
            sampleInfo = expInfo.loc[expInfo['Sample'] == sample]
            start      = sampleInfo['Run #'].min()
            end        = sampleInfo['Run #'].max()

            # make this play nice with unix
            sample_desc = sample.replace(' #','_')
            if mode == 'texture': sample_desc = sample_desc+'_texred'
            sample_num  = int(sample_desc.split('_')[1])
            smplRot = sample_rotations[str(sample_num)]

            if sample_num > 6:
                rot_phase = 'FCC'
            else:
                rot_phase = 'BCC'    

            # submit job
            futures[sample] = executor.submit(fit_WAND2, 
                                              sample_desc, 
                                              rootDir, 
                                              dataDir, 
                                              sampleInfo, 
                                              start, end, 
                                              **{'mode':mode,'pbaridx':pbaridx,'pbarcolor':next(tqdm_colors),
                                                  'smplRot':smplRot,'ranges':{'BCC':'0:0.5:5.5','FCC':'0:0.5:4'},'rot_phase':rot_phase,
                                                  'logfile':os.path.join(rootDir,'fit_results',sample_desc,'fit_log.txt')})
    
            time.sleep(0.05)

        while futures:

            # Wait for the next future to complete
            done, futures = cf.wait(
                futures, return_when = cf.FIRST_COMPLETED
                )             
            
            for fut in done:
                print("The outcome is {}".format(fut.result()))            
                
            # Schedule the next set of futures. We don't want more than N futures
            # in the pool at a time, to keep memory consumption down.
            
            # loop over slice of samples - based on n_tasks
            for pbaridx, sample in enumerate(itertools.islice(samples, len(done))):

                # get metadata for a given sample
                sampleInfo = expInfo.loc[expInfo['Sample'] == sample]
                start      = sampleInfo['Run #'].min()
                end        = sampleInfo['Run #'].max()

                # make this play nice with unix
                sample_desc = sample.replace(' #','_')
                if mode == 'texture': sample_desc = sample_desc+'_texred'
                sample_num  = int(sample_desc.split('_')[1])
                smplRot = sample_rotations[str(sample_num)]

                if sample_num > 6:
                    rot_phase = 'FCC'
                else:
                    rot_phase = 'BCC' 

                # submit job
                futures[sample] = executor.submit(fit_WAND2, 
                                                sample_desc, 
                                                rootDir, 
                                                dataDir, 
                                                sampleInfo, 
                                                start, end, 
                                                **{'mode':mode,'pbaridx':pbaridx,'pbarcolor':next(tqdm_colors),
                                                  'smplRot':smplRot,'ranges':{'BCC':'0:0.5:5.5','FCC':'0:0.5:4'},'rot_phase':rot_phase,
                                                  'logfile':os.path.join(rootDir,'fit_results',sample_desc,'fit_log.txt')})

def calcCaglioti(theta, U, V, W):

    """
    calculate Caglioti broadening

    Gaussian Profile ~ U*tan2θ + V*tanθ + W 

    """

    return np.sqrt( U * np.tan(np.deg2rad(theta))**2 + V * np.tan(np.deg2rad(theta)) + W ) / 2.355

def xye_reader(file):

    """
    read .xye files in

    three columns [2theta(deg) intensity esd]
    no headers
    """
    header = ['2theta','int','esd']

    try:
        datain = pd.read_fwf(file, widths=[17, 17, 17], header=None, names=header, converters={h:float for h in header})
    except:
        datain = pd.read_fwf(file, widths=[17, 17, 17], skiprows=9, header=None, names=header, converters={h:float for h in header})

    return datain

def fit_WAND2(desc, rootDir, dataDir, expInfo, start, end, mode, liveplot=True, sequential=False, poleFig_type='jul', smplRot=None, ranges=None, rot_phase=None, pbaridx=None, pbarcolor='WHITE', logfile=None):
    
    """
    peak-fit WAND2 data -> pole figures

    uses Pseudo-Voigt function in lmfit
    provides a live fit-plot status 

    inputs:
        desc: sample description (str)
        rootDir: root directory
        dataDir: reduced data directory
        expInfo: experiment metadata (DataFrame)
        start: start run#
        end: end run#
        pf_dat: pole figure storage (dict)
        tt_all: two-theta list

    """

    if logfile:
        # get dir, name
        logfile_dir, _ = os.path.split(logfile)
        # check if dir exists
        if logfile_dir == '': pass
        elif os.path.exists(logfile_dir): pass
        else: os.makedirs(logfile_dir)
        # setup logger
        from imp import \
            reload  # python 2.x don't need to import reload, use it directly
        reload(logging)
        logging.basicConfig(filename=logfile, level=logging.INFO, filemode='w+',format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%d-%b-%y %H:%M:%S')

    inst = 'WAND2'

    # caglioti   = {'U':0.5328423,
    #               'V':-1.3918402,
    #               'W':2.4130297}

    WAND2_inst = {'Lam':1.4863984427800867,
                    'SH/L':0.10901188216118053,
                    'V':-388.52153929753064,
                    'Zero':-0.06586769734319739,
                    'U':2284.34769527403,
                    'W':-150.9677549520271,
                    'Azimuth':0.0,
                    'Y':0.0,
                    'X':0.0,
                    'Z':0.0,
                    'Type':'PNC',
                    'Bank':1.0,
                    'Polariz.':0.0}

    WAND2_bkgd = {'Back;0': 1.0,
                    'Back;1': 0.0,
                    'nDebye': 0,
                    'nPeaks': 0,
                    'Back File': '',
                    'BF mult': 1.0,
                    'Type': 'chebyschev-1'}

    # check pole figure type
    assert poleFig_type in ['jul', 'xyz'], 'Invalid pole figure type choice; either \'xyz\' or \'jul\''

    ## fitting setup ##
    d_all    = []
    ref_all  = []
    cnt_all  = []
    name_all = []

    ## get phase data ##
    for pi, (pn, ph) in enumerate(phases.items()):

        for k,v in ph.d_spacing(dmin=lmbda/2).items():

            d_all.append(v[0])
            ref_all.append(v[-1])
            cnt_all.append(pi)

        name_all.append(pn)

    sort_idx = np.argsort(d_all)
    d_all    = [d_all[i] for i in sort_idx]
    ref_all  = [ref_all[i] for i in sort_idx]
    cnt_all  = [cnt_all[i] for i in sort_idx]
    tt_all   = [2*np.rad2deg(np.arcsin(lmbda/(2*d))) for d in d_all]
    
    # number of measured q
    if mode == 'auto': dataCnt = (end - start)+1
    elif mode == 'texture': dataCnt = ((end - start)+1)*3

    ## setup pole fig dictionary ##
    pf_dat = {}
    for i,(d,ref,pi,tt) in enumerate(zip(d_all,ref_all,cnt_all,tt_all)):
        
        pf_dat[i+1] = {}
        pf_dat[i+1]['phase'] = name_all[pi]
        pf_dat[i+1]['ref']   = ''.join(map(str,ref))
        pf_dat[i+1]['data']  = np.zeros(( dataCnt, 5 ))
        pf_dat[i+1]['tt']    = tt
        pf_dat[i+1]['lattice'] = phases[name_all[pi]].lattice
        pf_dat[i+1]['lattice_type'] = phases[name_all[pi]].get_type()

        # for PF Δk index
        # will iterate +1 on each insertion
        # to account for variable # of points for each PF 
        # (shouldn't be the case in CW?)
        pf_dat[i+1]['pole_cnt'] = 0
        
        # setup flag if it was fit or not
        pf_dat[i+1]['fit'] = False    

    # where to store
    poleFig_path   = os.path.join(rootDir,'pole_figs',desc)
    fitResult_path = os.path.join(rootDir,'fit_results',desc,'params')
    fitImage_path  = os.path.join(rootDir,'fit_results',desc,'figures')

    if not os.path.exists(fitResult_path): os.makedirs(fitResult_path)
    if not os.path.exists(fitImage_path): os.makedirs(fitImage_path)
    if not os.path.exists(poleFig_path): os.makedirs(poleFig_path)

    # progress bar setup
    if pbaridx is None:
        refine_pbar = tqdm(np.linspace(start,end,num=(end-start)+1,dtype=int),
                           desc=desc,total=dataCnt,colour=pbarcolor, dynamic_ncols=True)
    else:
        refine_pbar = tqdm(np.linspace(start,end,num=(end-start)+1,dtype=int),
                           desc=desc,position=pbaridx,total=dataCnt,colour=pbarcolor, dynamic_ncols=True)       

    border = "="*80
    clear_border = _term_move_up() + "\r" + " "*len(border) + "\r"

    # if liveplot:
    #     # ## figure setup
    #     fig = plt.figure(figsize=(12.8,4.8),constrained_layout=True)
    #     gs  = fig.add_gridspec(5,4)
    #     ax1 = fig.add_subplot(gs[:4,:2])
    #     ax2 = fig.add_subplot(gs[:4,2:])
    #     ax3 = fig.add_subplot(gs[4,:2])

    # plt.pause(0.05)

    # counter..
    ri = 0

    for run in refine_pbar:
        
        t0 = time.time()

        runInfo = expInfo[expInfo['Run #'] == run] 
        chi     = 90 - runInfo['Chi'].to_numpy()[0] ## 90 - chi for original projection
        phi     = 360 - runInfo['Phi'].to_numpy()[0]
        omega   = 0 # not used            

        ## these will vary with eta tilt
        # set up a list
        diffData = []
        etaAll   = []

        if mode == 'auto':
            ## no in/out of plane reduction
            # loop over chi to get iterable
            for e in [0]:

                dataFileName = 'HB2C_{}.xye'.format(run)
                dataFile = os.path.join(dataDir,dataFileName)
                diffData.append(xye_reader(dataFile))
                etaAll.append(e) # not used

        elif mode == 'texture':
            ## in/out of plane reduced
            # loop over chi to get iterable
            for e in [-5,0,5]:

                # slightly different file name
                dataFileName = 'HB2C_{}_eta_{}-0.xye'.format(run,e)
                dataFile = os.path.join(dataDir,dataFileName)
                diffData.append(xye_reader(dataFile))
                etaAll.append(e) # not used 

        inner_iter = zip(diffData,etaAll)

        t1 = time.time()
    
        logging.info('setup time:{}'.format(t1-t0))

        for dd,eta in inner_iter:

            t2 = time.time()

            label = '{}_e{}'.format(run,eta)

            # get 2theta range
            tt_ran = dd['2theta'].to_numpy()
            # get intensity
            meas_int = dd['int'].to_numpy()

            # weights - for Rwp
            weights  = 1 / (dd['esd']**2)
            
            # find what peaks are present
            tt_mask  = (tt_all >= min(tt_ran)) * (tt_all <= max(tt_ran))
            tt_pres = list(itertools.compress(tt_all,tt_mask))
            # only these are present
            tt_pres_num = list(itertools.compress(range(len(tt_all)),tt_mask))
            # adjust index
            tt_pres_num = [v+1 for v in tt_pres_num]
        
            # num of peaks
            num_peaks = len(tt_pres_num)
            
            # setup model
            model = LinearModel()
            # model = ConstantModel()
            for i in tt_pres_num:
                # add individual peaks
                model = model + PseudoVoigtModel(prefix='p{}_'.format(i))
                
            # initialize params
            params = model.make_params()
            
            # # adjust background slope
            # params['slope'].set(value = )
            # params['intercept'].set(value = meas_int[-1] / 2.7)

            # guess the background
            I_bkgd = np.median(meas_int)

            # params['c'].set(value = I_bkgd)

            params['slope'].set(value=1E-3)
            params['intercept'].set(value = I_bkgd)

            # set peak initial parameters
            for i in tt_pres_num:
                
                # set center
                params['p{}_center'.format(i)].set(value = pf_dat[i]['tt'],
                                                    min = pf_dat[i]['tt'] - 0.2,
                                                    max = pf_dat[i]['tt'] + 0.2)
                
                # set amplitude
                I_meas = max(meas_int[np.argmin( np.abs( tt_ran - pf_dat[i]['tt'] ) )],1E-3)
                params['p{}_amplitude'.format(i)].set(I_meas, min=0.0)
                
                # if isinstance(caglioti,dict):
                #     # set lims on FWHM
                #     sig_guess = caglioti(pf_dat[i]['tt']/2,
                #                         caglioti['U'],
                #                         caglioti['V'],
                #                         caglioti['W'])
                # else:
                #     sig_guess = 0.2

                sig_guess = calcCaglioti(pf_dat[i]['tt']/2,
                                         caglioti['U'],
                                         caglioti['V'],
                                         caglioti['W'])

                params['p{}_sigma'.format(i)].set(value=sig_guess, min=0.5*sig_guess, max=2.0*sig_guess)

            if sequential:
                # skip on first run
                if run == start: pass
                else: 
                    priorFitResult = os.path.join(fitResult_path,
                                                'fitParams_{}.json'.format(prev_label))
                    with open(priorFitResult,'r') as f_in:
                        params = params.load(f_in)

                    # set peak parameters
                    for i in tt_pres_num:
                        
                        # set amplitude
                        I_meas = max(meas_int[np.argmin( np.abs( tt_ran - pf_dat[i]['tt'] ) )],1E-3)
                        params['p{}_amplitude'.format(i)].set(I_meas, min=0)

            t3 = time.time()

            logging.info('model setup time:{}'.format(t3-t2))

            # fit model   
            init       = model.eval(params, x=tt_ran)
            try:
                out = model.fit(meas_int,
                                params,
                                x=tt_ran,                                         
                                fit_kws={'gtol':1E-2,
                                         'xtol':1E-2,
                                         'ftol':1E-2},
                                method='least_squares') 
            except:
                continue
            comps      = out.eval_components(x=tt_ran)

            t4 = time.time()
            
            logging.info('model fit time:{}'.format(t4-t3))

            out_pars   = out.params.copy()
            n_boot     = 5000
            II         = {}
            II_esd     = {}

            # # Get uncertainty estimate for integrated intensity (?)
            # for comp in out.model.components:
            #     if 'linear' in comp.name: continue
            #     # Get the names and params
            #     comp_par_names = comp._param_names
            #     comp_pars = []
            #     for par_name in comp_par_names:
            #         par = out_pars[par_name]
            #         if par.stderr is None:
            #             comp_pars.append(np.ones(n_boot)*par.value)
            #             # tqdm.write(str(par))
            #         else:
            #             try:
            #                 comp_pars.append(norm.rvs(loc=par.value,scale=par.stderr,size=n_boot))
            #             except ValueError:
            #                 comp_pars.append(np.ones(n_boot)*par.value)
            #     comp_pars = np.asarray(comp_pars).T
            
            #     comp_II = []

            #     for n in range(n_boot):
            #         # Evaluate the new set
            #         calc = comp.func(tt_ran,amplitude=comp_pars[n,0],center=comp_pars[n,1],sigma=comp_pars[n,2],fraction=comp_pars[n,3])
            #         comp_II.append(np.trapz(y = calc,
            #                                 x = tt_ran,
            #                                 dx = tt_ran[1]-tt_ran[0]))
                
            #     comp_II = removeOutliers(comp_II, 1.5)
            #     II[comp.prefix]     = np.mean(comp_II)
            #     II_esd[comp.prefix] = np.std(comp_II)

            # Get uncertainty estimate for integrated intensity - fast way, just use cov
            for comp in out.model.components:
                if 'linear' in comp.name: continue
                elif 'constant' in comp.name: continue
                comp_par_names = comp._param_names
                # II[comp.prefix]     = np.mean(out.params[comp_par_names[0]].stderr)
                esd = out.params[comp_par_names[0]].stderr
                if esd is None:
                    II_esd[comp.prefix] = 0.0
                elif np.isnan(esd) == False:
                    II_esd[comp.prefix] = esd
                else:
                    II_esd[comp.prefix] = 0.0

            prev_label = label

            # calculate weighted R (fit quality)
            rwp = np.sum( weights * out.residual**2 ) / np.sum( weights * meas_int**2 )
            # write to console
            logging.info(out.fit_report())
            if pbaridx is None:
                # running with concurrent.futures.. don't print all this
                logging.info('--- Run #{} | Eta:{} ---'.format(run,eta))
                if out.success: logging.info('Fit was successful!')
                logging.info('Rwp : {:3.2f}%'.format(rwp*100))
                logging.info(border)
                refine_pbar.update(1)
            else:
                refine_pbar.set_description('{}|Run #{}|Eta:{}|Rwp:{:3.2f}%'.format(desc,run,eta,rwp*100))
                refine_pbar.refresh()
                refine_pbar.update(1)
            
            # save fit params for posterity
            fitResult = os.path.join(fitResult_path,'fitParams_{}.json'.format(label))
            with open(fitResult,'w') as f_out:
                out.params.dump(f_out)            

            t5 = time.time()

            logging.info('model output time:{}'.format(t5-t4))

            # store peak intensity
            for i in tt_pres_num:
            
                # get q counter        
                pole_cnt = pf_dat[i]['pole_cnt']
                # get 2theta
                tt = out.params['p{}_center'.format(i)].value
                # get projection (q)
                q = rotate_project_q(tt/2, omega, chi, phi, 270 - eta)
                # store it
                pf_dat[i]['data'][pole_cnt,0] = q[0]
                pf_dat[i]['data'][pole_cnt,1] = q[1]
                pf_dat[i]['data'][pole_cnt,2] = q[2]
                
                # tell me it's fit
                pf_dat[i]['fit'] = True

                # tell me what type to output
                pf_dat[i]['type'] = poleFig_type

                # calc integ. int
                II = np.trapz(y = comps['p{}_'.format(i)],
                              x = tt_ran,
                              dx = tt_ran[1]-tt_ran[0])
            
                # store integ. int
                pf_dat[i]['data'][pole_cnt,3] = II
                pf_dat[i]['data'][pole_cnt,4] = II_esd['p{}_'.format(i)]

                ## counter for Δk
                pf_dat[i]['pole_cnt'] += 1
            
            t6 = time.time()

            logging.info('pole fig save time:{}'.format(t6-t5))

            if liveplot is True:
            
                # ## figure setup
                fig = plt.figure(figsize=(12.8,4.8),constrained_layout=True)
                gs  = fig.add_gridspec(5,4)
                ax1 = fig.add_subplot(gs[:4,:2])
                ax2 = fig.add_subplot(gs[:4,2:])
                ax3 = fig.add_subplot(gs[4,:2])
            
                ## print result plot      
                ax1.plot(tt_ran, meas_int, 'b')
                ax1.plot(tt_ran, init, 'k--', label='initial fit')
                ax1.plot(tt_ran, out.best_fit, 'r-', label='best fit')
                ax3.plot(tt_ran, out.best_fit - meas_int, 'g-')
                ax2.plot(tt_ran, meas_int, 'b')
            
                for i in tt_pres_num:
                    
                    ax2.plot(tt_ran, comps['p{}_'.format(i)]+comps['linear'], '--', label='Peak {}_{}'.format(pf_dat[i]['phase'],pf_dat[i]['ref']))
                            
                # housekeeping
                ax1.legend(loc='best')
                if num_peaks < 7: ax2.legend(loc='best')
                # ax2.legend(loc='best')
                ax1.set_ylabel('Intensity')
                ax1.set_xlabel('2θ (degrees)')
                ax2.set_ylabel('Intensity')
                ax2.set_xlabel('2θ (degrees)')
                ax3.set_ylabel('Difference')
                ax3.set_xlabel('2θ (degrees)')

                ax2.set_ylim(top=0.45*np.max(meas_int))

                ## save fit image for posterity
                plt.savefig(os.path.join(fitImage_path,'fit_{}.png'.format(label)))
                plt.close()
            
                ri += 1 

            t7 = time.time()

            logging.info('plot save time:{}'.format(t7-t6))

    ## close out
    if liveplot: plt.close()
    
    # export the pole figures
    export_pfs(inst, desc, pf_dat, poleFig_path)
    # # export the MTEX file
    # write_MTEX(desc, pf_dat, poleFig_path, smplSym='1', smplRot=smplRot, ranges=ranges, rot_phase=rot_phase)

    return pf_dat

def load_WAND2(desc, rootDir, dataDir, expInfo, start, end, mode, exportPFs=False, poleFig_type='jul', smplRot=None, ranges=None, rot_phase=None):
    
    """
    load a prior fit from .json files
    
    peak-fit WAND2 data -> pole figures
    uses Pseudo-Voigt function in lmfit
    
    inputs:
        desc: sample description (str)
        rootDir: root directory
        dataDir: reduced data directory
        expInfo: experiment metadata (DataFrame)
        start: start run#
        end: end run#
        pf_dat: pole figure storage (dict)
        tt_all: two-theta list
     
    """

    inst = 'WAND2'

    # check pole figure type
    assert poleFig_type in ['jul', 'xyz'], 'Invalid pole figure type choice; either \'xyz\' or \'jul\''

    ## fitting setup ##
    d_all    = []
    ref_all  = []
    cnt_all  = []
    name_all = []

    ## get phase data ##
    for pi, (pn, ph) in enumerate(phases.items()):

        for k,v in ph.d_spacing(dmin=lmbda/2).items():

            d_all.append(v[0])
            ref_all.append(v[-1])
            cnt_all.append(pi)

        name_all.append(pn)

    sort_idx = np.argsort(d_all)
    d_all    = [d_all[i] for i in sort_idx]
    ref_all  = [ref_all[i] for i in sort_idx]
    cnt_all  = [cnt_all[i] for i in sort_idx]
    tt_all   = [2*np.rad2deg(np.arcsin(lmbda/(2*d))) for d in d_all]
    
    # number of measured q
    if mode == 'auto': dataCnt = (end - start)+1
    elif mode == 'texture': dataCnt = ((end - start)+1)*3

    ## setup pole fig dictionary ##
    pf_dat = {}
    for i,(d,ref,pi,tt) in enumerate(zip(d_all,ref_all,cnt_all,tt_all)):
        
        pf_dat[i+1] = {}
        pf_dat[i+1]['phase'] = name_all[pi]
        pf_dat[i+1]['ref']   = ''.join(map(str,ref))
        pf_dat[i+1]['data']  = np.zeros(( dataCnt, 5 ))
        pf_dat[i+1]['tt']    = tt
        pf_dat[i+1]['lattice'] = phases[name_all[pi]].lattice
        pf_dat[i+1]['lattice_type'] = phases[name_all[pi]].get_type()

        # for PF Δk index
        # will iterate +1 on each insertion
        # to account for variable # of points for each PF 
        # (shouldn't be the case in CW?)
        pf_dat[i+1]['pole_cnt'] = 0
        
        # setup flag if it was fit or not
        pf_dat[i+1]['fit'] = False        
    
    # where to store
    poleFig_path   = os.path.join(rootDir,'pole_figs',desc)
    fitResult_path = os.path.join(rootDir,'fit_results',desc,'params')
    fitImage_path  = os.path.join(rootDir,'fit_results',desc,'figures')

    if not os.path.exists(fitResult_path): os.makedirs(fitResult_path)
    if not os.path.exists(fitImage_path): os.makedirs(fitImage_path)
    if not os.path.exists(poleFig_path): os.makedirs(poleFig_path)

    # progress bar setup
    refine_pbar = tqdm(np.linspace(start,end,num=(end-start)+1,dtype=int),
                       desc=desc)   

    liveplot = False

    if liveplot:
        # ## figure setup
        fig = plt.figure(figsize=(12.8,4.8),constrained_layout=True)
        gs  = fig.add_gridspec(5,4)
        ax1 = fig.add_subplot(gs[:4,:2])
        ax2 = fig.add_subplot(gs[:4,2:])
        ax3 = fig.add_subplot(gs[4,:2])

    plt.pause(0.05)

    # counter..
    ri = 0           

    for run in refine_pbar:
        
        runInfo = expInfo[expInfo['Run #'] == run] 
        chi     = 90 - runInfo['Chi'].to_numpy()[0] ## 90 - chi for original projection
        phi     = 360 - runInfo['Phi'].to_numpy()[0]
        omega   = 0 # not used            

        ## these will vary with eta tilt
        # set up a list
        diffData = []
        etaAll   = []

        if mode == 'auto':
            ## no in/out of plane reduction
            # loop over chi to get iterable
            for e in [0]:

                dataFileName = 'HB2C_{}.xye'.format(run)
                dataFile = os.path.join(dataDir,dataFileName)
                diffData.append(xye_reader(dataFile))
                etaAll.append(e) # not used

        elif mode == 'texture':
            ## in/out of plane reduced
            # loop over chi to get iterable
            for e in [-5,0,5]:

                # slightly different file name
                dataFileName = 'HB2C_{}_eta_{}-0.xye'.format(run,e)
                dataFile = os.path.join(dataDir,dataFileName)
                diffData.append(xye_reader(dataFile))
                etaAll.append(e) # not used 

        inner_iter = zip(diffData,etaAll)
    
        for dd,eta in inner_iter:
        
            label = '{}_e{}'.format(run,eta)

            # get 2theta range
            tt_ran = dd['2theta'].to_numpy()
            # get intensity
            meas_int = dd['int'].to_numpy()
            # weights - for Rwp
            weights  = 1 / (dd['esd']**2)
            
            # find what peaks are present
            tt_mask  = (tt_all >= min(tt_ran)) * (tt_all <= max(tt_ran))
            tt_pres = list(itertools.compress(tt_all,tt_mask))
            # only these are present
            tt_pres_num = list(itertools.compress(range(len(tt_all)),tt_mask))
            # adjust index
            tt_pres_num = [v+1 for v in tt_pres_num]
        
            # num of peaks
            num_peaks = len(tt_pres_num)
            
            # setup model
            # model = LinearModel()
            model = ConstantModel()
            for i in tt_pres_num:
                # add individual peaks
                model = model + PseudoVoigtModel(prefix='p{}_'.format(i))
                
            # initialize params
            params = model.make_params()
            
            # load old parameters from json files
            fitResult = os.path.join(fitResult_path,
                                    'fitParams_{}.json'.format(label))
            
            with open(fitResult,'r') as f_in:
                try:
                    params = params.load(f_in)
                    out    = ModelResult(model,params)
                    comps  = out.eval_components(x=tt_ran)
                    best   = out.eval(x=tt_ran)
                except:
                    pass  

            out_pars   = out.params.copy()
            n_boot     = 5000
            II         = {}
            II_esd     = {}

            # # Get uncertainty estimate for integrated intensity (?)
            # for comp in out.model.components:
            #     if 'linear' in comp.name: continue
            #     elif 'constant' in comp.name: continue
            #     # Get the names and params
            #     comp_par_names = comp._param_names
            #     comp_pars = []
            #     for par_name in comp_par_names:
            #         par = out_pars[par_name]
            #         if par.stderr is None:
            #             comp_pars.append(np.ones(n_boot)*par.value)
            #             # tqdm.write(str(par))
            #         else:
            #             try:
            #                 comp_pars.append(norm.rvs(loc=par.value,scale=par.stderr,size=n_boot))
            #             except ValueError:
            #                 comp_pars.append(np.ones(n_boot)*par.value)
            #     comp_pars = np.asarray(comp_pars).T
            #     tt_ran2 = np.tile(tt_ran, [5000,1])
            #     calc    = comp.func(tt_ran2, comp_pars[:,0][:,None],comp_pars[:,1][:,None],comp_pars[:,2][:,None],comp_pars[:,3][:,None])
            #     comp_II = np.trapz(calc, x=tt_ran2, dx=tt_ran[1]-tt_ran[0])


            #     # for n in range(n_boot):
            #     #     # Evaluate the new set
            #     #     calc = comp.func(tt_ran,amplitude=comp_pars[n,0],center=comp_pars[n,1],sigma=comp_pars[n,2],fraction=comp_pars[n,3])
            #     #     comp_II.append(np.trapz(y = calc,
            #     #                             x = tt_ran,
            #     #                             dx = tt_ran[1]-tt_ran[0]))
                
            #     comp_II = removeOutliers(comp_II, 1.5)
            #     II[comp.prefix]     = np.mean(comp_II)
            #     II_esd[comp.prefix] = np.std(comp_II)

            #     # if np.mean(comp_II) > 1E6:
            #     #     comp_pars = []
            #     #     for par_name in comp_par_names:
            #     #         comp_pars.append(out_pars[par_name].value) 
            #     #         print(f'{out_pars[par_name].value} - {out_pars[par_name].stderr}')
            #     #     calc = comp.func(tt_ran,amplitude=comp_pars[0],center=comp_pars[1],sigma=comp_pars[2],fraction=comp_pars[3])
            #     #     # fig, ax = plt.subplots(1)
            #     #     # ax.plot(tt_ran, calc)
            #     #     # ax.plot(tt_ran, meas_int)
            #     #     print(f'{comp_pars[1]/2} - {run} - {chi} - {phi} - {eta}')
            #     #     return comp_II, np.trapz(y = calc,x = tt_ran,dx = tt_ran[1]-tt_ran[0])

            # store peak intensity
            for i in tt_pres_num:
            
                # get q counter        
                pole_cnt = pf_dat[i]['pole_cnt']
                # get 2theta
                tt = out.params['p{}_center'.format(i)].value
                # get projection (q)
                q = rotate_project_q(tt/2, omega, chi, phi, 270 - eta)
                print(out.params[f'p{i}_center'])
                print(pf_dat[i]['ref'])
                print(f'{tt/2} - {omega} - {chi} - {phi} - {270-eta}')
                print(q)

                # store it
                pf_dat[i]['data'][pole_cnt,0] = q[0]
                pf_dat[i]['data'][pole_cnt,1] = q[1]
                pf_dat[i]['data'][pole_cnt,2] = q[2]
                
                # tell me it's fit
                pf_dat[i]['fit'] = True

                # tell me what type to output
                pf_dat[i]['type'] = poleFig_type

                # calc integ. int
                II = np.trapz(y = comps['p{}_'.format(i)],
                              x = tt_ran,
                              dx = tt_ran[1]-tt_ran[0])
            
                # # store integ. int
                # pf_dat[i]['data'][pole_cnt,3] = II
                # pf_dat[i]['data'][pole_cnt,4] = II_esd['p{}_'.format(i)]

                # store integ. int
                pf_dat[i]['data'][pole_cnt,3] = II
                pf_dat[i]['data'][pole_cnt,4] = 0.0
                
                ## counter for Δk
                pf_dat[i]['pole_cnt'] += 1

    #         if liveplot is True:
            
    #             if ri > 0:
    #                 ax1.clear()
    #                 ax2.clear()
    #                 ax3.clear()
            
    #             ## print result plot      
    #             ax1.plot(tt_ran, meas_int, 'b')
    #             # ax1.plot(tt_ran, init, 'k--', label='initial fit')
    #             ax1.plot(tt_ran, best, 'r-', label='best fit')
    #             ax3.plot(tt_ran, best - meas_int, 'g-')
    #             ax2.plot(tt_ran, meas_int, 'b')
            
    #             for i in tt_pres_num:
                    
    #                 ax2.plot(tt_ran, comps['p{}_'.format(i)]+comps['constant'], '--', label='Peak {}_{}'.format(pf_dat[i]['phase'],pf_dat[i]['ref']))
                            
    #             # housekeeping
    #             ax1.legend(loc='best')
    #             if num_peaks < 7: ax2.legend(loc='best')
    #             # ax2.legend(loc='best')
    #             ax1.set_ylabel('Intensity')
    #             ax1.set_xlabel('2θ (degrees)')
    #             ax2.set_ylabel('Intensity')
    #             ax2.set_xlabel('2θ (degrees)')
    #             ax3.set_ylabel('Difference')
    #             ax3.set_xlabel('2θ (degrees)')
            
    #             plt.pause(0.05)       
    #             plt.show()    
            
    #             ## save fit image for posterity
    #             plt.savefig(os.path.join(fitImage_path,'fit_{}.pdf'.format(label)))
            
    #             ri += 1 

    # ## close out
    # if liveplot: plt.close()

    # export the pole figures
    if exportPFs: export_pfs(inst, desc, pf_dat, poleFig_path)
    
    # write_MTEX(desc, pf_dat, poleFig_path, smplSym='1', smplRot=smplRot, ranges=ranges, rot_phase=rot_phase)

    return pf_dat

rootDir = '/' #

instName = 'WAND2'
mode     = 'auto' ## 'auto' = autoreduced | 'texture' = in/out of plane separate
lmbda    = 1.488 #angstroms

""" 
needs columns: 
[Time, Run #, Chi, Phi, Sample ]  
"""

# reduced XYE files
tex_dataDir     = os.path.join(rootDir,'texturereduced')
aut_dataDir     = os.path.join(rootDir,'autoreduced')
expInfoFile = os.path.join(rootDir,'expSummary_IPTS-22335.csv')
expInfo     = pd.read_csv(expInfoFile)

phases = {'FCC':Cell.cubic(3.5938,lattice_type='F'),
          'BCC':Cell.cubic(2.8728,lattice_type='I')}

caglioti   = {'U':0.5328423,
              'V':-1.3918402,
              'W':2.4130297}

# sample rotations
rotations = {'9':[('z',105), ('y',90)],
             '6':[('z',0),   ('y',90)],
             '3':[('z',25),  ('y',90)]} 

## setup run sets for fitting
samples = pd.unique(expInfo['Sample'])

realSamples = []

for sample in samples:
   
    if 'Steel' in sample:
    
        realSamples.append(sample)

# if __name__ == '__main__':
#     # run through texture-reduced
#     runMP(realSamples, expInfo, rootDir, tex_dataDir, mode='texture', sample_rotations=rotations)
#     # run through auto-reduced
#     runMP(realSamples, expInfo, rootDir, aut_dataDir, mode='auto', sample_rotations=rotations)

for sample in samples:
   
    if 'Steel' in sample:
    
        # get metadata for a given sample
        sampleInfo = expInfo.loc[expInfo['Sample'] == sample]
        start      = sampleInfo['Run #'].min()
        end        = sampleInfo['Run #'].max()
            
        # make this play nice with unix
        sample_desc = sample.replace(' #','_')
        if mode == 'texture': sample_desc = sample_desc+'_texred'
        sample_num  = int(sample_desc.split('_')[1])
        smplRot = rotations[str(sample_num)]

        if sample_desc == 'Steel_6':

            if sample_num > 6:
                rot_phase = 'FCC'
            else:
                rot_phase = 'BCC'

            # run fit routine
            fit_WAND2(sample_desc, rootDir, aut_dataDir, sampleInfo, start, end, 'auto', smplRot=smplRot, 
                            ranges={'BCC':'0:0.5:5.5','FCC':'0:0.5:4'},rot_phase=rot_phase, logfile=os.path.join(rootDir,'fit_results',sample_desc,'fit_log.txt'))

            sample_desc = sample_desc+'_texred'
            fit_WAND2(sample_desc, rootDir, tex_dataDir, sampleInfo, start, end, 'texture', smplRot=smplRot,
                            ranges={'BCC':'0:0.5:5.5','FCC':'0:0.5:4'},rot_phase=rot_phase, logfile=os.path.join(rootDir,'fit_results',sample_desc,'fit_log.txt'))
