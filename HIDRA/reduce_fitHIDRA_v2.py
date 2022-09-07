# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:21:37 2019

@author: npk
"""

"""
single peak fits -> pole figure
 
from HIDRA h5 files
"""

import os, sys, itertools, h5py, time
import numpy as np
import numpy.ma as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from lmfit.model import ModelResult
from lmfit.models import PseudoVoigtModel, ConstantModel
from tqdm import tqdm
from tqdm.utils import _term_move_up
from datetime import datetime
from pyFAI.calibrant import Cell
from scipy.spatial.transform import Rotation as R

import concurrent.futures as cf
import multiprocessing as mp

from ..fit_utils import rotate_project_q, export_pfs, write_MTEX 

now = datetime.now()

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
                                          '#999999'])})

"""
h5 file structure

** - data is used

*.h5
├── instrument
│   ├── calibration
│   ├── monochomator setting
│   │   └── wave length **
│   ├── efficiency calibration
│   └── geometry setup 
├── raw data
│   ├── sub-runs
│   └── logs
│       ├── chi    **
│       ├── phi    **
│       └── 2thetaSetpoint **
├── mask
├── reduced diffraction data
│   ├── eta_-5.0 **
│   ├── eta_0.0  **
│   ├── eta_5.0  **
│   └── 2theta   **
└── peaks

pf data-array format:

   [:,:,n] pf # 
  ↗
  → [:,n,:] phi, chi, integ. int 
 ↓
 [n,:,:]
 #
 u
 n
 i
 q
 u
 e
 
 Δ
 k
 
"""
def removeOutliers(x, outlierConstant):
    a = np.array(x)
    upper_quartile = np.percentile(a, 75)
    lower_quartile = np.percentile(a, 25)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile - IQR, upper_quartile + IQR)
    resultList = []
    for y in a.tolist():
        if y >= quartileSet[0] and y <= quartileSet[1]:
            resultList.append(y)
    return resultList

def fit_HIDRA(runNumber, rootDir, dataDir, phases, mode='texture', sequential=False, liveplot=True, exportPFs=False, pfType='jul', smplRot=None, pbaridx=None, pbarcolor='WHITE', ranges=None, rot_phase=None):

    """
    peak-fit HIDRA data -> pole figs

    """

    # define instrument
    inst  = 'HIDRA'

    # load in .h5 file
    if mode == 'auto': 
        fname = 'HB2B_{}.h5'.format(runNumber)
        desc  = '{}_aut_uq'.format(runNumber)
    elif mode == 'texture': 
        fname = 'HB2B_{}_texture.h5'.format(runNumber)
        desc  = '{}_tex_uq'.format(runNumber)
    else: raise ValueError('mode not recognized..')
    
    exp_h5 = h5py.File(os.path.join(dataDir,fname), 'r')

    # read wavelength
    lmbda = exp_h5['instrument/monochromator setting/wave length'][()][0]

    # read angular data
    chi       = exp_h5['raw data/logs/chi'][()]
    phi       = exp_h5['raw data/logs/phi'][()]
    omega     = exp_h5['raw data/logs/omega'][()]
    two_theta = exp_h5['reduced diffraction data/2theta'][()]

    # number of measured patterns (for loop)
    meas_num  = len(phi)

    # read intensity data
    if mode == 'auto': #no eta slice
        
        # get from raw data/logs/2thetaSetpoint
        num_det_pos   = len(np.unique(exp_h5['raw data/logs/2thetaSetpoint'][()]))
        num_eta_slice = 1

        eta_zero = np.nan_to_num(exp_h5['reduced diffraction data/main'][()])
        max_int  = np.max(eta_zero)

    elif mode == 'texture': #should have eta slices

        num_det_pos   = len(np.unique(exp_h5['raw data/logs/2thetaSetpoint'][()]))
        num_eta_slice = 3

        eta_neg5 = np.nan_to_num(exp_h5['reduced diffraction data/eta_-5.0'][()])
        eta_zero = np.nan_to_num(exp_h5['reduced diffraction data/eta_0.0'][()])
        eta_pos5 = np.nan_to_num(exp_h5['reduced diffraction data/eta_5.0'][()])

        max_int = np.max([np.max(eta) for eta in [eta_neg5, eta_zero, eta_pos5]])

    # close the h5 file
    exp_h5.close()

    # number of measured q
    rot_num   = int((meas_num/num_det_pos)*num_eta_slice)

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
    
    ## setup pole fig dictionary ##
    pfd = {}
    for i,(d,ref,pi,tt) in enumerate(zip(d_all,ref_all,cnt_all,tt_all)):
        
        pfd[i+1] = {}
        pfd[i+1]['phase'] = name_all[pi]
        pfd[i+1]['ref']   = ''.join(map(str,ref))
        pfd[i+1]['data']  = np.zeros(( rot_num, 5 ))
        pfd[i+1]['tt']    = tt
        pfd[i+1]['lattice'] = phases[name_all[pi]].lattice
        pfd[i+1]['lattice_type'] = phases[name_all[pi]].get_type()

        # for PF Δk index
        # will iterate +1 on each insertion
        # to account for variable # of points for each PF 
        # (shouldn't be the case in CW?)
        pfd[i+1]['pole_cnt'] = 0
        
        # setup flag if it was fit or not
        pfd[i+1]['fit'] = False

    # where to store
    poleFig_path   = os.path.join(rootDir,'pole_figs',desc)
    fitResult_path = os.path.join(rootDir,'fit_results',desc,'params')
    fitImage_path  = os.path.join(rootDir,'fit_results',desc,'figures')

    if not os.path.exists(fitResult_path): os.makedirs(fitResult_path)
    if not os.path.exists(fitImage_path): os.makedirs(fitImage_path)
    if not os.path.exists(poleFig_path): os.makedirs(poleFig_path)

    # progress bar setup
    if pbaridx is None:
        refine_pbar = tqdm(range(meas_num),desc=desc)
    else:
        refine_pbar = tqdm(range(meas_num),desc=desc, position=pbaridx)
        
    border = "="*80
    clear_border = _term_move_up() + "\r" + " "*len(border) + "\r"

    liveplot = False

    # ## figure setup
    # if liveplot is True:
        
    #     fig = plt.figure(figsize=(12.8,4.8),constrained_layout=True)
    #     gs  = fig.add_gridspec(5,4)
    #     ax1 = fig.add_subplot(gs[:4,:2])
    #     ax2 = fig.add_subplot(gs[:4,2:])
    #     ax3 = fig.add_subplot(gs[4,:2])
    #     plt.pause(0.05)

    k = 0

    ## loop over rotations
    for ri in refine_pbar:

        t0 = time.time()

        # easy to reference these later        
        o = omega[ri]
        c = 90 - chi[ri]
        p = 360 - phi[ri]
        
        if mode == 'auto': inner_iter = zip([eta_zero],[0])
        elif mode == 'texture': inner_iter = zip([eta_neg5, eta_zero, eta_pos5],[-5, 0, 5])
        # inner_iter = zip([eta_neg5, eta_zero, eta_pos5],[-5, 0, 5])

        # loop over data
        for meas_int,eta in inner_iter:

            # refine_pbar.write('\n')

            t2 = time.time()

            # if o*2 < 90: continue

            counter = 0

            label = 'tt{}_o{}_c{}_p{}_e{}'.format(round(o*2),round(o),round(c),round(p),270 - eta)
        
            # get mask on invalid data on edges
            valid_data = ma.masked_where(meas_int[ri,:]==0,meas_int[ri,:])
            valid      = ~valid_data.mask

            # get 2theta range of measurement
            tt_ran = two_theta[ri,valid]

            # get weights
            weights = 1 / meas_int[ri,valid]**2
            # get intensity
            inten = meas_int[ri,valid]
            
            # find what peaks are present
            tt_mask  = (tt_all >= min(tt_ran)) * (tt_all <= max(tt_ran))
            tt_pres = list(itertools.compress(tt_all,tt_mask))
            # only these are present
            tt_pres_num = list(itertools.compress(range(len(tt_all)),tt_mask))
            # adjust index
            tt_pres_num = [v+1 for v in tt_pres_num]
            
            # num of peaks
            num_peaks = len(tt_pres_num)

            # setup lmfit model
            model = ConstantModel()
            for i in tt_pres_num:
                # add individual peaks
                model = model + PseudoVoigtModel(prefix='p{}_'.format(i))
                
            ## initialize params
            params = model.make_params()
            
            # guess the background
            I_bkgd = np.median(inten)
            params['c'].set(value = I_bkgd)

            # set peak initial parameters
            for i in tt_pres_num:
                

                
                pk_loc = pfd[i]['tt']
                pk_loc_lo = pfd[i]['tt'] - 0.5
                pk_loc_hi = pfd[i]['tt'] + 0.5

                loi = np.argmin( np.abs( tt_ran - pk_loc_lo ) )
                hii = np.argmin( np.abs( tt_ran - pk_loc_hi ) )

                I_guess = (np.max(inten[loi:hii]) - I_bkgd)/2
                if I_guess < 4:
                    I_guess = 1E-2
                TT_guess = tt_ran[np.argmax(inten[loi:hii])+loi]

                # set center
                params['p{}_center'.format(i)].set(value = TT_guess,
                                                    min = TT_guess - 0.5,
                                                    max = TT_guess + 0.5)
                # set amplitude
                # print(f'{i} - {pk_loc}:{TT_guess} - {I_guess}')
                # print(f'{pk_loc_lo} - {pk_loc_hi}')
                # print(f'{i} - {inten[loi:hii]}')
                params['p{}_amplitude'.format(i)].set(I_guess, min=0)
                
                # set lims on FWHM
                params['p{}_sigma'.format(i)].set(value=0.2,min=0,max=0.35)

            # setup file to save parameters (.json)
            fitResult = os.path.join(fitResult_path,'fitParams_{}.json'.format(label))

            if sequential:
                # skip on first run
                if counter == 0: pass
                else: 
                    priorFitResult = os.path.join(fitResult_path,
                                                  'fitParams_{}.json'.format(prev_label))
                    with open(priorFitResult,'r') as f_in:
                        params = params.load(f_in)
            
            # fit model

            t3 = time.time()

            # refine_pbar.write('model setup time:{}'.format(t3-t2))   

            init       = model.eval(params, x=tt_ran)
            out        = model.fit(meas_int[ri, valid],
                            params,
                            x=tt_ran,                                         
                            fit_kws={'gtol':1E-3,
                                        'xtol':1E-3,
                                        'ftol':1E-3},
                            method='least_squares') 

            comps      = out.eval_components(x=tt_ran)

            t4 = time.time()
            
            # refine_pbar.write('model fit time:{}'.format(t4-t3))

            out_pars   = out.params.copy()
            n_boot     = 100
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
            #     tt_ran2 = np.tile(tt_ran, [n_boot,1])
            #     calc    = comp.func(tt_ran2, comp_pars[:,0][:,None],comp_pars[:,1][:,None],comp_pars[:,2][:,None],comp_pars[:,3][:,None])
            #     comp_II = np.trapz(calc, x=tt_ran2, dx=tt_ran[1]-tt_ran[0])


            #     # comp_pars = np.asarray(comp_pars).T

            #     # comp_II = []

            #     # for n in range(n_boot):
            #     #     # Evaluate the new set
            #     #     calc = comp.func(tt_ran,amplitude=comp_pars[n,0],center=comp_pars[n,1],sigma=comp_pars[n,2],fraction=comp_pars[n,3])
            #     #     comp_II.append(np.trapz(y = calc,
            #     #                             x = tt_ran,
            #     #                             dx = tt_ran[1]-tt_ran[0]))
                
            #     comp_II = removeOutliers(comp_II, 1.5)
            #     II[comp.prefix]     = np.mean(comp_II)
            #     II_esd[comp.prefix] = np.std(comp_II)

            #     # esd = out.params[comp_par_names[0]].stderr
            #     # print(f'{II[comp.prefix]} - {II_esd[comp.prefix]} - {esd} | {comp.prefix}')

            # Get uncertainty estimate for integrated intensity - fast way, just use cov
            for comp in out.model.components:
                if 'linear' in comp.name: continue
                elif 'constant' in comp.name: continue
                comp_par_names = comp._param_names
                # II[comp.prefix]     = np.mean(out.params[comp_par_names[0]].stderr)
                esd = out.params[comp_par_names[0]].stderr
                if esd is None:
                    II_esd[comp.prefix] = 0.0
                elif np.isnan(esd) is False:
                    II_esd[comp.prefix] = esd
                else:
                    II_esd[comp.prefix] = 0.0

            prev_label = label

            # calculate weighted R (fit quality)
            rwp = np.sum( weights * out.residual**2 ) / np.sum( weights * inten**2 )

            # write to console
            # this goes fast.. only print if there's a problem
            if not out.success: 
                refine_pbar.write(clear_border + '--- ω:{} | χ:{} | φ:{} | η:{} ---'.format(int(o),int(c),int(p),int(eta)))
                refine_pbar.update()
                refine_pbar.write(clear_border + 'Fit was not successful!')
                refine_pbar.update()
                refine_pbar.write(clear_border + 'Rwp : {:3.2f}%'.format(rwp*100))
                refine_pbar.update()
                refine_pbar.write(border)
                refine_pbar.update()

            # save fit params for posterity
            with open(fitResult,'w') as f_out:
                out.params.dump(f_out)   

            t5 = time.time()

            # refine_pbar.write('model output time:{}'.format(t5-t4))

            # store peak intensity
            for i in tt_pres_num:

                # get q counter
                pole_cnt = pfd[i]['pole_cnt']

                # get 2theta
                tt = out.params['p{}_center'.format(i)].value

                # get projection (q)
                q = rotate_project_q(tt/2, o, c, p, 270 - eta) #was 360 - p

                # store it
                pfd[i]['data'][pole_cnt,0] = q[0]
                pfd[i]['data'][pole_cnt,1] = q[1]
                pfd[i]['data'][pole_cnt,2] = q[2]     

                # tell me it's fit
                pfd[i]['fit'] = True

                # tell me what type to output
                pfd[i]['type'] = pfType
                
                # integrate
                II = np.trapz(y = comps['p{}_'.format(i)],
                              x = tt_ran,
                              dx = tt_ran[1]-tt_ran[0])

                # # store integ. int
                # pfd[i]['data'][pole_cnt,3] = II
                
                # store integ. int
                pfd[i]['data'][pole_cnt,3] = II
                pfd[i]['data'][pole_cnt,4] = II_esd['p{}_'.format(i)]

                ## counter for Δk
                pfd[i]['pole_cnt'] += 1
            
            # too fast to plot live
            if liveplot is True:

                # ## figure setup
                fig = plt.figure(figsize=(12.8,4.8),constrained_layout=True)
                gs  = fig.add_gridspec(5,4)
                ax1 = fig.add_subplot(gs[:4,:2])
                ax2 = fig.add_subplot(gs[:4,2:])
                ax3 = fig.add_subplot(gs[4,:2])

                # if k > 0:
                #     ax1.clear()
                #     ax2.clear()
                #     ax3.clear()
            
                ## print result plot      
                ax1.plot(tt_ran, inten, 'b')
                ax1.plot(tt_ran, init, 'k--', label='initial fit')
                ax1.plot(tt_ran, out.best_fit, 'r-', label='best fit')
                ax3.plot(tt_ran, out.best_fit - inten, 'g-')
                ax2.plot(tt_ran, inten, 'b')
            
                for i in tt_pres_num:
                    
                    ax2.plot(tt_ran, comps['p{}_'.format(i)], '--', label='Peak {}_{}'.format(pfd[i]['phase'],pfd[i]['ref']))
            
                # housekeeping
                ax1.legend(loc='best')
                if num_peaks < 7: ax2.legend(loc='best')
                ax1.set_ylim(0,max_int+50)
                ax2.set_ylim(0,max_int+50)
                ax1.set_ylabel('Intensity')
                ax1.set_xlabel('2θ (degrees)')
                ax2.set_ylabel('Intensity')
                ax2.set_xlabel('2θ (degrees)')
                ax3.set_ylabel('Difference')
                ax3.set_xlabel('2θ (degrees)')

                ax2.set_ylim(top=0.20*np.max(meas_int))
            
                # plt.pause(0.05)       
                # plt.show()    

                ## save fit image for posterity
                # plt.savefig(os.path.join(fitImage_path,'fit_{}'.format(label)),dpi=300)
                plt.close()

                k += 1

            t6 = time.time()

            # refine_pbar.write('plot save time:{}'.format(t6-t5))

    ## close out
    if liveplot: plt.close()
    
    # export the pole figures
    export_pfs(inst, desc, pfd, poleFig_path)

    # # write the MTEX file
    write_MTEX(desc, pfd, poleFig_path, smplSym='1', smplRot=smplRot, ranges=ranges, rot_phase=rot_phase)    

def fit_HIDRA_MP( runs, rootDir, dataDir, reduction, n_tasks=None ):

    """
    run fitting routine - using concurrent.futures
    
    can't pass Cell object thru concurrent.futures..

    """

    # https://alexwlchan.net/2019/10/adventures-with-concurrent-futures/

    # number of tasks to be scheduled at once
    if n_tasks is None:
        # play it safe...
        n_tasks = mp.cpu_count() - 2
    
    runs = iter(runs)
    
    assert reduction.lower() in ['auto', 'texture'], 'reduction must either be: \'auto\' or \'texture\''

    with cf.ProcessPoolExecutor(max_workers=n_tasks) as executor:

        futures = dict(
            ( executor.submit(fit_HIDRA, 
                              run, 
                              rootDir, 
                              dataDir, 
                              **{'reduction':reduction, 'pbaridx':pbaridx} ), run )
            for pbaridx, run in enumerate(itertools.islice(runs, n_tasks))
            )
    
        while futures:
            # Wait for the next future to complete
            done, futures = cf.wait(
                futures, return_when = cf.FIRST_COMPLETED
                )             
            
            # for fut in done:
            #     print("The outcome is {}".format(fut.result()))            
                
            # Schedule the next set of futures. We don't want more than N futures
            # in the pool at a time, to keep memory consumption down.
            
            for pbaridx, run in enumerate(itertools.islice(runs, len(done))):
                futures.add(
                    executor.submit(fit_HIDRA,
                                    run, 
                                    rootDir, 
                                    dataDir, 
                                    **{'reduction':reduction, 'pbaridx':pbaridx+n_tasks} )
                    ) 

def load_HIDRA(runNumber, rootDir, dataDir, phases, mode, exportPFs=True, pfType='jul', smplRot=None, pbaridx=None, pbarcolor='WHITE', ranges=None, rot_phase=None):
    
    """
    load a prior fit from .json files
    
    peak-fit HIDRA data -> pole figures
    uses Pseudo-Voigt function in lmfit
    
    inputs:
        desc:      sample description (str)
        rootDir:   root directory
        dataDir:   reduced data directory
        expInfo:   experiment metadata (DataFrame)
        pfd:       pole figure dictionary (setup_pfs)
        start:     start run #
        end:       end run #
        reduction: 'auto' (IP only) or 'texture' (IP & +/-5deg OOP)
        mode:      'fit' (perform fitting) or 'load' (load prior fit from .json)
    kwargs (optional):
        liveplot:     Will plot results live and save image for each fit
        sequential:   Will copy prior fit params as first guess on next fit
        exportPFs:    Will export pfs in directory (overwrite)
        pfType:       Either juelich or cartesian XYZ output format
        smplRot:      Sample rotation
        pbaridx:      Position of pbar (used in mp)
        pbarcolor:    Color of pbar (used in mp)
    
    peak-fit HIDRA data -> pole figs

    """

    # define instrument
    inst  = 'HIDRA'

    # load in .h5 file
    if mode == 'auto': 
        fname = 'HB2B_{}.h5'.format(runNumber)
        desc  = '{}_aut_uq'.format(runNumber)
    elif mode == 'texture': 
        fname = 'HB2B_{}_texture.h5'.format(runNumber)
        desc  = '{}_tex_uq'.format(runNumber)
    else: raise ValueError('reduction not recognized..')
    
    # will need to open this later
    exp_h5 = h5py.File(os.path.join(dataDir,fname), 'r')

    # read wavelength
    lmbda = exp_h5['instrument/monochromator setting/wave length'][()][0]

    # read angular data
    chi       = exp_h5['raw data/logs/chi'][()]
    phi       = exp_h5['raw data/logs/phi'][()]
    omega     = exp_h5['raw data/logs/omega'][()]
    two_theta = exp_h5['reduced diffraction data/2theta'][()]

    # number of measured patterns (for loop)
    meas_num  = len(phi)

    # read intensity data
    if mode == 'auto': #no eta slice
        
        # get from raw data/logs/2thetaSetpoint
        num_det_pos   = len(np.unique(exp_h5['raw data/logs/2thetaSetpoint'][()]))
        num_eta_slice = 1

        eta_zero = np.nan_to_num(exp_h5['reduced diffraction data/main'][()])
        max_int  = np.max(eta_zero)

    elif mode == 'texture': #should have eta slices

        num_det_pos   = len(np.unique(exp_h5['raw data/logs/2thetaSetpoint'][()]))
        num_eta_slice = 3

        eta_neg5 = np.nan_to_num(exp_h5['reduced diffraction data/eta_-5.0'][()])
        eta_zero = np.nan_to_num(exp_h5['reduced diffraction data/eta_0.0'][()])
        eta_pos5 = np.nan_to_num(exp_h5['reduced diffraction data/eta_5.0'][()])

        max_int = np.max([np.max(eta) for eta in [eta_neg5, eta_zero, eta_pos5]])

    # close the h5 file
    exp_h5.close()

    # number of measured q
    rot_num   = int((meas_num/num_det_pos)*num_eta_slice)

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
    
    ## setup pole fig dictionary ##
    pfd = {}
    for i,(d,ref,pi,tt) in enumerate(zip(d_all,ref_all,cnt_all,tt_all)):
        
        pfd[i+1] = {}
        pfd[i+1]['phase'] = name_all[pi]
        pfd[i+1]['ref']   = ''.join(map(str,ref))
        pfd[i+1]['data']  = np.zeros(( rot_num, 5 ))
        pfd[i+1]['tt']    = tt
        pfd[i+1]['lattice'] = phases[name_all[pi]].lattice
        pfd[i+1]['lattice_type'] = phases[name_all[pi]].get_type()

        # for PF Δk index
        # will iterate +1 on each insertion
        # to account for variable # of points for each PF 
        # (shouldn't be the case in CW?)
        pfd[i+1]['pole_cnt'] = 0
        
        # setup flag if it was fit or not
        pfd[i+1]['fit'] = False

    # where to store
    poleFig_path   = os.path.join(rootDir,'pole_figs',desc)
    fitResult_path = os.path.join(rootDir,'fit_results',desc,'params')
    fitImage_path  = os.path.join(rootDir,'fit_results',desc,'figures')

    if not os.path.exists(fitResult_path): os.makedirs(fitResult_path)
    if not os.path.exists(fitImage_path): os.makedirs(fitImage_path)
    if not os.path.exists(poleFig_path): os.makedirs(poleFig_path)

    # progress bar setup
    if pbaridx is None:
        refine_pbar = tqdm(range(meas_num),desc=desc)
    else:
        refine_pbar = tqdm(range(meas_num),desc=desc, position=pbaridx)
        
    border = "="*80
    clear_border = _term_move_up() + "\r" + " "*len(border) + "\r"

    # ## figure setup
    # if liveplot is True:
        
    #     fig = plt.figure(figsize=(12.8,4.8),constrained_layout=True)
    #     gs  = fig.add_gridspec(5,4)
    #     ax1 = fig.add_subplot(gs[:4,:2])
    #     ax2 = fig.add_subplot(gs[:4,2:])
    #     ax3 = fig.add_subplot(gs[4,:2])
    #     plt.pause(0.05)

    # k = 0

    ## loop over rotations
    for ri in refine_pbar:

        # easy to reference these later        
        o = omega[ri]
        c = 90 - chi[ri]
        p = 360 - phi[ri]
        
        if mode == 'auto': inner_iter = zip([eta_zero],[0])
        elif mode == 'texture': inner_iter = zip([eta_neg5, eta_zero, eta_pos5],[-5, 0, 5])
        # inner_iter = zip([eta_neg5, eta_zero, eta_pos5],[-5, 0, 5])

        # loop over data
        for meas_int,eta in inner_iter:

            counter = 0

            label = 'tt{}_o{}_c{}_p{}_e{}'.format(round(o*2),round(o),round(c),round(p),270 - eta)
        
            # get mask on invalid data on edges
            valid_data = ma.masked_where(meas_int[ri,:]==0,meas_int[ri,:])
            valid      = ~valid_data.mask

            # get 2theta range of measurement
            tt_ran = two_theta[ri,valid]

            # get weights
            weights = 1 / np.sqrt(meas_int[ri,valid])
            
            # find what peaks are present
            tt_mask  = (tt_all >= min(tt_ran)) * (tt_all <= max(tt_ran))
            tt_pres = list(itertools.compress(tt_all,tt_mask))
            # only these are present
            tt_pres_num = list(itertools.compress(range(len(tt_all)),tt_mask))
            # adjust index
            tt_pres_num = [v+1 for v in tt_pres_num]
            
            # num of peaks
            num_peaks = len(tt_pres_num)

            # setup lmfit model
            model = ConstantModel()
            for i in tt_pres_num:
                # add individual peaks
                model = model + PseudoVoigtModel(prefix='p{}_'.format(i))
                
            ## initialize params
            params = model.make_params()

            # setup file to save parameters (.json)
            fitResult = os.path.join(fitResult_path,'fitParams_{}.json'.format(label))

            with open(fitResult,'r') as f_in:
                try:
                    params = params.load(f_in)
                    out    = ModelResult(model,params)
                    comps  = out.eval_components(x=tt_ran)
                except:
                    pass 

            # # calculate weighted R (fit quality)
            # rwp = np.sum( weights * out.residual**2 ) / np.sum( weights * meas_int[ri,valid]**2 )

            # # write to console
            # # this goes fast.. only print if there's a problem
            # if not out.success: 
            #     refine_pbar.write(clear_border + '--- ω:{} | χ:{} | φ:{} | η:{} ---'.format(int(o),int(c),int(p),int(eta)))
            #     refine_pbar.update()
            #     refine_pbar.write(clear_border + 'Fit was not successful!')
            #     refine_pbar.update()
            #     refine_pbar.write(clear_border + 'Rwp : {:3.2f}%'.format(rwp*100))
            #     refine_pbar.update()
            #     refine_pbar.write(border)
            #     refine_pbar.update()

            # store peak intensity
            for i in tt_pres_num:

                # get q counter
                pole_cnt = pfd[i]['pole_cnt']

                # get 2theta
                # tt = out.params['p{}_center'.format(i)].value
                tt  = pfd[i]['tt']
                ref = pfd[i]['ref']
                
                # print(out.params['p{}_center'.format(i)])
                # print(f'{tt} - {ref}')
                # print(f'{tt_ran[:10]}')
                # print('--- ω:{} | χ:{} | φ:{} | η:{} ---'.format(int(o),int(c),int(p),int(eta)))
                # raise Exception

                # get projection (q)
                q = rotate_project_q(tt/2, o, c, 360 + p, 270 - eta) #was 360 - p

                # pfd[i]['offset'][pole_cnt,0] = (tt/2) - o

                # store it
                pfd[i]['data'][pole_cnt,0] = q[0]
                pfd[i]['data'][pole_cnt,1] = q[1]
                pfd[i]['data'][pole_cnt,2] = q[2]     

                # tell me it's fit
                pfd[i]['fit'] = True

                # tell me what type to output
                pfd[i]['type'] = pfType
                
                # integrate
                II = np.trapz(y = comps['p{}_'.format(i)],
                              x = tt_ran,
                              dx = tt_ran[1]-tt_ran[0])

                # store integ. int
                pfd[i]['data'][pole_cnt,3] = II
                pfd[i]['data'][pole_cnt,4] = 0.0
                
                ## counter for Δk
                pfd[i]['pole_cnt'] += 1
    
    # export the pole figures
    if exportPFs: export_pfs(inst, desc, pfd, poleFig_path)

    # # # write the MTEX file
    # write_MTEX(desc, pfd, poleFig_path, smplSym='1', smplRot=smplRot, ranges=ranges, rot_phase=rot_phase)    
   
    return pfd 

rootDir = '../Data/HIDRA/' #
tex_dataDir = os.path.join(rootDir,'texturereduce') #texreduce

phases = {'FCC':Cell.cubic(3.5938,lattice_type='F'),
          'BCC':Cell.cubic(2.8728,lattice_type='I')}

smplRot = {'1600':[('z',97), ('y',90), ('x',0)], 
           '1601':[('z',93), ('y',90), ('x',0)], 
           '1602':[('z',90), ('y',90), ('x',0)]} 

runs   = [1600,1601,1602]

for run in runs:

    if run in [1600]:
        rot_phase = 'FCC'
    else:
        rot_phase = 'BCC'

    # texture reduced
    pfd = load_HIDRA(run,rootDir,tex_dataDir,phases,reduction='texture',
                        exportPFs=False,smplRot=smplRot[str(run)],pfType='jul',
                        ranges={'BCC':'0:0.5:5.5','FCC':'0:0.5:4'},rot_phase=rot_phase)

