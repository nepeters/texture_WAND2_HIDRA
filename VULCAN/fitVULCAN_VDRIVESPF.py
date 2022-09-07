# -*- coding: utf-8 -*-
# @Author: Nate Peterson
# @Date:   2022-02-14 15:37:52
# @Last Modified by:   Nate Peterson
# @Last Modified time: 2022-02-22 12:22:22


"""
Extract pole figures from a VDRIVESPF fit


"""

import os, sys, copy
import pandas as pd
import numpy as np
from datetime import datetime


# location of VDRIVESPF output file
spf_file = r"VDRIVE\VDriveSPF-172644-172763-bk1.txt"
spf_dir, spf_fname  = os.path.split(spf_file)

# location of output pole figs
pf_dir   = os.path.join(spf_dir,'pole_figs')
if not os.path.exists(pf_dir): os.makedirs(pf_dir)

# sample description
desc  = 'Steel_6'
phase = 'BCC' 
inst  = 'VULCAN'

# import the fit results
spf_results = pd.read_csv(spf_file,sep='\t')

# store angles
meas_O = spf_results['O'] - spf_results['O'].iloc[0] + 45
# meas_O = spf_results['O']
meas_HROT = spf_results['HROT'] - spf_results['HROT'].iloc[0]
# meas_HROT = spf_results['HROT']
meas_VROT = spf_results['VROT'] #not used..

## get poles
# all of them
poles = [pol for pol in spf_results.columns if pol.startswith('I_V_')]
# unique
names = []
banks = []
for pol in poles:
    # name_bank
    _, _, name, bank = pol.split('_')
    if name in names: pass
    else:
        names.append(name)
    if bank in banks: pass
    else:
        banks.append(bank)
poles = copy.deepcopy(names)
# convert to int for index
banks = [int(bank) for bank in banks]

## focus_ew = 0
# 6 banks
if 6 in banks:
    all_tt   = [-90.0, -90.0, -90.0, 90.0, 90.0, 90.0]
    all_eta  = [-11.6, 0.0, 11.6, -11.6, 0.0, 11.6]  
elif 7 in banks:
    raise ValueError
# 2 banks
else:
    all_tt  = [-90.0, 90.0]
    all_eta = [0.0, 0.0]   

for pol in poles:

    # data array
    data = np.zeros((len(banks)*len(meas_O), 4))
    di   = 0

    for bid, tt, eta in zip(banks, all_tt, all_eta):
        # measured intensity
        meas_int = spf_results['I_V_{}_{}'.format(pol, bid)] 
        # loopw
        for int, omega, phi in zip(meas_int, meas_O, meas_HROT):
            # get q
            q = rotate_project_q(tt/2, omega, 0, phi, eta)
            data[di,:3] = q
            # data[di,:2] = to_jul(*q)
            data[di,3]  = int
            
            di += 1

    now = datetime.now()
    header  = 'x, y, z, integ. int\n'   
    # header  = 'alpha beta intensity\n'        
    pf_name = '{}_{}_{}{}'.format(desc,phase,pol,'.pf')   

    with open(os.path.join(pf_dir,pf_name),'w') as pf_out:
        pf_out.write('{}\n'.format(inst))
        pf_out.write('Single peak fit using VDRIVESPF, extracted at {}\n'.format(now))
        pf_out.write(header)
        np.savetxt(pf_out,data)