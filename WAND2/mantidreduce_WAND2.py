
# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def reduce_texture(eta, nexusFile, run_number):
    eta_step = 5
    eta_max = 8
    eta_min = -8
    reduced_data = []
    for eta_roi in [-5, 0, 5]:
        MaskSpectra(InputWorkspace=nexusFile,
                    InputWorkspaceIndexType='WorkspaceIndex',
                    InputWorkspaceIndexSet=np.where(eta > (eta_roi + eta_step / 2))[0],
                    OutputWorkspace='test.nxs')
        MaskSpectra(InputWorkspace='test.nxs',
                    InputWorkspaceIndexType='WorkspaceIndex',
                    InputWorkspaceIndexSet=np.where(eta < (eta_roi - eta_step / 2))[0],
                    OutputWorkspace='test.nxs')
 
        reduced_data.append(WANDPowderReduction(InputWorkspace='test.nxs',
                            Target='Theta',
                            NumberBins=2000,
                            OutputWorkspace='{}_eta_{}'.format(run_number, eta_roi)))
                        
    return reduced_data

# Fill this in with a specific IPTS
IPTS = 11111

runNumbers = range(1, 2 + 1)

# This may also need to be updated to a more current measurement
VanFile = '/HFIR/HB2C/IPTS-23858/nexus/HB2C_376206.nxs.h5'
nexusFile = 'HB2C_{}.nxs'.format(376206)
LoadWAND(Filename='/HFIR/HB2C/IPTS-23858/nexus/HB2C_376206.nxs.h5', OutputWorkspace=nexusFile)

instrument = mtd[nexusFile].getInstrument()
detIDs = mtd[nexusFile].getNumberHistograms()
x, y, z = np.zeros(detIDs), np.zeros(detIDs), np.zeros(detIDs)
for i in range(detIDs):
    x[i], y[i], z[i] = instrument.getDetector(i).getPos()
eta = np.arctan2(y, x) * 180. / np.pi

van_data = reduce_texture(eta, nexusFile, 376206)

for runNumber in runNumbers:
    nexusFile = 'HB2C_{}.nxs'.format(runNumber)
    LoadWAND(Filename='/HFIR/HB2C/IPTS-22335/nexus/HB2C_{}.nxs.h5'.format(runNumber), OutputWorkspace=nexusFile)

    texture_data = reduce_texture(eta, nexusFile, runNumber)

    Divide(LHSWorkspace=texture_data[0], RHSWorkspace=van_data[0], OutputWorkspace='normalized', AllowDifferentNumberSpectra=True)
    outfile = '/HFIR/HB2C/IPTS-{}/shared/norm_texture/HB2C_{}_eta_{}.xye'.format(IPTS, runNumber, -5)
    SaveFocusedXYE('normalized', outfile, Format='MAUD')
    Divide(LHSWorkspace=texture_data[1], RHSWorkspace=van_data[1], OutputWorkspace='normalized', AllowDifferentNumberSpectra=True)
    outfile = '/HFIR/HB2C/IPTS-{}/shared/norm_texture/HB2C_{}_eta_{}.xye'.format(IPTS, runNumber, 0)
    SaveFocusedXYE('normalized', outfile, Format='MAUD')
    Divide(LHSWorkspace=texture_data[2], RHSWorkspace=van_data[2], OutputWorkspace='normalized', AllowDifferentNumberSpectra=True)
    outfile = '/HFIR/HB2C/IPTS-{}/shared/norm_texture/HB2C_{}_eta_{}.xye'.format(IPTS, runNumber, 5)
    SaveFocusedXYE('normalized', outfile, Format='MAUD')
    
    DeleteWorkspace(nexusFile)
    DeleteWorkspace(texture_data[0])
    DeleteWorkspace(texture_data[1])
    DeleteWorkspace(texture_data[2])
