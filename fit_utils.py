# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 2021

@author: Nate Peterson
"""

"""
Utility routines for WAND and HIDRA pole figure reduction / fit

"""

import numpy
from datetime import datetime
import os
from lmfit.model import ModelResult
import sys
from datetime import datetime

## Coordinate transforms ##

def to_jul(h1, h2, h3):

    """
    converts h1,h2,h3 to alpha,beta 
    using convention in Chapter 8.3 in Bob He Two-Dimensional X-ray Diffraction
    """

    h_length = numpy.sqrt(numpy.square(h1) + numpy.square(h2))

    # alpha =  numpy.rad2deg(numpy.arcsin(numpy.abs(h3)))
    alpha = numpy.rad2deg(numpy.arccos(h_length))
    beta = numpy.rad2deg(numpy.arccos(h1 / h_length))

    if type(h2) == numpy.float64:
        if h2 < 0:
            beta = beta * -1
    else:
        beta[h2<0] = beta[h2<0] * -1

    return (90 - alpha, beta)

def rotate_project_q(theta: float, omega: float, chi: float, phi: float, eta: float):
    """
    Projection of angular dependent data onto pole sphere. Analytical solution taken from
    Chapter 8.3 in Bob He Two-Dimensional X-ray Diffraction

    _______________________
    :param two_theta:
    :param omega:
    :param chi:
    :param phi:
    :return: 3-tuple as the projection (h1, h2, h3)
    """

    sp = numpy.sin(numpy.deg2rad(phi))
    sw = numpy.sin(numpy.deg2rad(omega))
    sc = numpy.sin(numpy.deg2rad(chi))
    sg = numpy.sin(numpy.deg2rad(eta))
    st = numpy.sin(numpy.deg2rad(theta))

    cp = numpy.cos(numpy.deg2rad(phi))
    cw = numpy.cos(numpy.deg2rad(omega))
    cc = numpy.cos(numpy.deg2rad(chi))
    cg = numpy.cos(numpy.deg2rad(eta))
    ct = numpy.cos(numpy.deg2rad(theta))

    h1 = st*(sp*sc*sw+cp*cw) + ct*cg*sp*cc - ct*sg*(sp*sc*cw-cp*sw)
    h2 = -st*(cp*sc*sw-sp*cw) - ct*cg*cp*cc + ct*sg*(cp*sc*cw+sp*sw)
    h3 = st*cc*sw - ct*sg*cc*cw - ct*cg*sc

    return h1, h2, h3

## Import / export ##

def load_params(params, model, x, fitResult):

    """
    load in lmfit parameters from .json file
    evaluate model over x
    calculate components (peaks)

    """

    # load old parameters from json files
    with open(fitResult,'r') as f_in:
        try:
            params = params.load(f_in)
            out    = ModelResult(model,params)
            comps  = out.eval_components(x=x)
        except:
            raise FileNotFoundError('param file wouldn\'t load.. check if it exists') 

    return params,out,comps

def export_pfs(inst, desc, pfd, pf_path):

    """
    write data to pole figures readable by MTEX

    inputs:
        desc: sample description
        pfd: pole figure data (dict)
        pf_path: save path (str)

    """
            
    now = datetime.now()
    
    # loop over PFs
    for pi, pfd in pfd.items():
        
        if pfd['fit']:

            if pfd['type'] == 'jul': 

                # if pfd['data'].shape[1] == 4:
                
                #     ab = to_jul(pfd['data'][:,0],pfd['data'][:,1],pfd['data'][:,2])

                #     exten   = '.jul'
                #     header  = 'alpha beta intensity\n'
                #     dat_out = numpy.array((ab[0],ab[1],pfd['data'][:,-1])).T

                # else:

                exten   = '.jul'
                header  = 'alpha beta intensity\n'
                dat_out = pfd['data']    

            elif pfd['type'] == 'xyz': 
                
                exten   = '.pf'
                header  = 'x, y, z, integ. int\n'
                dat_out = pfd['data']

            pf_name = '{}_{}_{}{}'.format(desc,
                                          pfd['phase'],
                                          pfd['ref'],
                                          exten)   
    
            with open(os.path.join(pf_path,pf_name),'w') as pf_out:
                pf_out.write('{}\n'.format(inst))
                pf_out.write(header)
                pf_out.write('Single peak fit using lmfit routine (P-V profile) at {}\n'.format(now))
                numpy.savetxt(pf_out,dat_out, fmt='%3.5f')
        else:
            # print('skipping {}_{}.. not fit'.format(pfd['phase'],pfd['ref']))
            pass

def write_MTEX(desc, pfd, pf_path, smplSym='1', smplRot=None, ranges=None, rot_phase=None):

    """
    write an MTEX MATLAB script file for import

    """

    ## we will print the date for records
    now = datetime.now()

    # get unique phases
    phases = list(set([pfd['phase'] for pi,pfd in pfd.items()]))
    # will store the pole figure type; either jul or xyz for the phase
    # should be consistent
    fmt_type = {}

    # reorder phases if rot_phase is defined
    # this is the phase to run the centerSpecimen algo on
    if rot_phase is None: pass
    else:
        if rot_phase in phases:
            phases.pop(phases.index(rot_phase))
            phases.insert(0,rot_phase)
    
    if ranges is None:
        ranges = {}
        for phase in phases:
            ranges[phase] = '0:0.5:6'

    mtexFile = 'analyzePFs_{}.m'.format(desc)

    with open(os.path.join(pf_path,mtexFile),'w') as f_out:
        f_out.write('%% Import Script for PoleFigure Data\n\n')
        f_out.write('% created from reduce_fit python script @ {}\n'.format(now))
        f_out.write('% for MTEX version 5.30\n\n')
        
        f_out.write('%% Specify Crystal and Specimen Symmetries\n\n')
        f_out.write('% crystal symmetry\n')
        
        for pi,phase in enumerate(phases):
            i = 1
            printed = False
            while printed is False:
                # okay we have the phase
                if pfd[i]['phase'] == phase and pfd[i]['fit'] is True:
                    #get laue group
                    if pfd[i]['lattice'] == 'cubic': # only should see this..
                        f_out.write('CS_{} = crystalSymmetry(\'m-3m\', [1 1 1]); % lattice-type:{}\n'.format(pi+1,pfd[i]['lattice_type']))
                    else:
                        f_out.write('CS_{} = crystalSymmetry(\'1\');\n'.format(pi+1))
                    # write specimenSymmetry (none)
                    f_out.write('SS_{} = specimenSymmetry(\'1\');\n'.format(pi+1))
                    printed = True
                else:
                    i += 1

        f_out.write('\n% plotting convention\n')
        f_out.write('setMTEXpref(\'xAxisDirection\',\'north\');\n')
        f_out.write('setMTEXpref(\'zAxisDirection\',\'outOfPlane\');\n\n')

        if smplRot is None:
            f_out.write('rot1 = rotation(\'axis\',zvector,\'angle\',0*degree);\n')
            f_out.write('rot2 = rotation(\'axis\',yvector,\'angle\',0*degree);\n')
            f_out.write('rot3 = rotation(\'axis\',xvector,\'angle\',0*degree);\n')
        else:
            rots = []
            for ri, rot in enumerate(smplRot):
                f_out.write('rot{} = rotation(\'axis\',{}vector,\'angle\',{}*degree);\n'.format(ri+1, rot[0].lower(), rot[1]))
                rots.append('rot{}'.format(ri+1))

        rots.reverse()
        f_out.write('rot = {};\n\n'.format(' * '.join(rots)))

        f_out.write('%% Specify File Names\n\n')
        f_out.write('% path to files\n')

        # write out file names for import
        for pi,phase in enumerate(phases):
            # we can loop over all.. need to write all
            f_out.write('fname_{} = {{...\n'.format(pi+1))
            for i, pf in pfd.items():
                
                # okay we have the phase
                if pf['phase'] == phase and pf['fit'] is True:
                    
                    # set file extension
                    if pf['type'] == 'xyz': exten = '.pf'
                    elif pf['type'] == 'jul': exten = '.jul'
                    else: raise TypeError('pf type: {} not recognized..'.format(str(pf['type'])))
                    fmt_type[pi] = pf['type']         
                    
                    # write a single line for each pole figure
                    f_out.write('  \'.\{}_{}_{}{}\',...\n'.format(desc,
                                                                  pf['phase'],
                                                                  pf['ref'],
                                                                  exten))
            f_out.write('  };\n\n')

        f_out.write('\n%% Specify Miller Indice\n\n')
        # write out miller indicies
        for pi,phase in enumerate(phases):
            # we can loop over all.. need to write all
            f_out.write('h_{} = {{...\n'.format(pi+1))
            for i, pf in pfd.items():
                # okay we have the phase
                if pf['phase'] == phase and pf['fit'] is True:
                    # write a single line for each pole figure
                    f_out.write('  Miller({},{},{},CS_{}),...\n'.format(pf['ref'][0],
                                                                        pf['ref'][1],
                                                                        pf['ref'][2],
                                                                        pi+1))
            f_out.write('  };\n\n')  

        f_out.write('%% Import the Data\n\n')
        f_out.write('% create a Pole Figure variable containing the data\n\n')  
        # write out miller indicies
        for pi,phase in enumerate(phases):
            # we can loop over all.. need to write all
            f_out.write('% _{} is phase: {}\n\n'.format(pi+1, phase))
            
            if fmt_type[pi] == 'xyz':
                
                f_out.write('pf_{0} = PoleFigure.load(fname_{0},h_{0},CS_{0},SS_{0},\'interface\',\'generic\',...\n'.format(pi+1))
                f_out.write('  \'ColumnNames\', { \'x\' \'y\' \'z\' \'Intensity\'});\n')
            
            elif fmt_type[pi] == 'jul':

                f_out.write('pf_{0} = PoleFigure.load(fname_{0},h_{0},CS_{0},SS_{0},\'interface\',\'juelich\');\n'.format(pi+1))
            
            else: raise TypeError('pf type not recognized..')

            f_out.write('odf_{0} = calcODF(pf_{0});\n'.format(pi+1))
            if smplRot is None: pass
            else: f_out.write('odf_{0} = rotate(odf_{0},rot);\n'.format(pi+1))
            f_out.write('save(\'odf_{0}_{1}\',\'odf_{2}\');\n\n'.format(phase, desc, pi+1))

            if smplRot is None: pass
            else:
                f_out.write('% write out rotated pole figures\n')
                f_out.write('rotPF_{0} = rotate(pf_{0},rot);\n'.format(pi+1))
                f_out.write('export(rotPF_{},\'{}\',\'degree\')\n\n'.format(pi+1,phase))

            if rot_phase is None:
                f_out.write('[odf_{0},cen_rot_{0},~,~] = centerSpecimen(odf_{0},xvector,\'Fourier\');\n\n'.format(pi+1))
            else:
                if rot_phase in phases:
                    if phase == rot_phase:
                        f_out.write('[odf_{0},cen_rot_{0},~,~] = centerSpecimen(odf_{0},xvector,\'Fourier\');\n\n'.format(pi+1))
                    else:
                        f_out.write('odf_{0} = rotate(odf_{0},cen_rot_{1});\n\n'.format(pi+1, phases.index(rot_phase)+1))
                else: raise Exception('rot_phase; {} not in phase list; {}'.format(rot_phase,phases))

            f_out.write('figure\n')
            f_out.write('plotSection(odf_{},\'phi2\',45*degree)\n'.format(pi+1))
            f_out.write('xticks([0 45 90 135 180 225 270 315 360])\n')
            if phase == 'FCC':
                f_out.write('setColorRange([0 5],\'current\')\n')
            elif phase == 'BCC':
                f_out.write('setColorRange([0 14],\'current\')\n')
            f_out.write('mtexColorbar(\'title\',\'{}\')\n'.format(phase))
            f_out.write('saveFigure(\'{}_phi2_45.jpg\')\n'.format(phase))
            f_out.write('close\n\n')

        f_out.write('%% Plot recalculated pfs\n\n')
        # write out miller indicies
        for pi,phase in enumerate(phases):
            # we can loop over all.. need to write all
            f_out.write('% _{} is phase: {}\n\n'.format(pi+1, phase))
            
            f_out.write('nh_{0} = length(h_{0});\n'.format(pi+1))
            f_out.write('repf_{0} = calcPoleFigure(odf_{0},h_{0});\n'.format(pi+1))
            f_out.write('figure\n')
            # f_out.write('plot(repf_{0}{{nh_{0}-1:nh_{0}}},\'contourf\',\'colorrange\',\'equal\');\n'.format(pi+1))
            f_out.write('plot(repf_{0}{{nh_{0}-1:nh_{0}}},\'contourf\',{1});\n'.format(pi+1,ranges[phase]))
            f_out.write('mtexColorbar(\'title\',\'{}\')\n'.format(phase)) 
            f_out.write('saveFigure(\'{}_pfs.jpg\')\n'.format(phase))
            f_out.write('close\n\n') 

            f_out.write('figure\n')
            f_out.write('plot(rotPF_{0}{{nh_{0}}})\n'.format(pi+1))
            f_out.write('saveFigure(\'{0}_rawpfs.jpg\')\n'.format(phase))            
            f_out.write('close\n\n') 

        # volume fractions of components
        f_out.write('%% Volume fractions of texture components\n\n')
        f_out.write('BFOri = [[35.3 45.0 0.0 ];...\n')
        f_out.write('           [33.6 47.7 5.0 ];...\n')
        f_out.write('            [32.1 51.0 10.0];...\n')
        f_out.write('            [31.1 54.7 15.0];...\n')
        f_out.write('            [31.3 59.1 20.0];...\n')
        f_out.write('            [35.2 64.2 25.0];...\n')
        f_out.write('            [46.1 69.9 30.0];...\n')
        f_out.write('            [49.5 76.2 35.0];...\n')
        f_out.write('            [51.8 83.0 40.0];...\n')
        f_out.write('            [54.7 90.0 45.0];...\n')
        f_out.write('            [90.0 35.3 45.0];...\n')
        f_out.write('            [80.2 35.4 50.0];...\n')
        f_out.write('            [73.0 35.7 55.0];...\n')
        f_out.write('            [66.9 36.2 60.0];...\n')
        f_out.write('            [61.2 37.0 65.0];...\n')
        f_out.write('            [55.9 38.0 70.0];...\n')
        f_out.write('            [50.7 39.2 75.0];...\n')
        f_out.write('            [45.6 40.8 80.0];...\n')
        f_out.write('            [40.5 42.7 85.0];...\n')
        f_out.write('            [35.3 45.0 90.0]]*degree;\n\n')

        f_out.write('beta = orientation(\'Euler\', BFOri(:,1), BFOri(:,2), BFOri(:,3),crystalSymmetry(\'m-3m\'),specimenSymmetry(\'1\'),\'bunge\');\n')        
        f_out.write('alpha = fibre.alpha(crystalSymmetry(\'m-3m\'),specimenSymmetry(\'1\'),\'full\');\n')
        f_out.write('gamma = fibre.gamma(crystalSymmetry(\'m-3m\'),specimenSymmetry(\'1\'),\'full\');\n\n')

        for pi, phase in enumerate(phases):
            
            if phase == 'BCC':
                f_out.write('vol_alp = 100 * volume(odf_{0},alpha,10*degree);\n'.format(pi+1))
                f_out.write('vol_gam = 100 * volume(odf_{0},gamma,10*degree);\n'.format(pi+1))
                f_out.write('fprintf(\' BCC volume fractions | alpha fibre: %f6, gamma fibre: %f6\\n\',vol_alp,vol_gam)\n\n')
            elif phase == 'FCC':
                f_out.write('vol_bet = 100 * volume(odf_{0},beta,10*degree);\n'.format(pi+1))
                f_out.write('fprintf(\' FCC volume fractions | beta fibre: %f6\\n\',vol_bet)\n\n')

        # texture index
        f_out.write('%% texture strength (index) calculation\n\n')
        index_key = '| {0}: %f6 '
        indexes     = []
        index_prnts = []
        for pi, phase in enumerate(phases):
            f_out.write('ti_{0} = textureindex(odf_{0});\n'.format(pi+1))
            indexes.append('ti_{0}'.format(pi+1))
            index_prnts.append(index_key.format(phase))
        f_out.write('\n')
        f_out.write('fprintf(\' Texture index {0}'.format(''.join(index_prnts)))
        f_out.write('\\n\',{});'.format(','.join(indexes)))
        f_out.write('\n\n')
        
        # elasticity check
        f_out.write('%% betaBrass elasticity homogenization \n')
        f_out.write('C11 = 52;\n')
        f_out.write('C12 = 34;\n')
        f_out.write('C13 = C12;\n')
        f_out.write('C23 = C12;\n')
        f_out.write('C22 = 52;\n')
        f_out.write('C33 = C22;\n')
        f_out.write('C44 = 173;\n')
        f_out.write('C55 = C44;\n')
        f_out.write('C66 = C44;\n')
        f_out.write('T_cs = crystalSymmetry(\'m-3m\', [1 1 1], \'mineral\',...\n')
        f_out.write('    \'Copper\', \'color\', \'light blue\');\n\n')

        f_out.write('T_ijkl = stiffnessTensor(...\n')
        f_out.write('    [[  C11     C12    C13    0.0     0.0    0.0];...\n')
        f_out.write('    [   C12     C22    C23    0.0     0.0    0.0];...\n')
        f_out.write('    [   C13     C23    C33    0.0     0.0    0.0];...\n')
        f_out.write('    [   0.0     0.0    0.0    C44     0.0    0.0];...\n')
        f_out.write('    [   0.0     0.0    0.0    0.0     C55    0.0];...\n')
        f_out.write('    [   0.0     0.0    0.0    0.0     0.0    C66]],T_cs);\n\n')

        # write out miller indicies
        for pi,phase in enumerate(phases):

            f_out.write('fodf_{0} = FourierODF(odf_{0},4);\n'.format(pi+1))
            f_out.write('[TVoigt_{0}, TReuss_{0}, THill_{0}] = calcTensor(fodf_{0},T_ijkl);\n'.format(pi+1))    
            f_out.write('E_{0}  = THill_{0}.YoungsModulus;\n\n'.format(pi+1)) 

            f_out.write('figure\n')
            f_out.write('plot(E_{0},\'contourf\',\'complete\',\'upper\',\'colorrange\',[115 145]);\n'.format(pi+1))
            f_out.write('mtexColorbar\n')
            f_out.write('saveFigure(\'{}_betabrassE.jpg\')\n'.format(phase))
            f_out.write('close\n\n') 

            f_out.write('figure\n')
            f_out.write('plotPDF(fodf_{0},h_{0}{{end}},\'contourf\',\'colorrange\',[0 2])\n'.format(pi+1))
            f_out.write('mtexColorbar\n')
            f_out.write('saveFigure(\'{}_fodf4_pfs.jpg\')\n'.format(phase))
            f_out.write('close\n\n') 
