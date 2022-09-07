# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:21:37 2019

@author: npk
"""

"""
Create MAUD file from WAND2 xye files

group by chi tilt?
"""

import os,time
from datetime import datetime
import numpy as np
import pandas as pd
from tqdm import tqdm
from tqdm.utils import _term_move_up

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

def writeMAUD_WAND2(desc, rootDir, dataDir, expInfo, runs, rune, caglioti, mode, wimv_cell_size=5, smplRot=None, fullPar=False, bkgd=None):
    
    """
    write out a .cif file that can be read in by MAUD
    
    group data by tilts?

    """

    now     = datetime.now() 
    tz      = time.tzname[0]
    timestr = "%a %b %d %H:%M:%S {} %Y".format(tz)

    # for tqdm
    border = "="*80
    clear_border = _term_move_up() + "\r" + " "*len(border) + "\r"    

    assert wimv_cell_size in [5,10,15], 'wimv_cell_size must either be 5, 10, 15 (deg)'

    # set path to maud files
    maudFile_path = os.path.join(rootDir,'maud','samples_texred_{}deg'.format(wimv_cell_size),desc)
    if not os.path.exists(maudFile_path): os.makedirs(maudFile_path)

    # init dict
    diffData = {}

    if smplRot is None:
        smplRot = {'omega':0.0,
                      'chi':0.0,
                      'phi':0.0}
    elif isinstance(smplRot,dict): pass
    else: raise ValueError('smplRot must be dict')

    # we will group by chi tilt
    chiTilts = pd.unique(expInfo['Chi'])
    if fullPar is True: maudFile = os.path.join(maudFile_path,'HB2C_MAUD_{}_out.par'.format(desc)) 
    else: maudFile = os.path.join(maudFile_path,'HB2C_MAUD_{}_out'.format(desc))

    # print sample name to console
    print('\n --- {} ---\n'.format(desc))

    # progress bar setup
    refine_pbar = tqdm(np.linspace(runs,rune,num=(rune-runs)+1,dtype=int),
                       desc=desc)

    # get diffData
    for run in refine_pbar:

        if mode == 'auto':
            ## no in/out of plane reduction
            dataFileName = 'HB2C_{}.xye'.format(run)
            dataFile = os.path.join(dataDir,dataFileName)            
            diffData[run] = xye_reader(dataFile)

        elif mode == 'texture':
            ## in/out of plane reduced
            # the run# will have _ChiAng appended
            diffData[run] = {}
            for eta in [-5,0,5]:
                dataFileName = 'HB2C_{}_eta_{}-0.xye'.format(run,eta)
                dataFile     = os.path.join(dataDir,dataFileName)
                diffData[run][eta] = xye_reader(dataFile)

    # write maud file
    with open(maudFile,'w') as fout:

        if fullPar is True:
            fout.write('data_global\n')
            fout.write('_audit_creation_date \'{}\'\n'.format(now.strftime(timestr)))
            fout.write('_audit_creation_method \'Maud, version 2.80\'\n')
            fout.write('_audit_update_record \'{}\'\n'.format(now.strftime(timestr)))
            fout.write('_computing_structure_refinement \'Maud, version 2.80\'\n')
            fout.write('_refine_ls_R_factor_all 0.0\n')
            fout.write('_refine_ls_wR_factor_all 0.0\n')
            fout.write('_refine_ls_goodness_of_fit_all 0.0\n')
            fout.write('_publ_contact_author_name \'Luca Lutterotti\'\n')
            fout.write('_publ_section_title \'Put a title here\'\n')
            fout.write('_pd_proc_ls_extract_int \'end of iteration\'\n')
            fout.write('_pd_proc_ls_texture_comp \'end of iteration\'\n')
            fout.write('_computing_reduce_memory_occ true\n')
            fout.write('_pd_proc_ls_theoretical_weight false\n')
            fout.write('_pd_proc_ls_extract_pos \'end of iteration\'\n')
            fout.write('_pd_proc_ls_strain_comp \'end of iteration\'\n')
            fout.write('_pd_proc_ls_extract_Fhkl \'end of iteration\'\n')
            fout.write('_pd_proc_ls_Fhkl_comp \'end of iteration\'\n')
            fout.write('_pd_proc_ls_weight_scheme sqrt\n')
            fout.write('_refine_ls_weighting_scheme WgtSS\n')
            fout.write('_refine_ls_WSS_factor 0.0\n')
            fout.write('_maud_store_spectra_with_analysis true\n')
            fout.write('_riet_remove_phases_under 0.001\n')
            fout.write('_riet_refine_cell_over 0.1\n')
            fout.write('_riet_refine_sizestrain_over 0.1\n')
            fout.write('_riet_refine_crystal_structure_over 0.1\n')
            fout.write('_riet_refine_texture_over 0.15\n')
            fout.write('_riet_refine_strain_over 0.25\n')
            fout.write('_pd_proc_ls_interpolation_comp \'end of iteration\'\n')
            fout.write('_maud_analysis_name HB2C_MAUD_{}\n'.format(desc))
            fout.write('_maud_analysis_prog_number 0\n')
            fout.write('_maud_store_structure_factors_with_analysis true\n')
            fout.write('_maud_store_texture_factors_with_analysis true\n\n')
            fout.write('#subordinateObject_Marqardt Least Squares\n\n')
            fout.write('_computing_refinement_algorithm \'Marqardt Least Squares\'\n\n')
            fout.write('_refine_ls_number_iteration 5\n')
            fout.write('_riet_refine_ls_precision 0.00000001\n')
            fout.write('_riet_refine_ls_derivative_step 0.0001\n')
            fout.write('_riet_refine_ls_lambda 0.01\n')
            fout.write('_riet_refine_ls_double_derivative false\n\n')
            fout.write('#end_subordinateObject_Marqardt Least Squares\n\n\n\n')
            fout.write('data_Sample_Sample_x\n')
            fout.write('_pd_spec_description \'Sample description\'\n')
            fout.write('_riet_thin_film_phase_refinement films\n')
            fout.write('_reflectivity_layer_qc_from_phases false\n')
            fout.write('_pd_spec_orientation_omega {} #min 0.0 #max 0.0\n'.format(smplRot['omega']))
            fout.write('_pd_spec_orientation_chi {} #min 0.0 #max 0.0\n'.format(smplRot['chi']))
            fout.write('_pd_spec_orientation_phi {} #min 0.0 #max 0.0\n'.format(smplRot['phi']))
            fout.write('_riet_par_spec_displac_x 0 #min 0.0 #max 0.0\n')
            fout.write('_riet_par_spec_displac_y 0 #min 0.0 #max 0.0\n')
            fout.write('_riet_par_spec_displac_z 0 #min 0.0 #max 0.0\n')
            fout.write('_pd_spec_size_axial 0 #min 0.0 #max 100.0\n')
            fout.write('_pd_spec_size_equat 0 #min 0.0 #max 100.0\n')
            fout.write('_pd_spec_size_thick 0 #min 0.0 #max 100.0\n')
            fout.write('_pd_spec_size_radius 0 #min 0.0 #max 100.0\n')
            fout.write('_pd_spec_size_radius_y 0 #positive #min 0.0 #max 100.0\n\n')
            fout.write('#subordinateObject_No shape\n\n')
            fout.write('_pd_spec_shape \'No shape\'\n\n\n')
            fout.write('#end_subordinateObject_No shape\n\n\n')
            fout.write('#subordinateObject_None Layer workout\n\n')
            fout.write('_riet_layer_solution_method \'None Layer workout\'\n\n\n')
            fout.write('#end_subordinateObject_None Layer workout\n\n\n')
            fout.write('#subordinateObject_No precession\n\n')
            fout.write('_riet_spec_precession_error \'No precession\'\n\n\n')
            fout.write('#end_subordinateObject_No precession\n\n\n')
            fout.write('#subordinateObject_layer1\n\n')
            fout.write('_riet_spec_layer_id \'layer1\'\n\n')
            fout.write('_reflectivity_layer_qc_from_phases false\n')
            fout.write('_riet_par_spec_layer_thickness 1.0E7 #positive #min 0.1 #max 1.0E8\n')
            fout.write('_reflectivity_layer_critical_qc 0.04 #positive #min 0.0010 #max 0.5\n')
            fout.write('_reflectivity_layer_absorption 2.0E-7 #positive #min 1.0E-8 #max 1.0E-5\n')
            fout.write('_reflectivity_layer_roughness 2.0 #positive #min 0.01 #max 100.0\n\n')
            fout.write('loop_\n')
            fout.write('_pd_phase_atom_%\n')
            fout.write('  0.5(0.0) #positive #min 0.0 #max 1.0\n')
            fout.write('  0.5(0.0) #positive #min 0.0 #max 1.0\n\n\n')
            fout.write('#end_subordinateObject_layer1\n\n\n')
        else: pass

        # loop over chi tilts
        for chi in chiTilts:

            # get all files that share the same chi tilt
            chiTiltFiles = expInfo[expInfo['Chi'] == chi] 
                
            # dataset name in .cif file
            datasetString = 'data_dataset_chi{}'.format(chi)
            
            refine_pbar.write(clear_border + 'Writing... {}'.format(datasetString))
            refine_pbar.update()

            # header for dataset
            fout.write('#subordinateObject_{}\n'.format(datasetString))
            fout.write('_pd_meas_dataset_id \'{}\'\n'.format(datasetString))
            fout.write('_pd_meas_datetime_initiated \'Date/time meas\'\n')
            fout.write('_pd_meas_info_author_name ?\n')
            fout.write('_riet_meas_datafile_format ?\n')
            fout.write('_pd_proc_ls_background_function ?\n')
            fout.write('_pd_proc_ls_profile_function ?\n')
            fout.write('_pd_proc_ls_peak_cutoff 30\n')

            ## set 2theta min & max ##
            fout.write('_pd_proc_2theta_range_min {}\n'.format(twotheta_min))
            fout.write('_pd_proc_2theta_range_max {}\n'.format(twotheta_max))

            fout.write('_pd_proc_2theta_range_inc ?\n')
            fout.write('_diffrn_ambient_pressure ?\n')
            fout.write('_diffrn_ambient_temperature ?\n')
            fout.write('_riet_lorentz_restricted true\n')
            fout.write('_riet_par_background_interpolated false\n')
            fout.write('_riet_par_background_interpolation_range 10\n')
            fout.write('_riet_meas_dataset_compute true\n')
            fout.write('_riet_meas_datafile_replace false\n')
            fout.write('_riet_meas_dataset_random_texture false\n')
            fout.write('_maud_background_add_automatic false\n')
            fout.write('_maud_interpolated_background_iterations 10\n')
            fout.write('_pd_proc_ls_datafile_weight 1.0\n')
            fout.write('_riet_meas_dataset_no_strain false\n')
            fout.write('_riet_par_background_exp_shift 0.0\n')
            fout.write('_riet_par_background_exp_thermal_shift 0.0\n')
            fout.write('_pd_spec_orientation_omega 0.0\n')
            fout.write('_pd_spec_orientation_chi 0.0\n')
            fout.write('_pd_spec_orientation_phi 0.0\n')
            fout.write('_pd_meas_orientation_omega_offset 0.0\n')
            fout.write('_pd_meas_orientation_chi_offset 0.0\n')
            fout.write('_pd_meas_orientation_phi_offset 0.0\n')
            fout.write('_pd_meas_orientation_eta_offset 0.0\n')
            fout.write('_riet_par_spec_displac_x 0.0\n')
            fout.write('_riet_par_spec_displac_y 0.0\n')
            fout.write('_riet_par_spec_displac_z 0.0\n\n')
            
            ## set bkgd function polynomial order
            fout.write('loop_\n')
            fout.write('_riet_par_background_pol\n')
            if bkgd is None:
                fout.write('  0.0\n')
                fout.write('  0.0\n')
            else:
                fout.write(f'  {bkgd[chi][0]:3.4f}\n')
                fout.write(f'  {bkgd[chi][1]:3.4f}\n')
                
            fout.write('#subordinateObject_Le Bail\n\n')
            fout.write('_riet_intensity_extraction \'Le Bail\'\n\n')
            
            fout.write('_riet_lebail_iteration_max 5\n')
            fout.write('_riet_lebail_error_max 0.0050\n')
            fout.write('_riet_lebail_range_factor 0.05\n')
            fout.write('_riet_lebail_use_bkg true\n')
            fout.write('_riet_lebail_use_hkl true\n')
            fout.write('_riet_lebail_summation_delta 1.0E-4\n\n')
            
            fout.write('#end_subordinateObject_Le Bail\n\n\n')


            fout.write('#subordinateObject_none pe\n\n')
            
            fout.write('_riet_position_extraction \'none pe\'\n\n\n')
            

            fout.write('#end_subordinateObject_none pe\n\n\n')
            
            
            fout.write('#subordinateObject_none reflectivity\n\n')
            
            fout.write('_reflectivity_model_type \'none reflectivity\'\n\n\n')
            
            
            fout.write('#end_subordinateObject_none reflectivity\n\n')
            
            fout.write('#subordinateObject_'+instName+'\n\n')
        
            fout.write('_diffrn_measurement_device '+instName+'\n\n')
        
            fout.write('_diffrn_measurement_device_type '+instName+'\n')
            fout.write('_maud_optional_intensity_factor 1.0\n')
            fout.write('_pd_proc_intensity_incident 1.0 #positive\n\n\n')
            # fout.write('loop_\n')
            # fout.write('_riet_par_2-theta_offset\n')
            # fout.write('  0\n\n\n')     
            
            fout.write('#subordinateObject_none cal\n\n')
            fout.write('_inst_intensity_calibration \'none cal\'\n\n\n')
                
            fout.write('#end_subordinateObject_none cal\n\n\n')
        
            fout.write('#subordinateObject_Instrument disalignment\n\n')
            
            fout.write('_inst_angular_calibration \'Instrument disalignment\'\n\n\n')
        
            fout.write('loop_\n')
            fout.write('_riet_par_2-theta_offset\n')
            fout.write('  0\n\n')

            fout.write('#end_subordinateObject_Instrument disalignment\n\n\n')
        
            fout.write('#subordinateObject_Debye-Scherrer\n\n')
        
            fout.write('_pd_instr_geometry \'Debye-Scherrer\'\n\n')
        
            fout.write('_diffrn_radiation_monochromator Filtered\n')
            fout.write('_pd_instr_2theta_monochr_post 0\n')
            fout.write('_pd_instr_dist_src/samp 710\n')
            fout.write('_pd_instr_monochr_pre_spec Monochromator\n')
            fout.write('_pd_instr_2theta_monochr_pre 51.5\n')
            fout.write('_pd_instr_divg_eq_src/samp 0.0\n')
            fout.write('_pd_instr_divg_slit_auto false\n')
            fout.write('_pd_instr_divg_ax_src/samp 0.0\n')
            fout.write('_diffrn_radiation_polarisn_norm 0\n')
            fout.write('_diffrn_radiation_polarisn_ratio 0\n\n')
        
            fout.write('#end_subordinateObject_Debye-Scherrer\n\n\n')
        
        
            fout.write('#subordinateObject_2Theta\n\n')
        
            fout.write('_diffrn_measurement_method \'2Theta\'\n\n\n')
        
        
            fout.write('#end_subordinateObject_2Theta\n\n\n')
        
        
            fout.write('#subordinateObject_Neutron\n\n')
        
            fout.write('_diffrn_radiation_type \'Neutron\'\n\n')
        
            fout.write('_diffrn_dynamical_scattering_correction true\n')
            fout.write('_diffrn_dynamical_scattering_correction_crystallite false\n\n')
        
            fout.write('#subordinateObject_Added_by_default\n\n')
        
            fout.write('_diffrn_radiation_wavelength_id \'Added_by_default\'\n\n')
        
            fout.write('_diffrn_radiation_wavelength {}\n'.format(wavelength))
            fout.write('_diffrn_radiation_wavelength_wt 1.0\n\n')
        
            fout.write('#end_subordinateObject_Added_by_default\n\n\n')
        
        
            fout.write('#end_subordinateObject_Neutron\n\n\n')
        
        
            fout.write('#subordinateObject_Curved Position Sensitive\n\n')
        
            fout.write('_diffrn_radiation_detector \'Curved Position Sensitive\'\n\n\n')
        
            fout.write('#end_subordinateObject_Curved Position Sensitive\n\n\n')
        
        
            fout.write('#subordinateObject_Caglioti PV\n\n')
        
            fout.write('_diffrn_inst_broadening \'Caglioti PV\'\n\n')

            fout.write('_riet_caglioti_d_dep true\n')
            fout.write('_riet_asymmetry_tan_dep true\n')
            fout.write('_riet_omega/chi_broadening_convoluted false\n')
            fout.write('_riet_par_asymmetry_truncation 0.4\n')
            fout.write('_riet_par_asymmetry_reciprocal false\n\n')

            # fout.write('loop_\n')
            # fout.write(' _riet_par_asymmetry_value\n')
            # fout.write('  {}\n'.format(caglioti['as1']))
            # fout.write('  {}\n\n'.format(caglioti['as2']))

            fout.write('loop_\n')
            fout.write('_riet_par_caglioti_value\n')
            fout.write('  {}\n'.format(caglioti['U']))
            fout.write('  {}\n'.format(caglioti['V']))
            fout.write('  {}\n\n'.format(caglioti['W']))

            fout.write('loop_\n')
            fout.write('_riet_par_gaussian_value\n')
            fout.write('  {}\n'.format(caglioti['g1']))
            fout.write('  {}\n\n\n'.format(caglioti['g2']))  
        
            fout.write('#end_subordinateObject_Caglioti PV\n\n\n')


            fout.write('#subordinateObject_none abs\n\n')

            fout.write('_exptl_absorpt_correction_type \'none abs\'\n\n\n')


            fout.write('#end_subordinateObject_none abs\n\n\n') 


            fout.write('#end_subordinateObject_{}\n\n\n'.format(instName)) 


            fout.write('#subordinateObject_none fluorescence\n\n')

            fout.write('_fluorescence_model_type \'none fluorescence\'\n\n\n')

            fout.write('#end_subordinateObject_none fluorescence\n\n\n')


            # fout.write('#subordinateObject_Basic diffraction\n\n')

            # fout.write('_diffraction_model_type \'Basic diffraction\'\n\n')

            # fout.write('#end_subordinateObject_Basic diffraction\n\n')
                        
            ## loop over files (subset of dataFrame)
            for index,run in chiTiltFiles.iterrows(): 

                if mode == 'auto':
                    ## eta only at zero
                    inner_iter = zip([diffData[run['Run #']]],[0])

                if mode == 'texture':
                    ## eta at -5, 0 and 5
                    inner_iter = zip([diffData[run['Run #']][-5],
                                      diffData[run['Run #']][0],
                                      diffData[run['Run #']][5]],[-5,0,5])

                for data,eta in inner_iter:

                    if mode == 'auto': datasetName = 'HB2C_{}.xye'.format(run['Run #'])
                    elif mode == 'texture': datasetName = 'HB2C_{}_eta_{}-0.xye'.format(run['Run #'],eta)
                    fout.write('#subordinateObject_{}\n\n'.format(datasetName))
                    fout.write('_riet_meas_datafile_name \'{}\'\n\n'.format(datasetName))
                    
                    fout.write('_riet_meas_datafile_format ?\n')
                    fout.write('_pd_meas_angle_omega 0.0\n')
                    fout.write('_pd_meas_angle_chi {}\n'.format(float(90 - chi)))
                    fout.write('_pd_meas_angle_phi {}\n'.format(float(360 - run['Phi'])))
                    fout.write('_pd_meas_angle_eta {}\n'.format(float(eta)))
                    fout.write('_pd_meas_angle_2theta 0.0\n')
                    fout.write('_pd_meas_energy_kev 0.0\n')
                    fout.write('_riet_meas_datafile_compute true\n')
                    fout.write('_riet_meas_datafile_fitting false\n')
                    fout.write('_pd_meas_detector_id none\n')
                    fout.write('_pd_meas_step_count_time 1.0\n')
                    fout.write('_pd_meas_units_of_intensity counts\n')
                    fout.write('_riet_meas_datafile_as_background false\n')
                    fout.write('_riet_meas_data_group_count 1\n')
                    fout.write('_riet_datafile_type 0\n')
                    fout.write('_riet_datafile_save_custom ?\n')
                    fout.write('_pd_meas_image_id -1\n')
                    fout.write('_riet_background_interpolated_manual false\n')
                    fout.write('_pd_meas_datetime_initiated 1984-01-01T00:00:00\n')
                    fout.write('_pd_proc_ls_datafile_weight 1.0\n')
                    fout.write('_riet_use_count_time false\n')
                    fout.write('_riet_chebyshev_polynomial_background ?\n')
                    fout.write('_pd_meas_counts_monitor 1.0 #positive\n')
                    fout.write('_riet_par_spec_displac_x_2R 0\n')
                    fout.write('_riet_par_spec_displac_z_2R 0\n')
                    fout.write('_riet_par_fluo_diffr_scale 1.0 #positive\n\n\n')

                    ## get data out
                    TT    = data['2theta']
                    Data  = data['int']
                    Error = data['esd'] 

                    fout.write('#custom_object_intensity_data\n')
                    fout.write('_pd_meas_number_of_points {}\n'.format(len(TT)))
                    fout.write('_riet_meas_datafile_calibrated false\n')
                    fout.write('_riet_meas_datafile_dspacing_based false\n')
                    fout.write('_riet_meas_datafile_energy_dispersive false\n')
                    fout.write('_riet_meas_datafile_constant_step true\n')
                    fout.write('loop_\n')
                    fout.write('_pd_meas_position _pd_meas_intensity_total _pd_meas_intensity_sigma\n')              
            
                    for index in range(len(TT)):
                        fout.write(' {:f} {:f} {:f}\n'.format(TT[index], Data[index], Error[index]))
                        
                    fout.write('\n')    
                    fout.write('#end_custom_object_intensity_data\n\n\n')

                    if fullPar is True:
                        fout.write('#custom_object_texture_factors\n')
                        fout.write('_rita_texture_points_number 1_rita_texture_radiations_number 1\n\n')
                        fout.write('data_phase_Fe(BCC)\n\n')
                        fout.write('loop_\n')
                        fout.write('_refln_index_h _refln_index_k _refln_index_l _rita_texture_factor_meas _rita_texture_factor_calc _rita_texture_point_index _rita_texture_radiation_index\n')
                        fout.write('1 1 0 0.0 0.0 0 0\n')
                        fout.write('2 0 0 0.0 0.0 0 0\n')
                        fout.write('2 1 1 0.0 0.0 0 0\n')
                        fout.write('2 2 0 0.0 0.0 0 0\n')
                        fout.write('3 1 0 0.0 0.0 0 0\n')
                        fout.write('2 2 2 0.0 0.0 0 0\n')
                        fout.write('3 2 1 0.0 0.0 0 0\n')
                        fout.write('4 0 0 0.0 0.0 0 0\n\n\n')
                        fout.write('data_phase_Fe(FCC)\n\n')
                        fout.write('loop_\n')
                        fout.write('_refln_index_h _refln_index_k _refln_index_l _rita_texture_factor_meas _rita_texture_factor_calc _rita_texture_point_index _rita_texture_radiation_index\n')
                        fout.write('1 1 1 0.0 0.0 0 0\n')
                        fout.write('2 0 0 0.0 0.0 0 0\n')
                        fout.write('2 2 0 0.0 0.0 0 0\n')
                        fout.write('3 1 1 0.0 0.0 0 0\n')
                        fout.write('2 2 2 0.0 0.0 0 0\n')
                        fout.write('4 0 0 0.0 0.0 0 0\n')
                        fout.write('3 3 1 0.0 0.0 0 0\n')
                        fout.write('4 2 0 0.0 0.0 0 0\n')
                        fout.write('4 2 2 0.0 0.0 0 0\n\n\n')
                        fout.write('#end_custom_object_texture_factors\n\n\n')
            
                    fout.write('#end_subordinateObject_{}\n\n'.format(datasetName))

            fout.write('#custom_object_Fhkl\n\n')

            fout.write('data_phase_Fe(BCC)\n\n')

            fout.write('loop_\n')
            fout.write('_refln_index_h _refln_index_k _refln_index_l _refln_F_squared_meas _refln_F_squared_calc _refln_F_squared_sigma\n') 
            fout.write('1 1 0 0.0 0.0 0\n')
            fout.write('2 0 0 0.0 0.0 0\n')
            fout.write('2 1 1 0.0 0.0 0\n')
            fout.write('2 2 0 0.0 0.0 0\n')
            fout.write('3 1 0 0.0 0.0 0\n')
            fout.write('2 2 2 0.0 0.0 0\n')
            fout.write('3 2 1 0.0 0.0 0\n')
            fout.write('4 0 0 0.0 0.0 0\n\n\n')


            fout.write('data_phase_Fe(FCC)\n\n')

            fout.write('loop_\n')
            fout.write('_refln_index_h _refln_index_k _refln_index_l _refln_F_squared_meas _refln_F_squared_calc _refln_F_squared_sigma\n') 
            fout.write('1 1 1 0.0 0.0 0\n')
            fout.write('2 0 0 0.0 0.0 0\n')
            fout.write('2 2 0 0.0 0.0 0\n')
            fout.write('3 1 1 0.0 0.0 0\n')
            fout.write('2 2 2 0.0 0.0 0\n')
            fout.write('4 0 0 0.0 0.0 0\n')
            fout.write('3 3 1 0.0 0.0 0\n')
            fout.write('4 2 0 0.0 0.0 0\n')
            fout.write('4 2 2 0.0 0.0 0\n\n\n')

            fout.write('#end_custom_object_Fhkl\n\n\n')


            fout.write('#end_subordinateObject_{}\n'.format(datasetString))

        if fullPar is True:

            fout.write('#subordinateObject_Fe(BCC)\n\n')

            fout.write('_pd_phase_name \'Fe(BCC)\'\n\n')

            fout.write('_chemical_name_common Fe(BCC)\n')
            fout.write('_chemical_formula_sum Fe2\n')
            fout.write('_symmetry_cell_setting cubic\n')
            fout.write('_symmetry_space_group_name_H-M Im-3m\n')
            fout.write('_cell_formula_units_Z 2\n')
            fout.write('_refine_ls_d_res_low 0\n')
            fout.write('_refine_ls_d_res_high 5000\n')
            fout.write('_reflns_d_resolution_low 0.7\n')
            fout.write('_reflns_d_resolution_high 50\n')
            fout.write('_maud_sg_centering_type P\n')
            fout.write('_chemical_name_mineral Unknown\n')
            fout.write('_chemical_name_systematic Unknown\n\n')

            fout.write('_cell_length_a 2.8688653 #positive #min 2.0 #max 4.0\n')
            fout.write('_cell_length_b 2.8688653 #positive #min 2.0 #max 4.0\n')
            fout.write('_cell_length_c 2.8688653 #positive #min 2.0 #max 4.0\n')
            fout.write('_cell_angle_alpha 90.00000000 #positive #min 90.0 #max 120.0\n')
            fout.write('_cell_angle_beta 90.00000000 #positive #min 90.0 #max 120.0\n')
            fout.write('_cell_angle_gamma 90.00000000 #positive #min 90.0 #max 120.0\n')
            fout.write('_riet_par_strain_thermal 0 #min -0.1 #max 0.1\n')
            fout.write('_exptl_absorpt_cryst_size 0 #positive #min 0.001 #max 100.0\n')
            fout.write('_riet_par_phase_scale_factor 1.0 #positive #min 0.0 #max 100.0\n\n\n')


            fout.write('#subordinateObject_E-WIMV\n')

            fout.write('_pd_proc_ls_pref_orient_corr \'E-WIMV\'\n\n')

            fout.write('_rita_generate_symmetry none\n')
            fout.write('_rita_wimv_sum_coincidence true\n')
            fout.write('_rita_wimv_iteration_max 10\n')
            fout.write('_rita_wimv_exponent 0.01\n')
            fout.write('_rita_wimv_refl_min_int 0.001\n')
            if wimv_cell_size == 15:
                # for 15deg res
                fout.write('_rita_wimv_odf_resolution 15.0\n')
            elif wimv_cell_size == 10:
                # for 10deg res
                fout.write('_rita_wimv_odf_resolution 10.0\n')
            elif wimv_cell_size == 5:
                # for 5deg res
                fout.write('_rita_wimv_odf_resolution 5.0\n')            
            fout.write('_rita_wimv_tube_projection true\n')
            fout.write('_rita_wimv_store_ang_conv true\n')
            fout.write('_rita_wimv_odf_coverage_% 100.0\n')
            fout.write('_rita_odf_sharpness 1.0312783712561846\n')
            fout.write('_rita_wimv_phon_use ?\n')
            fout.write('_rita_wimv_tube_weight 0.5\n')
            fout.write('_rita_wimv_normalize_pole_figures true\n')
            fout.write('_rita_wimv_weigths_exponent 0.5\n')
            fout.write('_rita_odf_refinable true\n')
            fout.write('_rita_wimv_refl_min_dspacing 0.0\n\n')

            fout.write('#custom_object_odf\n')
            fout.write('loop_\n')
            fout.write('_rita_wimv_odf_values\n')
            if wimv_cell_size == 15:
                # for 15deg res
                for i in range(7):
                    for j in range(7):
                        if j > 5: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n\n')
                        else: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n')
            elif wimv_cell_size == 10:
                # for 10deg res
                for i in range(10):
                    for j in range(10):
                        if j > 8: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n\n')
                        else: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n')  
            elif wimv_cell_size == 5:          
                ## for 5deg res
                for i in range(19):
                    for j in range(19):
                        if j > 17: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 \n\n')
                        else: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 \n')                        

            fout.write('\n#end_custom_object_odf\n\n\n')


            fout.write('#end_subordinateObject_E-WIMV\n\n\n')


            fout.write('#subordinateObject_Delft\n\n')

            fout.write('_riet_size_strain_model \'Delft\'\n\n\n')


            fout.write('#end_subordinateObject_Delft\n\n\n')


            fout.write('#subordinateObject_Isotropic\n\n')

            fout.write('_riet_sizestrain_sym_model \'Isotropic\'\n\n')

            fout.write('_riet_par_cryst_size 1000.0 #positive #min 50.0 #max 5000.0\n')
            fout.write('_riet_par_rs_microstrain 0.0006 #positive #min 0.0 #max 0.005\n')

            fout.write('#end_subordinateObject_Isotropic\n\n\n')


            fout.write('#subordinateObject_none abm\n\n')

            fout.write('_riet_antiphase_boundary \'none abm\'\n\n\n')


            fout.write('#end_subordinateObject_none abm\n\n\n')


            fout.write('#subordinateObject_none pd\n\n')

            fout.write('_riet_planar_defect \'none pd\'\n\n\n')


            fout.write('#end_subordinateObject_none pd\n\n\n')


            fout.write('#subordinateObject_no magnetic\n\n')

            fout.write('_riet_magnetic_structure_model \'no magnetic\'\n\n\n')


            fout.write('#end_subordinateObject_no magnetic\n\n\n')


            fout.write('#subordinateObject_no strain\n\n')

            fout.write('_riet_par_strain_model \'no strain\'\n\n\n')


            fout.write('#end_subordinateObject_no strain\n\n\n')


            fout.write('#subordinateObject_No microabsorption\n\n')

            fout.write('_riet_micro_absorption_model \'No microabsorption\'\n\n\n')


            fout.write('#end_subordinateObject_No microabsorption\n\n\n')


            fout.write('#subordinateObject_Atomic Structure\n\n')

            fout.write('_riet_structure_model \'Atomic Structure\'\n\n')

            fout.write('_riet_structure_quantity_from_occupancy true\n')
            fout.write('_refine_ls_energy_weight 1.0\n')
            fout.write('_riet_structure_use_U_dimensionless false\n\n')

            fout.write('#subordinateObject_No force field\n\n')

            fout.write('_riet_structure_force_field \'No force field\'\n\n\n')


            fout.write('#end_subordinateObject_No force field\n\n\n')


            fout.write('#subordinateObject_Fe0\n\n')

            fout.write('_atom_site_label \'Fe0\'\n\n')

            fout.write('_atom_site_type_symbol dummy\n')
            fout.write('_atom_site_constraints ?\n')
            fout.write('_atom_type_number_in_cell 2.0\n')
            fout.write('_atom_site_calc_flag .\n')
            fout.write('_atom_site_occupancy 1 #positive #min 0.0 #max 1.0\n')
            fout.write('_atom_site_fract_x 0 #min 0.0 #max 1.0\n')
            fout.write('_atom_site_fract_y 0 #min 0.0 #max 1.0\n')
            fout.write('_atom_site_fract_z 0 #min 0.0 #max 1.0\n')
            fout.write('_atom_site_B_iso_or_equiv 0.3272 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_11 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_22 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_33 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_12 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_13 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_23 0 #min -1.0 #max 10.0\n\n')

            fout.write('#subordinateObject_Fe\n\n')

            fout.write('_rg_site_scatterer \'Fe\'\n\n')

            fout.write('_atom_type_symbol Fe\n')
            fout.write('_atom_site_occupancy 1.0 #positive #min 0.0 #max 1.0\n\n')

            fout.write('#end_subordinateObject_Fe\n\n\n')


            fout.write('#end_subordinateObject_Fe0\n\n\n')


            fout.write('#end_subordinateObject_Atomic Structure\n\n\n')


            fout.write('#subordinateObject_atomic standard model\n\n')

            fout.write('_riet_structure_factor_model \'atomic standard model\'\n\n\n')


            fout.write('#end_subordinateObject_atomic standard model\n\n\n')


            fout.write('#subordinateObject_Le Bail\n\n')

            fout.write('_riet_structure_factor_extractor \'Le Bail\'\n\n')

            fout.write('_riet_lebail_iteration_max 5\n')
            fout.write('_riet_lebail_error_max 0.005\n')
            fout.write('_riet_lebail_range_factor 0.05\n')
            fout.write('_riet_lebail_use_bkg true\n')
            fout.write('_riet_lebail_summation_delta 1.0E-4\n')
            fout.write('_riet_lebail_use_previous_factors true\n\n')

            fout.write('#end_subordinateObject_Le Bail\n\n\n')


            fout.write('#subordinateObject_None TDS\n\n')

            fout.write('_riet_tds_model \'None TDS\'\n\n\n')


            fout.write('#end_subordinateObject_None TDS\n\n\n')


            fout.write('#end_subordinateObject_Fe(BCC)\n\n\n')


            fout.write('#subordinateObject_Fe(FCC)\n\n')

            fout.write('_pd_phase_name \'Fe(FCC)\'\n\n')

            fout.write('_chemical_name_common Fe(FCC)\n')
            fout.write('_chemical_formula_sum Fe1\n')
            fout.write('_symmetry_cell_setting cubic\n')
            fout.write('_symmetry_space_group_name_H-M Fm-3m\n')
            fout.write('_cell_formula_units_Z 4\n')
            fout.write('_refine_ls_d_res_low 0\n')
            fout.write('_refine_ls_d_res_high 5000\n')
            fout.write('_reflns_d_resolution_low 0.7\n')
            fout.write('_reflns_d_resolution_high 50\n')
            fout.write('_maud_sg_centering_type P\n')
            fout.write('_chemical_name_mineral Unknown\n')
            fout.write('_chemical_name_systematic \'Fe(FCC)\'\n\n')

            fout.write('_cell_length_a 3.599214 #positive #min 3.0 #max 5.0\n')
            fout.write('_cell_length_b 3.599214 #positive #min 3.0 #max 5.0\n')
            fout.write('_cell_length_c 3.599214 #positive #min 3.0 #max 5.0\n')
            fout.write('_cell_angle_alpha 90. #positive #min 90.0 #max 120.0\n')
            fout.write('_cell_angle_beta 90. #positive #min 90.0 #max 120.0\n')
            fout.write('_cell_angle_gamma 90. #positive #min 90.0 #max 120.0\n')
            fout.write('_riet_par_strain_thermal 0 #min -0.1 #max 0.1\n')
            fout.write('_exptl_absorpt_cryst_size 0 #positive #min 0.001 #max 100.0\n')
            fout.write('_riet_par_phase_scale_factor 1.0 #positive #min 0.0 #max 100.0\n\n\n')


            fout.write('#subordinateObject_E-WIMV\n\n')

            fout.write('_pd_proc_ls_pref_orient_corr \'E-WIMV\'\n\n')

            fout.write('_rita_generate_symmetry none\n')
            fout.write('_rita_wimv_sum_coincidence true\n')
            fout.write('_rita_wimv_iteration_max 10\n')
            fout.write('_rita_wimv_exponent 0.01\n')
            fout.write('_rita_wimv_refl_min_int 0.001\n')
            if wimv_cell_size == 15:
                # for 15deg res
                fout.write('_rita_wimv_odf_resolution 15.0\n')
            elif wimv_cell_size == 10:
                # for 10deg res
                fout.write('_rita_wimv_odf_resolution 10.0\n')
            elif wimv_cell_size == 5:
                # for 5deg res
                fout.write('_rita_wimv_odf_resolution 5.0\n')            
            fout.write('_rita_wimv_tube_projection true\n')
            fout.write('_rita_wimv_store_ang_conv true\n')
            fout.write('_rita_wimv_odf_coverage_% 100.0\n')
            fout.write('_rita_odf_sharpness 1.0312783712561846\n')
            fout.write('_rita_wimv_phon_use ?\n')
            fout.write('_rita_wimv_tube_weight 0.5\n')
            fout.write('_rita_wimv_normalize_pole_figures true\n')
            fout.write('_rita_wimv_weigths_exponent 0.5\n')
            fout.write('_rita_odf_refinable true\n')
            fout.write('_rita_wimv_refl_min_dspacing 0.0\n\n')

            fout.write('#custom_object_odf\n')
            fout.write('loop_\n')
            fout.write('_rita_wimv_odf_values\n')
            if wimv_cell_size == 15:
                # for 15deg res
                for i in range(7):
                    for j in range(7):
                        if j > 5: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n\n')
                        else: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n')
            elif wimv_cell_size == 10:
                # for 10deg res
                for i in range(10):
                    for j in range(10):
                        if j > 8: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n\n')
                        else: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n')  
            elif wimv_cell_size == 5:          
                ## for 5deg res
                for i in range(19):
                    for j in range(19):
                        if j > 17: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 \n\n')
                        else: fout.write('1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 \n')                        


            fout.write('\n#end_custom_object_odf\n\n\n')


            fout.write('#end_subordinateObject_E-WIMV\n\n\n')


            fout.write('#subordinateObject_Delft\n\n')

            fout.write('_riet_size_strain_model \'Delft\'\n\n\n')


            fout.write('#end_subordinateObject_Delft\n\n\n')


            fout.write('#subordinateObject_Isotropic\n\n')

            fout.write('_riet_sizestrain_sym_model \'Isotropic\'\n\n')

            fout.write('_riet_par_cryst_size 1000.0 #positive #min 50.0 #max 5000.0\n')
            fout.write('_riet_par_rs_microstrain 0.0006 #positive #min 0.0 #max 0.005\n')

            fout.write('#end_subordinateObject_Isotropic\n\n\n')


            fout.write('#subordinateObject_none abm\n\n')

            fout.write('_riet_antiphase_boundary \'none abm\'\n\n\n')


            fout.write('#end_subordinateObject_none abm\n\n\n')


            fout.write('#subordinateObject_none pd\n\n')

            fout.write('_riet_planar_defect \'none pd\'\n\n\n')


            fout.write('#end_subordinateObject_none pd\n\n\n')


            fout.write('#subordinateObject_no magnetic\n\n')

            fout.write('_riet_magnetic_structure_model \'no magnetic\'\n\n\n')


            fout.write('#end_subordinateObject_no magnetic\n\n\n')


            fout.write('#subordinateObject_no strain\n\n')

            fout.write('_riet_par_strain_model \'no strain\'\n\n\n')


            fout.write('#end_subordinateObject_no strain\n\n\n')


            fout.write('#subordinateObject_No microabsorption\n\n')

            fout.write('_riet_micro_absorption_model \'No microabsorption\'\n\n\n')


            fout.write('#end_subordinateObject_No microabsorption\n\n\n')


            fout.write('#subordinateObject_Atomic Structure\n\n')

            fout.write('_riet_structure_model \'Atomic Structure\'\n\n')

            fout.write('_riet_structure_quantity_from_occupancy true\n')
            fout.write('_refine_ls_energy_weight 1.0\n')
            fout.write('_riet_structure_use_U_dimensionless false\n\n')

            fout.write('#subordinateObject_No force field\n\n')

            fout.write('_riet_structure_force_field \'No force field\'\n\n\n')


            fout.write('#end_subordinateObject_No force field\n\n\n')


            fout.write('#subordinateObject_Fe0\n\n')

            fout.write('_atom_site_label \'Fe0\'\n\n')

            fout.write('_atom_site_type_symbol dummy\n')
            fout.write('_atom_site_constraints ?\n')
            fout.write('_atom_type_number_in_cell 2.0\n')
            fout.write('_atom_site_calc_flag .\n')
            fout.write('_atom_site_occupancy 1 #positive #min 0.0 #max 1.0\n')
            fout.write('_atom_site_fract_x 0 #min 0.0 #max 1.0\n')
            fout.write('_atom_site_fract_y 0 #min 0.0 #max 1.0\n')
            fout.write('_atom_site_fract_z 0 #min 0.0 #max 1.0\n')
            fout.write('_atom_site_B_iso_or_equiv 0.5615 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_11 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_22 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_33 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_12 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_13 0 #min -1.0 #max 10.0\n')
            fout.write('_atom_site_aniso_B_23 0 #min -1.0 #max 10.0\n\n')

            fout.write('#subordinateObject_Fe\n\n')

            fout.write('_rg_site_scatterer \'Fe\'\n\n')

            fout.write('_atom_type_symbol Fe\n')
            fout.write('_atom_site_occupancy 1.0 #positive #min 0.0 #max 1.0\n\n')

            fout.write('#end_subordinateObject_Fe\n\n\n')


            fout.write('#end_subordinateObject_Fe0\n\n\n')


            fout.write('#end_subordinateObject_Atomic Structure\n\n\n')


            fout.write('#subordinateObject_atomic standard model\n\n')

            fout.write('_riet_structure_factor_model \'atomic standard model\'\n\n\n')


            fout.write('#end_subordinateObject_atomic standard model\n\n\n')


            fout.write('#subordinateObject_Le Bail\n\n')

            fout.write('_riet_structure_factor_extractor \'Le Bail\'\n\n')

            fout.write('_riet_lebail_iteration_max 5\n')
            fout.write('_riet_lebail_error_max 0.005\n')
            fout.write('_riet_lebail_range_factor 0.05\n')
            fout.write('_riet_lebail_use_bkg true\n')
            fout.write('_riet_lebail_summation_delta 1.0E-4\n')
            fout.write('_riet_lebail_use_previous_factors true\n')

            fout.write('#end_subordinateObject_Le Bail\n\n\n')


            fout.write('#subordinateObject_None TDS\n\n')

            fout.write('_riet_tds_model \'None TDS\'\n\n\n')


            fout.write('#end_subordinateObject_None TDS\n\n\n')


            fout.write('#end_subordinateObject_Fe(FCC)\n\n\n')

rootDir = '/' #

# general infor
instName   = 'WAND2'
mode       = 'texture' ## 'auto' = autoreduced | 'texture' = in/out of plane separate
wavelength = 1.4882562

# old: 0.5651981, -1.5459125, 2.665404
caglioti   = {'U':0.5328423,
              'V':-1.3918402,
              'W':2.4130297,
              'g1':0.021636,
              'g2':0.001913,
              'as1':0,
              'as2':0} # U,V,W

# min and max bounds for MAUD computation
twotheta_min = 0
twotheta_max = 0

rotations = {9:{'omega':-90.0,
                 'chi':0.0,
                 'phi':-105.0},
             6:{'omega':-90.0,
                 'chi':0.0,
                 'phi':0.0},
             3:{'omega':-90.0,
                 'chi':0.0,
                 'phi':-25.0},}

""" 
needs columns: 
[Time, Run #, Chi, Phi, Sample ]  
"""

dataDir     = os.path.join(rootDir,'texturereduced')
outDir      = os.path.join(rootDir,'maud','samples')
expInfoFile = os.path.join(rootDir,'expSummary.csv')
expInfo     = pd.read_csv(expInfoFile)

## setup run sets for fitting
samples = pd.unique(expInfo['Sample'])

for sample in samples:

    # get metadata for a given sample
    sampleInfo = expInfo.loc[expInfo['Sample'] == sample]
    start      = sampleInfo['Run #'].min()
    end        = sampleInfo['Run #'].max()    

    # make this play nice with unix
    sample_desc = sample.replace(' #','_')
    
    if 'Steel' in sample_desc: 
        sample_num  = int(sample_desc.split('_')[1])
        smplRot = rotations[sample_num]
        fullPar   = True
        writeMAUD_WAND2(sample_desc, rootDir, dataDir, sampleInfo, start, end, caglioti, mode, wimv_cell_size=5, smplRot=smplRot, fullPar=fullPar)
        writeMAUD_WAND2(sample_desc, rootDir, dataDir, sampleInfo, start, end, caglioti, mode, wimv_cell_size=10, smplRot=smplRot, fullPar=fullPar)
        writeMAUD_WAND2(sample_desc, rootDir, dataDir, sampleInfo, start, end, caglioti, mode, wimv_cell_size=15, smplRot=smplRot, fullPar=fullPar)
     
