#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script to fit data previously processed with pastis and stored in a pkl file.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.fit as fit
import mrs.sim as sim
import mrs.aliases as xxx
import mrs.log as log
import numpy as np
import pandas as pd
import os
import itertools
from datetime import datetime

import pdb

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.INFO)

# data to process is in here
db_filepath = "/home/tangir/crmbm/acq_db/brain.pkl"
fit_suffix_pkl = "_fit"

# display real-time fit and other stuff?
display_stuff = False
display_range_ppm_ws = [0, 6]
display_range_ppm_nows = [4, 6]

# remove residual water with HLSVD?
remove_residual_water_range_ppm = [4.3, 6]
# remove lipids with HLSVD?
remove_lipids_range_ppm = [0, 1.8]
                                                      
# fit ppm ranges (experimental)
fit_ppm_range_ws = [0, 4.2]
fit_ppm_range_nows = None

# %% select datasets via dataframe

df = pd.read_pickle(db_filepath)

#df = df.loc["b057a3e07319341e032e3c6c36ec4d83"]
#df = df.iloc[3]

#df = df.loc[df["reco_dataset_raw_data_name"].str.contains("314").fillna(False)]
#df = df.iloc[2]

# keep a dataframe type, even if one line
if(type(df) is pd.core.series.Series):
    df = df.to_frame().T

print(df)

# %% prepare simulation machine

# metabolite db for metabolites
meta_bs = sim.metabolite_basis_set()
# for water
meta_bs_nows = sim.metabolite_basis_set()

# %% data fit strategies

# list to save different strategies
fit_ws_list = []

# --- list of sequence to try ---
sequence_list = [None]
sequence_list.append(sim.mrs_seq_press)

# --- list of metabolite lists to try ---
metabolites_list_list = []

metabolites_list_list.append(np.sort([
    xxx.m_LipA,
    xxx.m_LipB,
    xxx.m_LipC,
    xxx.m_NAA,
    xxx.m_NAAG,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_PCr,
    xxx.m_GPC,
    xxx.m_PC,
    xxx.m_mI,
    xxx.m_Glu,
    xxx.m_Gln,
    xxx.m_Gsh,
    xxx.m_Asp,
    xxx.m_Tau,
    xxx.m_GABA,
    xxx.m_Glc]))

# --- without Glc, GABA, Asp ---
metabolites_list_list.append(np.sort([
    xxx.m_LipA,
    xxx.m_LipB,
    xxx.m_LipC,
    xxx.m_NAA,
    xxx.m_NAAG,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_PCr,
    xxx.m_GPC,
    xxx.m_PC,
    xxx.m_mI,
    xxx.m_Glu,
    xxx.m_Gln,
    xxx.m_Gsh,
    xxx.m_Tau]))

# --- without Tau ---
metabolites_list_list.append(np.sort([
    xxx.m_LipA,
    xxx.m_LipB,
    xxx.m_LipC,
    xxx.m_NAA,
    xxx.m_NAAG,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_PCr,
    xxx.m_GPC,
    xxx.m_PC,
    xxx.m_mI,
    xxx.m_Glu,
    xxx.m_Gln,
    xxx.m_Gsh]))

constraint_shrink_list = [0.5, 1]
constraint_type_list = ["lastallfree", "allfree", "2groups", "4groups", "alllinked"]

param_big_list = itertools.product(metabolites_list_list, sequence_list, constraint_shrink_list, constraint_type_list)

# create all fit strategies
for (this_met_list, this_seq, this_constraint_shrink, this_constraint_type) in param_big_list:

    # linklock: relations between fit parameters
    linklock = np.full([len(meta_bs), 4], 1)
    # link singlets by linewidth and phase
    linklock[this_met_list, :] = [0, 300, 200, 100]
    linklock[xxx.m_mI, :] = [0, -300, -200, -100]
    
    if(this_constraint_type == "allfree"):
        linklock[this_met_list, xxx.p_dd] = 0        
        
    elif(this_constraint_type == "2groups"):
        linklock[xxx.m_NAA, :] = [0, 400, 200, 100]
        linklock[xxx.m_NAAG, :] = [0, 400, 200, 100]
        linklock[xxx.m_PC, :] = [0, 400, 200, 100]
        linklock[xxx.m_GPC, :] = [0, 400, 200, 100]
        
        linklock[xxx.m_Cr_CH3, :] = [0, -400, 200, 100]
        #linklock[xxx.m_Cr_CH2, :] = [0, 0, 200, 100]
        
    elif(this_constraint_type == "4groups"):
        linklock[xxx.m_NAA, :] = [0, -400, 200, 100]
        linklock[xxx.m_NAAG, :] = [0, 400, 200, 100]
        
        linklock[xxx.m_PC, :] = [0, -500, 200, 100]
        linklock[xxx.m_GPC, :] = [0, 500, 200, 100]
        
        linklock[xxx.m_Cr_CH3, :] = [0, 0, 200, 100]
        #linklock[xxx.m_Cr_CH2, :] = [0, 0, 200, 100]
        
    elif(this_constraint_type == "alllinked"):
        pass
        
    elif(this_constraint_type == "lastallfree"):
        pass
    
    else:
        log.error("weird constraint type!")
    
    # leave water free
    #linklock[xxx.m_Water, :] = [0, 0, 0, 0]
    # leave MMs linewidth free
    linklock[xxx.m_All_MMs, :] = [0, 0, 0, 100]

    this_fit = fit.fit_pastis(meta_bs=meta_bs)
    this_fit.metabolites_area_integration = [xxx.m_NAA_CH3, xxx.m_Cr_CH3, xxx.m_Cho_CH3]
    this_fit.area_integration_peak_ranges = [0.1, 0.1, 0.1]
    this_fit.metabolites = this_met_list
    this_fit.sequence = this_seq

    # default fitting bounds from template
    this_fit.params_min = this_fit.params_min.set_default_min().add_macromolecules_min()
    this_fit.params_max = this_fit.params_max.set_default_max().add_macromolecules_max()

    # re-init params_init
    this_fit.params_init = (this_fit.params_min + this_fit.params_max) / 2.0

    # ranges for concentration
    this_fit.params_min[:, xxx.p_cm] = 0.0
    this_fit.params_max[:, xxx.p_cm] = 200.0
    this_fit.params_max[xxx.m_All_MMs, xxx.p_cm] = 1000.0
    #this_fit.params_max[xxx.m_Water, xxx.p_cm] = 100000.0

    # linewidth bounds for metabolites
    this_fit.params_min[:, xxx.p_dd] = 5

    # let the MMs go above
    this_fit.params_min[xxx.m_All_MMs, xxx.p_dd] = 150
    this_fit.params_max[xxx.m_All_MMs, xxx.p_dd] = 300
    # start with very narrow peaks / long-T2 FID
    this_fit.params_init[:, xxx.p_dd] = this_fit.params_min[:, xxx.p_dd] * 1.1

    # initial concentrations
    this_fit.params_init[:, xxx.p_cm] = 0.0
    this_fit.params_init[this_met_list, xxx.p_cm] = 0.1

    # frequency shifts for metabolites and MMs
    this_fit.params_init[:, xxx.p_df] = 0.0
    this_fit.params_min[:, xxx.p_df] = -10.0
    this_fit.params_max[:, xxx.p_df] = 10.0

    # phase bounds for all
    this_fit.params_min[:, xxx.p_dp] = -0.1
    this_fit.params_max[:, xxx.p_dp] = +0.1

    # only now, set the LL vector
    this_fit.params_linklock = linklock
    this_fit._set_unique_linklock()

    # numerical optimization parameters
    this_fit.optim_xtol = 1e-6
    this_fit.optim_ftol = 1e-6
    this_fit.optim_gtol = 1e-6
    this_fit.optim_ppm_range = fit_ppm_range_ws
    this_fit.fqn_noise_range = [-2, 0]

    # display parameters
    this_fit.display_enable = display_stuff
    this_fit.display_range_ppm = display_range_ppm_ws
    this_fit.display_frequency = 2

    # name
    this_fit.name = "pastis_" + str(len(this_met_list)) + "metabolites_" + str(this_seq) + "_" + str(this_constraint_shrink)

    # extra parameters (experimental)
    this_fit.__isfrozen = False
    this_fit.constraint_shrink = this_constraint_shrink
    this_fit.constraint_type = this_constraint_type
    this_fit.__isfrozen = True

    # store
    fit_ws_list.append(this_fit)

# --- lcmodel strategy ---
this_fit = fit.fit_lcmodel()
this_fit.name = "lcmodel"
this_fit.display_enable = display_stuff
# store
fit_ws_list.append(this_fit)

# %% water data fit strategy

fit_nows = fit.fit_pastis()

fit_nows.name = "water fit"
fit_nows.metabolites = metabolites_list_list[0]

# peak area integration
fit_nows.area_integration_peaks = [xxx.m_Water]
fit_nows.area_integration_peak_ranges = [0.5]

# min/max fitting bounds
fit_nows.params_min = fit_nows.params_min.set_default_water_min()
fit_nows.params_max = fit_nows.params_max.set_default_water_max()
fit_nows.params_max[xxx.m_Water, xxx.p_cm] = 1e6
# water linewidth
fit_nows.params_min[xxx.m_Water, xxx.p_dd] = 10.0
fit_nows.params_max[xxx.m_Water, xxx.p_dd] = 100.0
# water frequency shift
fit_nows.params_min[xxx.m_Water, xxx.p_df] = -50.0
fit_nows.params_max[xxx.m_Water, xxx.p_df] = +50.0

# initial fit values
fit_nows.params_init = fit_nows.params_min.copy()
fit_nows.params_init[xxx.m_Water, xxx.p_cm] = 1.0
fit_nows.params_init[xxx.m_Water, xxx.p_dd] = 11.0
fit_nows.params_init[xxx.m_Water, xxx.p_df] = 0.0
fit_nows.params_init[xxx.m_Water, xxx.p_dp] = 0.0

fit_nows.optim_ppm_range = fit_ppm_range_nows

# display
fit_nows.display_range_ppm = display_range_ppm_nows
fit_nows.display_enable = display_stuff

# %% prepare data to fit

# list to save fit results
fit_results_list = []

# for ETA estimation
start_time = datetime.today()

# browse though datasets
for this_row_i, (this_index, this_row) in enumerate(df.iterrows()):

    this_raw_data = this_row["reco_dataset_raw_data_obj"]
    this_dcm_data = this_row["reco_dataset_dcm_data_obj"]

    # get the data
    if(np.isnan(this_raw_data).any()):
        # no raw data? ok get the dicom
        this_data = this_dcm_data
    else:
        this_data = this_raw_data

    # display the data
    if(display_stuff):
        this_data.display_spectrum_1d()
        
    # %% pre-processing

    # removing water and any artefact > 5ppm (corrupts fit quality criteria)
    if(remove_residual_water_range_ppm is not None):
        this_data = this_data.correct_peak_removal_1d(50, remove_residual_water_range_ppm, display=display_stuff)

    # removing lipids
    if(remove_lipids_range_ppm is not None):
        this_data = this_data.correct_peak_removal_1d(50, remove_lipids_range_ppm, display=display_stuff)
        
    # recalibrating water spectrum
    this_data_ref = this_data.data_ref.copy()
    this_data_ref = this_data_ref.correct_freqshift_1d()

    # %% fit non water-suppressed data

    this_fit_nows = fit_nows.copy()
    this_fit_nows.data = this_data_ref
    this_fit_nows.run()
    this_fit_nows_df = this_fit_nows.to_dataframe('fit_nows_')

    this_linewidth_estimated = this_fit_nows.params_fit[xxx.m_Water, xxx.p_dd]

    # %% use various fit strategies
    for fit_ws_i, fit_ws in enumerate(fit_ws_list):

        # --- fit water-suppressed data ---
        this_fit_ws = fit_ws.copy()
        this_fit_ws.data = this_data.copy()

        if(type(this_fit_ws) is fit.fit_pastis):
            
            # get extra parameters
            this_constraint_shrink = this_fit_ws.constraint_shrink
    
            # prepare initial LW ranges according to water FWHM
            this_fit_ws.params_min[xxx.m_All_MBs, xxx.p_dd] = this_linewidth_estimated * 0.3
            this_fit_ws.params_max[xxx.m_All_MBs, xxx.p_dd] = this_linewidth_estimated * 1.7
            this_fit_ws.params_init[xxx.m_All_MBs, xxx.p_dd] = this_linewidth_estimated * 0.6

            # run the fit
            this_fit_ws.run()

            # number of fit runs
            n_iter = 2
            
            for i in range(n_iter):

                # check previous parameter ranges
                params_min_max_range = this_fit_ws.params_max[this_fit_ws.metabolites, :] - this_fit_ws.params_min[this_fit_ws.metabolites, :]
                # reduce these ranges
                params_min_max_range[:] = params_min_max_range[:] * this_constraint_shrink
                
                # set the min
                this_fit_ws.params_min[this_fit_ws.metabolites, :] = this_fit_ws.params_fit[this_fit_ws.metabolites, :] - params_min_max_range / 2.0
                # set the max
                this_fit_ws.params_max[this_fit_ws.metabolites, :] = this_fit_ws.params_fit[this_fit_ws.metabolites, :] + params_min_max_range / 2.0
                
                # with an exception for the concentration
                this_fit_ws.params_min[:, xxx.p_cm] = 0.0
                this_fit_ws.params_max[this_fit_ws.metabolites, xxx.p_cm] = this_fit_ws.params_fit[this_fit_ws.metabolites, xxx.p_cm] + params_min_max_range[:, xxx.p_cm] / 2.0
                
                # force minimal dd
                #this_fit_ws.params_min[xxx.m_All_MBs, xxx.p_dd] = this_linewidth_estimated * 0.3
                
                # release df
                this_fit_ws.params_linklock[this_fit_ws.metabolites, xxx.p_df] = 0
                
                if((i == (n_iter - 1)) and (this_fit_ws.constraint_type == "lastallfree")):
                   # release all dd for this last fit
                   this_fit_ws.params_linklock[this_fit_ws.metabolites, xxx.p_dd] = 0
                
                this_fit_ws.params_init[this_fit_ws.metabolites, :] = this_fit_ws.params_fit[this_fit_ws.metabolites, :]
                
                ind_init_below_min = (this_fit_ws.params_init[this_fit_ws.metabolites, :] <= this_fit_ws.params_min[this_fit_ws.metabolites, :])
                this_fit_ws.params_init[this_fit_ws.metabolites, :][ind_init_below_min] = this_fit_ws.params_min[this_fit_ws.metabolites, :][ind_init_below_min] + params_min_max_range[ind_init_below_min] * 0.01
                
                ind_init_above_max = (this_fit_ws.params_init[this_fit_ws.metabolites, :] >=  this_fit_ws.params_max[this_fit_ws.metabolites, :])
                this_fit_ws.params_init[this_fit_ws.metabolites, :][ind_init_above_max] = this_fit_ws.params_max[this_fit_ws.metabolites, :][ind_init_above_max] - params_min_max_range[ind_init_above_max] * 0.01

                # run the fit
                this_fit_ws.run()
        else:
            # if LCModel fit, just run it
            this_fit_ws.run()

        # --- save the fit results ---
        this_fit_ws_df = this_fit_ws.to_dataframe('fit_ws_')
        this_fit_df = pd.concat([this_fit_ws_df, this_fit_nows_df], axis=1)
        this_fit_df = this_fit_df.set_index(["fit_ws_fit_hash", "fit_ws_data_file_hash"])
        fit_results_list.append(this_fit_df)

        # --- store progress ---
        prct_done = (this_row_i * len(fit_ws_list) + fit_ws_i) / (len(df) * len(fit_ws_list)) * 100.0
        
        # ETA estimation
        if(prct_done > 0):
            current_time = datetime.today()
            eta_time = (100 - prct_done) * (current_time - start_time) / prct_done + start_time
            os.system("echo " + str(datetime.now()) + " " + os.path.basename(__file__) + "/" + os.path.basename(db_filepath) + ": %.1f%% ETA = %s >> progress.txt" % (prct_done, str(eta_time)))
        else:
            os.system("echo " + str(datetime.now()) + " " + os.path.basename(__file__) + "/" + os.path.basename(db_filepath) + ": %.1f%% >> progress.txt" % prct_done)

        
# append all the results and store
df_fit_results = pd.concat(fit_results_list, axis=0)
db_filepath_noext, _ = os.path.splitext(db_filepath)
df_fit_results.to_pickle(db_filepath_noext + fit_suffix_pkl + ".pkl")
