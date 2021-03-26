#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script to be called once some data was processed using the MRSData2 class and is is currently in memory.

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

import pdb

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

# data to process is in here
db_filepath = "/home/tangir/crmbm/acq_db/brain.pkl"

# display real-time fit and other stuff?
display_stuff = True

# %% select datasets via dataframe

df = pd.read_pickle(db_filepath)
df = df.iloc[0]
# df = df.loc[(df["patient"] == 319) & (df["study"] == 1)]
# df = df.loc[(df.index == "3695a02a4fbf0153215cf09eaa21658e")]

if(type(df) is pd.core.series.Series):
    df = df.to_frame().T

# %% prepare simulation machine

# metabolite db
meta_bs = sim.metabolite_basis_set()

# %% data fit strategies

# list to save different strategies
fit_ws_list = []

# list of sequence to try
# sequence_list = [None, sim.mrs_seq_press]
sequence_list = [sim.mrs_seq_press]

# --- easy strategy: singlets ---
metabolites2fit = np.sort([
    xxx.m_Ala,
    xxx.m_Asp,
    xxx.m_PCr,
    xxx.m_GABA,
    xxx.m_Gsh,
    xxx.m_Lac,
    xxx.m_NAAG,
    xxx.m_sI,
    xxx.m_Tau,
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_mI,
    xxx.m_Glu,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock = np.full([len(meta_bs), 4], 1)
# link singlets by linewidth and phase
linklock[metabolites2fit, :] = [0, 200, 300, 100]
linklock[xxx.m_Cr_CH3, :] = [0, -200, -300, -100]
# leave water free
linklock[xxx.m_Water, :] = [0, 0, 0, 0]
# leave Lipids linewidth free
linklock[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for this_sequence in sequence_list:
    this_fit = fit.fit_pastis(meta_bs=meta_bs)
    this_fit.name = "singlets" + "_" + str(this_sequence)
    this_fit.metabolites_area_integration = [xxx.m_NAA_CH3, xxx.m_Cr_CH3, xxx.m_Cho_CH3]
    this_fit.area_integration_peak_ranges = [0.1, 0.1, 0.1]
    this_fit.metabolites = metabolites2fit
    this_fit.params_linklock = linklock
    this_fit.sequence = this_sequence

    # default fitting bounds from template
    this_fit.params_min = this_fit.params_min.set_default_min().add_macromolecules_min()
    this_fit.params_max = this_fit.params_max.set_default_max().add_macromolecules_max()

    # re-init params_init
    this_fit.params_init = (this_fit.params_min + this_fit.params_max) / 2.0

    # ranges for concentration
    this_fit.params_min[:, xxx.p_cm] = 0.09
    this_fit.params_max[:, xxx.p_cm] = 300.0
    this_fit.params_max[xxx.m_All_MMs, xxx.p_cm] = 5000.0
    this_fit.params_max[xxx.m_Water, xxx.p_cm] = 100000.0

    # linewidth bounds for metabolites
    this_fit.params_min[:, xxx.p_dd] = 5

    # let the MMs go above
    this_fit.params_max[xxx.m_All_MMs, xxx.p_dd] = 400
    # with early minimal damping
    this_fit.params_init[:, xxx.p_dd] = this_fit.params_min[:, xxx.p_dd] * 1.1

    # initial concentrations
    this_fit.params_init[:, xxx.p_cm] = 0.0
    this_fit.params_init[metabolites2fit, xxx.p_cm] = 0.1

    # frequency shifts for metabolites
    # TODO: understand why I have a 0.05ppm systematic shift
    this_fit.params_init[:, xxx.p_df] = -25.0
    this_fit.params_min[:, xxx.p_df] = -30.0
    this_fit.params_max[:, xxx.p_df] = 15.0
    # a bit more for the lipids
    this_fit.params_min[:, xxx.p_df] = -55.0
    this_fit.params_max[:, xxx.p_df] = 25.0

    # phase bounds for metabolites and lipids
    this_fit.params_min[:, xxx.p_dp] = -0.1
    this_fit.params_max[:, xxx.p_dp] = +0.1

    # numerical optimization parameters
    this_fit.display_CRBs = False
    this_fit.display_subplots_types = [fit.fit_plot_type.BARGRAPH_CM, fit.fit_plot_type.BARGRAPH_DD, fit.fit_plot_type.BARGRAPH_DF, fit.fit_plot_type.BARGRAPH_DP]
    this_fit.display_enable = display_stuff
    this_fit.display_range_ppm = [0, 6]
    this_fit.display_frequency = 2

    # store
    fit_ws_list.append(this_fit)

# --- lcmodel strategy ---
this_fit = fit.fit_lcmodel()
this_fit.name = "lcmodel"
this_fit.display_enable = display_stuff
# store
fit_ws_list.append(this_fit)

# water concentration assumption
water_concentration = 55000.0  # mmol/kg

# %% water data fit strategy

fit_nows = fit.fit_pastis(meta_bs=meta_bs)

fit_nows.name = "water fit"
fit_nows.metabolites = metabolites2fit
fit_nows.area_integration_peaks = [xxx.m_Water]
fit_nows.area_integration_peak_ranges = [0.5]

# peak area integration

# min/max fitting bounds
fit_nows.params_min = fit_nows.params_min.set_default_water_min()
fit_nows.params_max = fit_nows.params_max.set_default_water_max()
# water concentration parameter
fit_nows.params_min[xxx.m_Water, xxx.p_cm] = 0
fit_nows.params_max[xxx.m_Water, xxx.p_cm] = 1000000.0
# water linewidth
fit_nows.params_min[xxx.m_Water, xxx.p_dd] = 10.0
fit_nows.params_max[xxx.m_Water, xxx.p_dd] = 500.0
# water frequency shift
fit_nows.params_min[xxx.m_Water, xxx.p_df] = -50.0
fit_nows.params_max[xxx.m_Water, xxx.p_df] = +50.0
# water phase
fit_nows.params_min[xxx.m_Water, xxx.p_dp] = -0.1
fit_nows.params_max[xxx.m_Water, xxx.p_dp] = +0.1

# initial fit values
fit_nows.params_init = fit_nows.params_min.copy()
fit_nows.params_init[xxx.m_Water, xxx.p_cm] = 1.0
fit_nows.params_init[xxx.m_Water, xxx.p_dd] = 100.0
fit_nows.params_init[xxx.m_Water, xxx.p_df] = 0.0
fit_nows.params_init[xxx.m_Water, xxx.p_dp] = 0.0

# numerical optimization parameters
fit_nows.optim_xtol = 1e-9
fit_nows.optim_ftol = 1e-9
fit_nows.optim_gtol = 1e-9

# display
fit_nows.display_range_ppm = [3, 7]
fit_nows.display_enable = display_stuff

# %% prepare data to fit

# list to save fit results
fit_results_list = []

# browse though datasets
for this_index, this_row in df.iterrows():

    this_raw_data = this_row["dataset_raw_data_obj"]
    this_dcm_data = this_row["dataset_dcm_data_obj"]

    # get the data
    if(this_raw_data is None):
        # no raw data? ok get the dicom
        this_data = this_dcm_data
    else:
        this_data = this_raw_data

    # display the data
    if(display_stuff):
        this_data.display_spectrum_1d()

    # removing water and any artefact > 5ppm (corrupts fit quality criteria)
    # this_data = this_data.correct_water_removal_1d(8, [4.5, 6], display_stuff)

    # filtering
    # this_data = this_data.correct_bandpass_filtering_1d([0, 6], np.ones, display_stuff)

    # reapodize to remove filtering artefact (will not affect linewidth because already apodized)
    # this_data = this_data.correct_apodization_nd(5.0, display_stuff)

    # %% fit non water-suppressed data

    this_fit_nows = fit_nows.copy()
    this_fit_nows.data = this_data.data_ref.copy()
    this_fit_nows.run()
    this_fit_nows_df = this_fit_nows.to_dataframe('nows_')

    # %% use various fit strategies
    for fit_ws in fit_ws_list:

        # --- fit water-suppressed data ---
        this_fit_ws = fit_ws.copy()
        this_fit_ws.data = this_data.copy()

        if(type(this_fit_ws) is fit.fit_pastis):
            # linewidth bounds for metabolites
            # use water linewidth for max
            this_fit_ws.params_max[:, xxx.p_dd] = this_fit_nows.params_fit[xxx.m_Water, xxx.p_dd] + 5

        # run the fit
        this_fit_ws.run()
        this_fit_ws_df = this_fit_ws.to_dataframe('ws_')

        # --- save the fit results ---
        this_fit_df = pd.concat([this_fit_ws_df, this_fit_nows_df], axis=1)
        this_fit_df = this_fit_df.set_index(["ws_fit_hash", "ws_data_file_hash"])
        fit_results_list.append(this_fit_df)


# append all the results and store
df_fit_results = pd.concat(fit_results_list, axis=0)
db_filepath_noext, _ = os.path.splitext(db_filepath)
df_fit_results.to_pickle(db_filepath_noext + "_fit.pkl")
