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
import mrs.db as db
import mrs.log as log
import numpy as np

import pdb

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = db.data_db("/home/tangir/crmbm/acq_db/brain.pkl")

# fast test: only last fit strategy only
only_one_fit_strategy = False

# display real-time fit and other stuff?
display_stuff = False

# %% select datasets via dataframe

df_sel = rdb.df_reco
# df_sel = rdb.df_reco.loc[(rdb.df_reco["patient"] == 319) & (rdb.df_reco["study"] == 1)]
# df_sel = rdb.df_reco.loc[(rdb.df_reco.index == "3695a02a4fbf0153215cf09eaa21658e")]

# %% fit strategies

# list to save different strategies
fit_stategies_list = []

# list of sequence to try
fit_strategies_seq_list = [sim.mrs_sequence, sim.mrs_seq_press]
#fit_strategies_seq_list = [sim.mrs_seq_press]

# init metabolite basis set and linklock
meta_bs = sim.metabolite_basis_set()
meta_bs.initialize()
linklock_arr_unset = np.full([len(meta_bs), 4], 1)

# --- easy strategy: singlets ---
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_Cr_CH3,
    xxx.m_Cho_CH3,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock_arr = linklock_arr_unset.copy()
# link singlets by linewidth and phase
linklock_arr[metabolites_fit, :] = [0, 200, 0, 100]
linklock_arr[xxx.m_Cr_CH3, :] = [0, -200, 0, -100]
# leave water free
linklock_arr[xxx.m_Water, :] = [0, 0, 0, 0]
# leave Lipids linewidth free
linklock_arr[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for s in fit_strategies_seq_list:
    suffix = "_" + str(s)
    this_strategy = fit.fit_stategy("singlets" + suffix, metabolites_fit, linklock_arr, s)
    fit_stategies_list.append(this_strategy)

# --- easy strategy: free singlets (amares) ---
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_Cr_CH3,
    xxx.m_Cho_CH3,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock_arr = linklock_arr_unset.copy()
# link singlets by phase only
linklock_arr[metabolites_fit, :] = [0, 0, 0, 100]
linklock_arr[xxx.m_Cr_CH3, :] = [0, 0, 0, -100]
# leave water free
linklock_arr[xxx.m_Water, :] = [0, 0, 0, 0]
# leave Lipids linewidth free
linklock_arr[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for s in fit_strategies_seq_list:
    suffix = "_" + str(s)
    this_strategy = fit.fit_stategy("singlets_free" + suffix, metabolites_fit, linklock_arr, s)
    fit_stategies_list.append(this_strategy)

# --- medium strategy: singlets, CH2s, mI ---
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_mI,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock_arr = linklock_arr_unset.copy()
# link metabolites by linewidth and phase
linklock_arr[metabolites_fit, :] = [0, 200, 0, 100]
linklock_arr[xxx.m_Cr_CH3, :] = [0, -200, 0, -100]
# leave water free
linklock_arr[xxx.m_Water, :] = [0, 0, 0, 0]
# leave Lipids linewidth free
linklock_arr[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for s in fit_strategies_seq_list:
    suffix = "_" + str(s)
    this_strategy = fit.fit_stategy("singlets_CH2s_mI" + suffix, metabolites_fit, linklock_arr, s)
    fit_stategies_list.append(this_strategy)

# --- medium strategy: free singlets, CH2s, mI ---
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_mI,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock_arr = linklock_arr_unset.copy()
# link metabolites by linewidth and phase
linklock_arr[metabolites_fit, :] = [0, 200, 0, 100]
linklock_arr[xxx.m_mI, :] = [0, -200, 0, -100]
# leave water free
linklock_arr[xxx.m_Water, :] = [0, 0, 0, 0]
# let the singlets free
linklock_arr[xxx.m_NAA_CH3, :] = [0, 0, 0, 100]
linklock_arr[xxx.m_Cr_CH3, :] = [0, 0, 0, 100]
linklock_arr[xxx.m_Cho_CH3, :] = [0, 0, 0, 100]
# leave Lipids linewidth free
linklock_arr[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for s in fit_strategies_seq_list:
    suffix = "_" + str(s)
    this_strategy = fit.fit_stategy("free_singlets_CH2s_mI" + suffix, metabolites_fit, linklock_arr, s)
    fit_stategies_list.append(this_strategy)

# --- difficult strategy: singlets, CH2s, mI, Glx, Tau ---
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_Gln,
    xxx.m_Glu,
    xxx.m_mI,
    xxx.m_Tau,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock_arr = linklock_arr_unset.copy()
# link metabolites by linewidth and phase
linklock_arr[metabolites_fit, :] = [0, 200, 0, 100]
linklock_arr[xxx.m_Cr_CH3, :] = [0, -200, 0, -100]
# leave water free
linklock_arr[xxx.m_Water, :] = [0, 0, 0, 0]
# leave Lipids linewidth free
linklock_arr[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for s in fit_strategies_seq_list:
    suffix = "_" + str(s)
    this_strategy = fit.fit_stategy("singlets_CH2s_mI_Glx_Tau" + suffix, metabolites_fit, linklock_arr, s)
    fit_stategies_list.append(this_strategy)

# --- difficult strategy: free singlets, CH2s, mI, Glx, Tau ---
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_Gln,
    xxx.m_Glu,
    xxx.m_mI,
    xxx.m_Tau,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock_arr = linklock_arr_unset.copy()
# link metabolites by linewidth and phase
linklock_arr[metabolites_fit, :] = [0, 200, 0, 100]
linklock_arr[xxx.m_mI, :] = [0, -200, 0, -100]
# leave water free
linklock_arr[xxx.m_Water, :] = [0, 0, 0, 0]
# let the singlets free
linklock_arr[xxx.m_NAA_CH3, :] = [0, 0, 0, 100]
linklock_arr[xxx.m_Cr_CH3, :] = [0, 0, 0, 100]
linklock_arr[xxx.m_Cho_CH3, :] = [0, 0, 0, 100]
# leave Lipids linewidth free
linklock_arr[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for s in fit_strategies_seq_list:
    suffix = "_" + str(s)
    this_strategy = fit.fit_stategy("free_singlets_CH2s_mI_Glx_Tau" + suffix, metabolites_fit, linklock_arr, s)
    fit_stategies_list.append(this_strategy)

# --- extra difficult strategy: free singlets, CH2s, mI, Glx, Tau ---
metabolites_fit = np.sort([
    xxx.m_Asp,
    xxx.m_GABA,
    xxx.m_Glc,
    xxx.m_Gsh,
    xxx.m_GPC,
    xxx.m_Asc,
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_Gln,
    xxx.m_Glu,
    xxx.m_mI,
    xxx.m_Tau,
    xxx.m_Water,
    xxx.m_Lip1,
    xxx.m_Lip2])

# linklock: relations between fit parameters
linklock_arr = linklock_arr_unset.copy()
# link metabolites by linewidth and phase
linklock_arr[metabolites_fit, :] = [0, 200, 0, 100]
linklock_arr[xxx.m_mI, :] = [0, -200, 0, -100]
# leave water free
linklock_arr[xxx.m_Water, :] = [0, 0, 0, 0]
# let the singlets free
linklock_arr[xxx.m_NAA_CH3, :] = [0, 0, 0, 100]
linklock_arr[xxx.m_Cr_CH3, :] = [0, 0, 0, 100]
linklock_arr[xxx.m_Cho_CH3, :] = [0, 0, 0, 100]
# leave Lipids linewidth free
linklock_arr[[xxx.m_Lip1, xxx.m_Lip2], :] = [0, 0, 0, 100]

# store this strategy
for s in fit_strategies_seq_list:
    suffix = "_" + str(s)
    this_strategy = fit.fit_stategy("free_singlets_hardest" + suffix, metabolites_fit, linklock_arr, s)
    fit_stategies_list.append(this_strategy)

water_concentration = 55000.0  # mmol/kg

if(only_one_fit_strategy):
    fit_stategies_list = [fit_stategies_list[-1]]

# %% prepare data to fit

# browse though datasets
for this_hash, this_dataset in zip(df_sel.index, df_sel["dataset"]):
    # get the data
    this_data = this_dataset["raw"]["data"]
    if(this_data is None):
        # no raw data? ok get the dicom
        this_data = this_dataset["dcm"]["data"]

    # display the data
    if(display_stuff):
        this_data.display_spectrum_1d()

    # removing water and any artefact > 5ppm (corrupts fit quality criteria)
    this_data = this_data.correct_water_removal_1d(8, [4.5, 6], display_stuff)

    # filtering
    #this_data = this_data.correct_bandpass_filtering_1d([0, 6], np.ones, display_stuff)

    # reapodize to remove filtering artefact (will not affect linewidth because already apodized)
    #this_data = this_data.correct_apodization_nd(5.0, display_stuff)

    # %% prepare simulation machine

    # metabolite db
    meta_bs = sim.metabolite_basis_set()
    meta_bs.initialize()

    # sequence

    # take sequence from dataset, usually sLASER
    seq = this_data.sequence

    # %% area integration

    # for non water-suppressed data
    prefittool = fit.prefit_tool(this_data.data_ref, seq)
    prefittool.area_integration_peak_ranges = [0.5]
    prefittool.area_integration_peaks = [xxx.m_Water]
    prefittool.display_enable = display_stuff
    prefittool.initialize()
    _, params_ref_area_pnorm = prefittool.run()

    # for water-suppressed data
    prefittool = fit.prefit_tool(this_data, seq)
    prefittool.area_integration_peak_ranges = [0.1, 0.1, 0.1]
    prefittool.area_integration_peaks = [xxx.m_NAA_CH3, xxx.m_Cr_CH3, xxx.m_Cho_CH3]
    prefittool.display_enable = display_stuff
    prefittool.initialize()
    _, params_area_pnorm = prefittool.run()

    # %% fit non water-suppressed data

    fittool = fit.fit_tool(this_data.data_ref, seq)
    # min/max fitting bounds
    fittool.params_min = fittool.params_min.set_default_water_min()
    fittool.params_max = fittool.params_max.set_default_water_max()
    # water concentration parameter
    fittool.params_min[xxx.m_Water, xxx.p_cm] = 0
    fittool.params_max[xxx.m_Water, xxx.p_cm] = 1000000.0
    # water linewidth
    fittool.params_min[xxx.m_Water, xxx.p_dd] = 10.0
    fittool.params_max[xxx.m_Water, xxx.p_dd] = 500.0
    # water frequency shift
    fittool.params_min[xxx.m_Water, xxx.p_df] = -50.0
    fittool.params_max[xxx.m_Water, xxx.p_df] = +50.0
    # water phase
    fittool.params_min[xxx.m_Water, xxx.p_dp] = -0.1
    fittool.params_max[xxx.m_Water, xxx.p_dp] = +0.1

    # initial fit values
    fittool.params_init = fittool.params_min.copy()
    fittool.params_init[xxx.m_Water, xxx.p_cm] = 1.0
    fittool.params_init[xxx.m_Water, xxx.p_dd] = 100.0
    fittool.params_init[xxx.m_Water, xxx.p_df] = 0.0
    fittool.params_init[xxx.m_Water, xxx.p_dp] = 0.0

    # numerical optimization parameters
    fittool.optim_xtol = 1e-9
    fittool.optim_ftol = 1e-9
    fittool.optim_gtol = 1e-9
    fittool.display_range_ppm = [3, 7]

    fittool.display_enable = display_stuff

    # run the fit
    fittool.initialize()
    [params_ref_fit, optim_result_ref] = fittool.run()

    # %% use various fit strategies
    for this_fit_strategy in fit_stategies_list:

        # --- sequence ---
        if(this_fit_strategy.sequence == sim.mrs_sequence):
            # take sequence from dataset, usually sLASER
            seq = this_data.sequence
        else:
            # custom sequence
            seq = this_fit_strategy.sequence(seq.te, seq.tr, seq.na, seq.ds, seq.nuclei, seq.npts, seq.fs, seq.f0)

        # init
        seq.initialize(meta_bs)
        metabolites_fit = this_fit_strategy.metabolites

        # --- fit water-suppressed data ---

        fittool = fit.fit_tool(this_data, seq)

        # default fitting bounds from muscle template
        fittool.params_min = fittool.params_min.set_default_min().add_macromolecules_min()
        fittool.params_max = fittool.params_max.set_default_max().add_macromolecules_max()

        # re-init params_init
        fittool.params_init = (fittool.params_min + fittool.params_max) / 2.0

        # ranges for concentration
        fittool.params_min[:, xxx.p_cm] = 0.09
        fittool.params_max[:, xxx.p_cm] = 300.0
        fittool.params_max[xxx.m_All_MMs, xxx.p_cm] = 5000.0
        fittool.params_max[xxx.m_Water, xxx.p_cm] = 100000.0

        # linewidth bounds for metabolites
        # use water linewidth for max
        fittool.params_min[:, xxx.p_dd] = 5
        fittool.params_max[:, xxx.p_dd] = params_ref_fit[xxx.m_Water, xxx.p_dd] + 5
        # let the MMs go above
        fittool.params_max[xxx.m_All_MMs, xxx.p_dd] = 400
        # with early minimal damping
        fittool.params_init[:, xxx.p_dd] = fittool.params_min[:, xxx.p_dd] * 1.1

        # initial concentrations
        fittool.params_init[:, xxx.p_cm] = 0.0
        fittool.params_init[metabolites_fit, xxx.p_cm] = 0.1

        # frequency shifts for metabolites
        # TODO: understand why I have a 0.05ppm systematic shift
        fittool.params_init[:, xxx.p_df] = -15.0
        fittool.params_min[:, xxx.p_df] = -35.0
        fittool.params_max[:, xxx.p_df] = 5.0
        # a bit more for the lipids
        fittool.params_min[:, xxx.p_df] = -55.0
        fittool.params_max[:, xxx.p_df] = 25.0

        # phase bounds for metabolites and lipids
        fittool.params_min[:, xxx.p_dp] = -0.1
        fittool.params_max[:, xxx.p_dp] = +0.1

        # linklock
        linklock_arr = this_fit_strategy.linklock
        fittool.params_init.linklock[:] = linklock_arr.copy()

        # numerical optimization parameters
        fittool.display_range_ppm = [0, 6]
        fittool.display_frequency = 2

        fittool.display_enable = display_stuff

        # run the fit
        fittool.initialize()
        [params_fit_final, optim_result] = fittool.run()

        # --- save the fit pipeline to the db ---
        fit_results = {"params_ref_area_pnorm": params_ref_area_pnorm,
                        "params_area_pnorm": params_area_pnorm,
                        "params_ref_fit": params_ref_fit,
                        "params_fit_final": params_fit_final,
                        "optim_result": optim_result}

        # save in memory only
        rdb.save_fit_results(this_hash, fit_results, fittool, this_fit_strategy, False)

# now write to pkl file
rdb.write_pickle_files()
