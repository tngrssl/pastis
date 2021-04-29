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

import pandas as pd
import numpy as np
import os

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.INFO)

# data to process is in here
db_filepath = "/home/tangir/crmbm/acq_db/lipids.pkl"

# %% retrieve data to fit

df = pd.read_pickle(db_filepath)
df = df.loc[df["reco_dataset_legend"].str.contains("T2 estimation 23/04/2021")]

# %% set metabolites to fit

# follow templates Lipids from xls file
meta_bs = sim.metabolite_basis_set("Lipids")

metabolites_list = np.sort([
                    xxx.m_Water,
                    xxx.m_LipA,
                    xxx.m_LipB,
                    xxx.m_LipC,
                    xxx.m_LipD,
                    xxx.m_LipE,
                    xxx.m_LipF,
                    xxx.m_LipG,
                    xxx.m_LipH,
                    xxx.m_LipI,
                    xxx.m_LipJ])

# %% fit each dataset

fit_results_list = []

for this_data in df["reco_dataset_dcm_data_obj"].to_list():

    fittool = fit.fit_pastis(data=this_data, meta_bs=meta_bs)

    fittool.name = "fitting T2w mrs"
    fittool.metabolites = metabolites_list
    fittool.metabolites_area_integration = metabolites_list
    fittool.display_range_ppm = [0, 6]
    # fittool.display_CRBs = False

    # default fitting bounds from muscle template
    fittool.params_min = fittool.params_min.set_default_min()
    fittool.params_max = fittool.params_max.set_default_max()

    # ranges for concentration
    fittool.params_min[:, xxx.p_cm] = -1
    fittool.params_min[metabolites_list, xxx.p_cm] = 0
    fittool.params_max[:, xxx.p_cm] = 5000.0
    # initial concentrations
    fittool.params_init[:, xxx.p_cm] = 0
    fittool.params_init[metabolites_list, xxx.p_cm] = 0.1

    # linewidth bounds for metabolites
    fittool.params_min[:, xxx.p_dd] = 10.0
    fittool.params_max[:, xxx.p_dd] = 300.0
    # initial damping
    fittool.params_init[:, xxx.p_dd] = 100

    # frequency shifts
    fittool.params_min[:, xxx.p_df] = -5.0
    fittool.params_max[:, xxx.p_df] = 5.0
    # water
    fittool.params_min[xxx.m_Water, xxx.p_df] = 10.0
    fittool.params_max[xxx.m_Water, xxx.p_df] = 20.0
    fittool.params_init[xxx.m_Water, xxx.p_df] = +13

    # phase shifts
    fittool.params_min[:, xxx.p_dp] = -0.1
    fittool.params_max[:, xxx.p_dp] = +0.1

    # linklock: global phase, frequency and damping, independant amplitudes
    fittool.params_linklock[:] = 1
    fittool.params_linklock[metabolites_list, :] = [0, 200, 0, 100]
    fittool.params_linklock[xxx.m_LipA, :] = [0, -200, 0, -100]

    # leave water free
    fittool.params_linklock[xxx.m_Water, :] = [0, 0, 0, 0]

    # run the fit
    fittool.run()

    # --- save the fit results ---
    fittool_df = fittool.to_dataframe('lip_')
    fittool_df = fittool_df.set_index(["lip_fit_hash", "lip_data_file_hash"])
    fit_results_list.append(fittool_df)

# append all the results and store
df_fit_results = pd.concat(fit_results_list, axis=0)
db_filepath_noext, _ = os.path.splitext(db_filepath)
df_fit_results.to_pickle(db_filepath_noext + "_fit.pkl")
