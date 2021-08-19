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
db_filepath = "/home/tangir/crmbm/acq_db/lipids_exvivo_bruker_500.pkl"

# %% retrieve data to fit

df = pd.read_pickle(db_filepath)

# 23/04/2021 steam / histo tests on foie gras sample
# df = df.loc[df["reco_dataset_legend"].str.contains("T2 estimation 23/04/2021")]

# 23/04/2021 steam test on foie gras sample
df = df.loc[df["reco_dataset_legend"].str.contains("STEAM2")]

# test fitting simulation
# s._noise_level = 0.01
# df["reco_dataset_dcm_data_obj"] = [s]

print(df)

# %% fit each dataset using independant lipid component model (simple)

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

fit_results_list = []

for this_data in df["reco_dataset_raw_data_obj"].to_list():

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
    fittool.params_min[:, xxx.p_dd] = 5.0
    fittool.params_max[:, xxx.p_dd] = 300.0
    # initial damping
    fittool.params_init[:, xxx.p_dd] = 100

    # frequency shifts
    fittool.params_min[:, xxx.p_df] = -5.0
    fittool.params_max[:, xxx.p_df] = 5.0
    # water
    #fittool.params_min[xxx.m_Water, xxx.p_df] = 10.0
    #fittool.params_max[xxx.m_Water, xxx.p_df] = 20.0
    #fittool.params_init[xxx.m_Water, xxx.p_df] = +13

    # phase shifts
    fittool.params_min[:, xxx.p_dp] = -0.1
    fittool.params_max[:, xxx.p_dp] = +0.1

    # linklock: global phase, frequency and damping, independant amplitudes
    fittool.params_linklock[:] = 1
    fittool.params_linklock[metabolites_list, :] =  [xxx.ll_FREE, xxx.ll_SLAVE1, xxx.ll_FREE, xxx.ll_SLAVE2]
    fittool.params_linklock[xxx.m_LipA, :] =        [xxx.ll_FREE, xxx.ll_MASTER1, xxx.ll_FREE, xxx.ll_MASTER2]

    # leave water free
    fittool.params_linklock[xxx.m_Water, :] =       [xxx.ll_FREE, xxx.ll_FREE, xxx.ll_FREE, xxx.ll_FREE]

    # run the fit
    fittool.run()

    # --- save the fit results ---
    fittool_df = fittool.to_dataframe('fit_')
    fittool_df = fittool_df.set_index(["fit_fit_hash", "fit_data_file_hash"])
    fit_results_list.append(fittool_df)

# append all the results and store
df_fit_results = pd.concat(fit_results_list, axis=0)
db_filepath_noext, _ = os.path.splitext(db_filepath)
df_fit_results.to_pickle(db_filepath_noext + "_fit.pkl")

# %% fit each dataset using crazy lipid model

# follow templates Lipids from xls file
meta_bs = sim.metabolite_basis_set("Lipids")

metabolites_list = np.sort([
                    xxx.m_Water,
                    xxx.m_LipA1,
                    xxx.m_LipB1,
                    xxx.m_LipB2,
                    xxx.m_LipB3,
                    xxx.m_LipC1,
                    xxx.m_LipD1,
                    xxx.m_LipD2,
                    xxx.m_LipE1,
                    xxx.m_LipF1 ])

                    # xxx.m_LipG1,
                    # xxx.m_LipH1,
                    # xxx.m_LipI1,
                    # xxx.m_LipJ1])

fit_results_list = []

for this_data in df["reco_dataset_raw_data_obj"].to_list():

    fittool = fit.fit_pastis(data=this_data, meta_bs=meta_bs)

    fittool.name = "crazy model fit"
    fittool.metabolites = metabolites_list
    fittool.metabolites_area_integration = metabolites_list
    fittool.display_range_ppm = [0, 6]

    # default fitting bounds from muscle template
    fittool.params_min = fittool.params_min.set_default_min()
    fittool.params_max = fittool.params_max.set_default_max()

    # ranges for concentration
    fittool.params_min[:, xxx.p_cm] = -1
    fittool.params_min[metabolites_list, xxx.p_cm] = 0
    fittool.params_max[:, xxx.p_cm] = 1000.0
    fittool.params_max[xxx.m_Water, xxx.p_cm] = 10000.0
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

    # fixed phase exception: negative peaks
    fittool.params_init[xxx.m_LipB2, xxx.p_dp] = np.pi
    fittool.params_init[xxx.m_LipD2, xxx.p_dp] = np.pi

    # linklock: global damping, phase locked
    fittool.params_linklock[:] = xxx.ll_FIXED
    fittool.params_linklock[metabolites_list, :] =  [xxx.ll_FREE, xxx.ll_SLAVE1, xxx.ll_FREE, xxx.ll_FIXED]
    fittool.params_linklock[xxx.m_LipA1, :] =       [xxx.ll_FREE, xxx.ll_MASTER1, xxx.ll_FREE, xxx.ll_FIXED]

    # crazy amplitudes linking here
    ##
    fittool.params_linklock[xxx.m_LipA1, xxx.p_cm] = xxx.ll_MASTER2
    fittool.params_linklock[xxx.m_LipC1, xxx.p_cm] = xxx.ll_SLAVE2
    fittool.params_linklock[xxx.m_LipE1, xxx.p_cm] = xxx.ll_SLAVE2
    # fittool.params_linklock[xxx.m_LipG1, xxx.p_cm] = xxx.ll_SLAVE2
    # fittool.params_linklock[xxx.m_LipH1, xxx.p_cm] = xxx.ll_SLAVE2
    # fittool.params_linklock[xxx.m_LipI1, xxx.p_cm] = xxx.ll_SLAVE2
    ## (CL-4) is alone and free
    fittool.params_linklock[xxx.m_LipB1, xxx.p_cm] = xxx.ll_FREE
    ## 2ndb group
    fittool.params_linklock[xxx.m_LipB2, xxx.p_cm] = xxx.ll_MASTER3
    fittool.params_linklock[xxx.m_LipD1, xxx.p_cm] = xxx.ll_SLAVE3
    # fittool.params_linklock[xxx.m_LipJ1, xxx.p_cm] = xxx.ll_SLAVE3
    ## 2nmidb group
    fittool.params_linklock[xxx.m_LipB3, xxx.p_cm] = xxx.ll_SLAVE4
    fittool.params_linklock[xxx.m_LipD2, xxx.p_cm] = xxx.ll_SLAVE4
    fittool.params_linklock[xxx.m_LipF1, xxx.p_cm] = xxx.ll_MASTER4

    ## link frequency shifts for positive-negative brother peaks
    fittool.params_linklock[xxx.m_LipB1, xxx.p_df] = xxx.ll_MASTER5
    fittool.params_linklock[xxx.m_LipB2, xxx.p_df] = xxx.ll_SLAVE5
    fittool.params_linklock[xxx.m_LipB3, xxx.p_df] = xxx.ll_SLAVE5
    #
    fittool.params_linklock[xxx.m_LipD1, xxx.p_df] = xxx.ll_MASTER6
    fittool.params_linklock[xxx.m_LipD2, xxx.p_df] = xxx.ll_SLAVE6

    # leave water free
    fittool.params_linklock[xxx.m_Water, :] = [xxx.ll_FREE, xxx.ll_FREE, xxx.ll_FREE, xxx.ll_FREE]

    # run the fit
    fittool.run()

    # quick test (should be done in notebook)
    print("--- Lipid fit results ---")

    # normalize concentrations using a estimated concentration known to be 1
    param_fit_res = fittool.params_fit.copy()
    param_fit_res[:, xxx.p_cm] = param_fit_res[:, xxx.p_cm] / param_fit_res[xxx.m_LipA1, xxx.p_cm]

    # CL, ndb, nmidb from normalized concentrations
    print("CL = %.2f" % (param_fit_res[xxx.m_LipB1, xxx.p_cm] + 4.0))
    print("ndb = %.2f" % (param_fit_res[xxx.m_LipB2, xxx.p_cm] / 2.0))
    print("nmidb = %.2f" % (param_fit_res[xxx.m_LipB3, xxx.p_cm] / 2.0))

    print("--- Lipid fit results, T2-corrected ---")

    # normalize concentrations using a estimated concentration known to be 1
    param_fit_res_T2c = fittool.params_fit.correct_T2s(this_data.te)
    param_fit_res_T2c[:, xxx.p_cm] = param_fit_res_T2c[:, xxx.p_cm] / param_fit_res_T2c[xxx.m_LipA1, xxx.p_cm]

    # CL, ndb, nmidb from normalized concentrations
    print("CL = %.2f" % (param_fit_res_T2c[xxx.m_LipB1, xxx.p_cm] + 4.0))
    print("ndb = %.2f" % (param_fit_res_T2c[xxx.m_LipB2, xxx.p_cm] / 2.0))
    print("nmidb = %.2f" % (param_fit_res_T2c[xxx.m_LipB3, xxx.p_cm] / 2.0))

    # --- save the fit results ---
    fittool_df = fittool.to_dataframe('fit_')
    fittool_df = fittool_df.set_index(["fit_fit_hash", "fit_data_file_hash"])
    fit_results_list.append(fittool_df)

# append all the results and store
df_fit_results = pd.concat(fit_results_list, axis=0)
db_filepath_noext, _ = os.path.splitext(db_filepath)
df_fit_results.to_pickle(db_filepath_noext + "_fit.pkl")
