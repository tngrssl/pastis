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

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.INFO)

# data to process is in here
db_filepath = "/home/tangir/crmbm/acq_db/muscle.pkl"

# %% retrieve data to process

df = pd.read_pickle(db_filepath)
df = df.loc[df["timestamp"] == df["timestamp"].max()]
data = df.iloc[0]["reco_dataset_raw_data_obj"]

# %% set metabolites to fit, water concentration, csv file

# some fitting parameters
metabolites_fit = np.sort([
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Alcar,
    xxx.m_Water])

data.data_ref = data.data_ref.correct_freqshift_1d([4, 6])

lipids_fit = np.sort([
    xxx.m_Lip1,
    xxx.m_Lip2,
    xxx.m_Lip3,
    xxx.m_Lip4,
    xxx.m_Lip5,
    xxx.m_Lip6])

water_concentration = 87700.0  # mmol/kg

# %% prepare simulation machine


basis_set_name = 'Muscle'
basis_set_pkl_file = "20210331_muscle.pkl"

# metabolite db
meta_bs = sim.metabolite_basis_set(basis_set_name=basis_set_name, non_coupled_only=True, one_proton_mode=True)

# sequence
seq = data.sequence
seq.db_file = "20210331_muscle.pkl"
seq.initialize(meta_bs)

# B0 field switch
gamma1H = 42.576
b0 = seq.f0 / gamma1H
b0_factor = b0 / 7.0

# %% fit non water-suppressed data

fittool_nows = fit.fit_tool(data.data_ref, seq)
# min/max fitting bounds
fittool_nows.params_min = fittool_nows.params_min.set_default_water_min()
fittool_nows.params_max = fittool_nows.params_max.set_default_water_max()
# water concentration parameter
fittool_nows.params_min[xxx.m_Water, xxx.p_cm] = 0
fittool_nows.params_max[xxx.m_Water, xxx.p_cm] = 1000000.0
# water linewidth
fittool_nows.params_min[xxx.m_Water, xxx.p_dd] = 10.0 * b0_factor
fittool_nows.params_max[xxx.m_Water, xxx.p_dd] = 500.0 * b0_factor
# water frequency shift
fittool_nows.params_min[xxx.m_Water, xxx.p_df] = -50.0 * b0_factor
fittool_nows.params_max[xxx.m_Water, xxx.p_df] = +50.0 * b0_factor
# water phase
fittool_nows.params_min[xxx.m_Water, xxx.p_dp] = -0.1
fittool_nows.params_max[xxx.m_Water, xxx.p_dp] = +0.1

# initial fit values
fittool_nows.params_init = fittool_nows.params_min.copy()
fittool_nows.params_init[xxx.m_Water, xxx.p_cm] = 1.0
fittool_nows.params_init[xxx.m_Water, xxx.p_dd] = 100.0 * b0_factor
fittool_nows.params_init[xxx.m_Water, xxx.p_df] = 0.0 * b0_factor
fittool_nows.params_init[xxx.m_Water, xxx.p_dp] = 0.0

# numerical optimization parameters
fittool_nows.optim_xtol = 1e-9
fittool_nows.optim_ftol = 1e-9
fittool_nows.optim_gtol = 1e-9
fittool_nows.display_range_ppm = [3, 7]

# run the fit
fittool_nows.run()

# %% fit water-suppressed data

fittool_ws = fit.fit_tool(data, seq)

# default fitting bounds from muscle template
fittool_ws.params_min = fittool_ws.params_min.set_default_min()
fittool_ws.params_max = fittool_ws.params_max.set_default_max()

# ranges for concentration
fittool_ws.params_min[:, xxx.p_cm] = 0.0
fittool_ws.params_max[:, xxx.p_cm] = 1000.0
fittool_ws.params_max[xxx.m_Water, xxx.p_cm] = 100000.0

# linewidth bounds for metabolites
fittool_ws.params_min[xxx.m_All_MBs, xxx.p_dd] = 5.0 * b0_factor
fittool_ws.params_max[xxx.m_All_MBs, xxx.p_dd] = 90.0 * b0_factor

# initial concentrations
fittool_ws.params_init[:, xxx.p_cm] = 0.0
fittool_ws.params_init[metabolites_fit, xxx.p_cm] = 0.1

# frequency shifts for metabolites
fittool_ws.params_min[xxx.m_All_MBs, xxx.p_df] = -20.0 * b0_factor
fittool_ws.params_max[xxx.m_All_MBs, xxx.p_df] = 20.0 * b0_factor

# phase bounds for metabolites and lipids
fittool_ws.params_min[:, xxx.p_dp] = -0.1
fittool_ws.params_max[:, xxx.p_dp] = +0.1

# linklock
fittool_ws.params_init.linklock[:] = 1
fittool_ws.params_init.linklock[metabolites_fit[0], :] = [0, -200, 0, -100]
for im in metabolites_fit[1:]:
    fittool_ws.params_init.linklock[im, :] = [0, 200, 0, 100]

# leave water phase free
fittool_ws.params_init.linklock[xxx.m_Water, :] = [0, 0, 0, 0]

# link Cr CH3 to CH2 ?
# fittool_ws.params_init.linklock[xxx.m_Cr_CH3, xxx.p_cm] = -6
# fittool_ws.params_init.linklock[xxx.m_Cr_CH2, xxx.p_cm] = +6

# lipids concentration init
fittool_ws.params_init[lipids_fit, xxx.p_cm] = 0.1

# max lipid cm
fittool_ws.params_max[lipids_fit, xxx.p_cm] = 1000.0

# linewidth bounds for metabolites
fittool_ws.params_min[lipids_fit, xxx.p_dd] = 70.0 * b0_factor
fittool_ws.params_max[lipids_fit, xxx.p_dd] = 225.0 * b0_factor
fittool_ws.params_init[lipids_fit, xxx.p_dd] = 71.0 * b0_factor

# lipids frequency shift init and bound values
fittool_ws.params_min[lipids_fit, xxx.p_df] = -40.0 * b0_factor
fittool_ws.params_max[lipids_fit, xxx.p_df] = +60.0 * b0_factor
fittool_ws.params_init[lipids_fit, xxx.p_df] = 0.0 * b0_factor

fittool_ws.params_init[xxx.m_Lip3, xxx.p_df] = +12
fittool_ws.params_init[xxx.m_Lip4, xxx.p_df] = +3
fittool_ws.params_init[xxx.m_Lip6, xxx.p_df] = -9

# lipids linklock
fittool_ws.params_init.linklock[xxx.m_All_MMs, :] = 1
fittool_ws.params_init.linklock[lipids_fit[0], :] = [0, 0, 0, 100]
for il in lipids_fit[1:]:
    fittool_ws.params_init.linklock[il, :] = [0, 0, 0, 100]

# add a "baseline" using vey large MM resonances
fittool_ws.params_min.add_macromolecules_min()
fittool_ws.params_max.add_macromolecules_max()
# take ony one resonance
fittool_ws.params_init.linklock[xxx.m_All_MMs, :] = [1, 1, 1, 1]
fittool_ws.params_init.linklock[xxx.m_Offset2, :] = [0, 0, 0, 0]

# huge linewidth
fittool_ws.params_min[xxx.m_Offset2, xxx.p_dd] = 300
fittool_ws.params_max[xxx.m_Offset2, xxx.p_dd] = 2000
fittool_ws.params_init[xxx.m_Offset2, xxx.p_dd] = 600

# not too strong
fittool_ws.params_max[xxx.m_Offset2, xxx.p_cm] = 200
fittool_ws.params_init[xxx.m_Offset2, xxx.p_cm] = 0.1

# numerical optimization parameters
fittool_ws.display_range_ppm = [0, 6]
fittool_ws.display_frequency = 2

# run the fit
fittool_ws.run()

# %% dump results to csv file

# save the fit results ---
this_fit_nows_df = fittool_nows.to_dataframe('nows_')
this_fit_ws_df = fittool_ws.to_dataframe('ws_')
this_fit_df = pd.concat([this_fit_nows_df, this_fit_ws_df], axis=1)
this_fit_df = this_fit_df.set_index(["ws_fit_hash", "ws_data_file_hash"])
this_fit_df.to_csv("fit_results.csv")
