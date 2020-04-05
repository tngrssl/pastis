#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script to be called once some data was processed using the MRSData2 class and is is currently in memory.

@author: Tangi Roussel
"""
# %% init
import mrs.fit as fit
import mrs.sim as sim
import mrs.aliases as xxx
import numpy as np
import matplotlib.pylab as plt

# constants
fit_cre_concentration = 7.5  # mmol/kg

# init
fit_metabolites = np.sort(fit_metabolites)
data = s
data_ref = s_ref

# %% prepare virtual sequence
seq = data.sequence
meta_bs = sim.metabolite_basis_set()()
seq.initialize(meta_bs)

# %% prefit data
prefittool = fit.prefit_pipeline(data, seq)
prefittool.display_fig_index = 500
prefittool.area_integration_peakwidth = 60.0
prefittool.area_integration_peak_search_range = 0.1
prefittool.area_integration_peaks = fit_metabolites
params_prefit = prefittool.run()

# %% prefit ref. data if any
if(data_ref is not None):
    prefittool_ref = fit.prefit_pipeline(data_ref, seq)
    prefittool.display_fig_index = 600
    prefittool_ref.area_integration_peakwidth = 120.0
    prefittool_ref.area_integration_peaks = [xxx.m_Water]
    params_prefit_ref = prefittool_ref.run()
else:
    params_prefit_ref = params_prefit.copy()

# calculate prefit absolute concentration
params_prefit_abs = params_prefit.correct_absolute(params_prefit_ref)
params_prefit_abs.print(False, False)

# %% fit reference scan
if(data_ref is not None):

    fittool = fit.fit_pipeline(data_ref, seq)
    # display
    fittool.display_subplots_types[3] = fit.fit_plot_type.BARGRAPH_DP

    # --- 1st fit ---
    # bounds
    fittool.params_min = fittool.params_min.set_default_water_min()
    fittool.params_max = fittool.params_max.set_default_water_max()
    # usual init
    fittool.params_init = (fittool.params_min + fittool.params_max) / 2.0
    fittool.params_init[xxx.m_Water, xxx.p_dd] = fittool.params_min[xxx.m_Water, xxx.p_dd] * 1.1
    fittool.params_init[xxx.m_Water, xxx.p_cm] = 0.1
    # go
    [params_ref_fit, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

    # --- 2nd fit ---
    # adapt max bounds for concentration
    fittool.params_max[xxx.m_Water, xxx.p_cm] = params_ref_fit[xxx.m_Water, xxx.p_cm] * 10.0
    # init
    fittool.params_init = (fittool.params_min + fittool.params_max) / 2.0
    fittool.params_init[xxx.m_Water, xxx.p_dd] = fittool.params_min[xxx.m_Water, xxx.p_dd] * 1.1
    # go
    [params_ref_fit, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()
else:
    params_ref_fit = params_prefit_abs.copy()

# use prefit results to get initial parameters
params_init_fit = params_prefit_abs.correct_relative(params_ref_fit)
if(data_ref is None):
    params_init_fit[fit_metabolites, xxx.p_cm] = 0.1

# %% fit
fittool = fit.fit_pipeline(data, seq)

# --- 1st fit ---
# min bound
fittool.params_min.set_default_human_brain_min()
fittool.params_min[:, xxx.p_cm] = 0.0
fittool.params_min[:, xxx.p_df] = -30.0

# max bound
fittool.params_max.set_default_human_brain_max()
fittool.params_max[:, xxx.p_cm] = 60.0
fittool.params_max[:, xxx.p_df] = +30.0

# init
fittool.params_init = (fittool.params_min + fittool.params_max) / 2.0
fittool.params_init[fit_metabolites, xxx.p_dd] = fittool.params_min[fit_metabolites, xxx.p_dd] * 1.1
fittool.params_init[fit_metabolites, xxx.p_df] = 0.0

# copy initial concentrations found by area integration
fittool.params_init[:, xxx.p_cm] = 0.0
fittool.params_init[fit_metabolites, xxx.p_cm] = params_init_fit[fit_metabolites, xxx.p_cm]

# linklock
fittool.params_init.linklock[:] = 1
fittool.params_init.linklock[fit_metabolites[0], :] = [0, -2, 0, -3]
for im in fit_metabolites[1:]:
    fittool.params_init.linklock[im, :] = [0, 2, 0, 3]

# set normalization parameters
fittool.display_normalized = False
fittool.display_normalized_abs_ref_params = params_ref_fit.copy()

# go
[params_fit_all, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

# %% normalize concentrations

# correct for T2 relaxation assuming some T2s found in literature
params_fit_all_linked_Ts = params_fit_all.correct_T2s(data.te).correct_T1s(data.tr)
params_ref_fit_Ts = params_ref_fit.correct_T2s(data_ref.te).correct_T1s(data_ref.tr)
# normalize concentrations to water
params_fit_all_linked_Ts_abs_water = params_fit_all_linked_Ts.correct_absolute(params_ref_fit_Ts)

# nice display
fig = plt.figure(700)
fig.clf()
fig.canvas.set_window_title("Final results")

# fit plot
ax = plt.subplot(1, 2, 1)
sim.disp_fit(ax, data, params_fit_all, seq, True, True)

# now the bargraph
ax = plt.subplot(2, 2, 2)
pZeroes = params_prefit_abs * 0.0

sim.disp_bargraph(ax,
                  [params_prefit_abs,
                   params_fit_all_linked_Ts_abs_water,],

                  [pZeroes,
                      pZeroes],

                  ['Area integration Abs.',
                      'Fit Ts-cor rel. Water'], False, False, False, False, xxx.p_cm, fit_metabolites, 0.1)
plt.ylim([0, 80])
