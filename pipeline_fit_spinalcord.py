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
import warnings
import pickle

# constants
fit_cre_concentration = 7.5  # mmol/kg

# init
fit_metabolites_prefit = np.sort(fit_metabolites_prefit)
data = s
data_ref = s_ref

# %% prepare virtual sequence
seq = data.sequence
meta_bs = sim.metabolite_basis_set()
seq.initialize(meta_bs)

# %% add metabolites according to threshold
if(fit_metabolite_concentration_threshold is not None):
    fit_metabolites_threshold = sim.params(meta_bs).get_MBs_above_concentration(fit_metabolite_concentration_threshold)
    fit_metabolites = np.unique(np.concatenate((np.array(fit_metabolites), fit_metabolites_threshold)))

# %% prefit data
prefittool = fit.prefit_pipeline(data, seq)
prefittool.display_fig_index = 500
prefittool.area_integration_peakwidth = 25.0
prefittool.area_integration_peak_search_range = 0.1
prefittool.area_integration_peaks = fit_metabolites_prefit
params_prefit = prefittool.run()

# %% prefit ref. data if any
if(data_ref is not None):
    prefittool_ref = fit.prefit_pipeline(data_ref, seq)
    prefittool.display_fig_index = 600
    prefittool_ref.area_integration_peakwidth = 50.0
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
    fittool.initialize()
    [params_ref_fit, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

    # --- 2nd fit ---
    # adapt max bounds for concentration
    fittool.params_max[xxx.m_Water, xxx.p_cm] = params_ref_fit[xxx.m_Water, xxx.p_cm] * 10.0
    # init
    fittool.params_init = (fittool.params_min + fittool.params_max) / 2.0
    fittool.params_init[xxx.m_Water, xxx.p_dd] = fittool.params_min[xxx.m_Water, xxx.p_dd] * 1.1
    # go
    fittool.initialize()
    [params_ref_fit, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()
else:
    params_ref_fit = params_prefit_abs.copy()

# use prefit results to get initial parameters
params_init_fit = params_prefit_abs.correct_relative(params_ref_fit)
if(data_ref is None):
    params_init_fit[fit_metabolites_prefit, xxx.p_cm] = 0.1

# %% fit the singlets
fittool = fit.fit_pipeline(data, seq)

# --- 1st fit ---
# min bound
fittool.params_min.set_default_min()
fittool.params_min[:, xxx.p_cm] = 0.0
# max bound
fittool.params_max.set_default_max()
fittool.params_max[:, xxx.p_cm] = 60.0
# init
fittool.params_init = (fittool.params_min + fittool.params_max) / 2.0
fittool.params_init[fit_metabolites_prefit, xxx.p_dd] = fittool.params_min[fit_metabolites_prefit, xxx.p_dd] * 1.1
# copy initial concentrations found by area integration
fittool.params_init[:, xxx.p_cm] = 0.0
fittool.params_init[fit_metabolites_prefit, xxx.p_cm] = params_init_fit[fit_metabolites_prefit, xxx.p_cm]

# linklock
fittool.params_init.linklock[:] = 1
fittool.params_init.linklock[fit_metabolites_prefit[0], :] = [0, -2, -3, -4]
for im in fit_metabolites_prefit[1:]:
    fittool.params_init.linklock[im, :] = [0, 2, 3, 4]

# set normalization parameters
fittool.display_normalized = False
fittool.display_normalized_abs_ref_params = params_ref_fit.copy()

# go
fittool.initialize()
[params_fit_singlets, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

# --- 2nd fit ---
# adapt max bounds for concentration
fittool.params_max[fit_metabolites_prefit, xxx.p_cm] = params_fit_singlets[fit_metabolites_prefit, xxx.p_cm] * 10.0
# init
fittool.params_init[fit_metabolites_prefit, :] = (fittool.params_min[fit_metabolites_prefit, :] + fittool.params_max[fit_metabolites_prefit, :]) / 2.0
fittool.params_init[fit_metabolites_prefit, xxx.p_cm] = fittool.params_min[fit_metabolites_prefit, xxx.p_cm] + 0.1 * (fittool.params_max[fit_metabolites_prefit, xxx.p_cm] - fittool.params_min[fit_metabolites_prefit, xxx.p_cm])
# go
fittool.initialize()
[params_fit_singlets, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

# %% fit all metabolites

# --- 1st fit ---
# init: copy initial concentrations found with last fit
fittool.params_init[:, xxx.p_cm] = 0.0
fittool.params_init[fit_metabolites, xxx.p_cm] = 0.5
fittool.params_init[fit_metabolites_prefit, :] = params_fit_singlets[fit_metabolites_prefit, :]

# linklock
fittool.params_init.linklock[:] = 1
fittool.params_init.linklock[fit_metabolites[0], :] = [0, -2, -3, -4]
for im in fit_metabolites[1:]:
    fittool.params_init.linklock[im, :] = [0, 2, 3, 4]

# link Cre s/m, Cho s/m, NAA s/m
fittool.params_init.linklock[xxx.m_Cho_CH3, xxx.p_cm] = -5
fittool.params_init.linklock[xxx.m_Cho_CH2, xxx.p_cm] = +5
fittool.params_init.linklock[xxx.m_Cr_CH3, xxx.p_cm] = -6
fittool.params_init.linklock[xxx.m_Cr_CH2, xxx.p_cm] = +6
fittool.params_init.linklock[xxx.m_NAA_CH3, xxx.p_cm] = -7
fittool.params_init.linklock[xxx.m_NAA_CH2, xxx.p_cm] = +7

# go
fittool.initialize()
[params_fit_all_freq_locked, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

# --- 2nd fit ---
# adapt max bounds for concentration
fittool.params_max[fit_metabolites, xxx.p_cm] = params_fit_all_freq_locked[fit_metabolites, xxx.p_cm] * 10.0
# init
fittool.params_init[fit_metabolites, :] = (fittool.params_min[fit_metabolites, :] + fittool.params_max[fit_metabolites, :]) / 2.0
fittool.params_init[fit_metabolites, xxx.p_cm] = fittool.params_min[fit_metabolites, xxx.p_cm] + 0.1 * (fittool.params_max[fit_metabolites, xxx.p_cm] - fittool.params_min[fit_metabolites, xxx.p_cm])
# go
fittool.initialize()
[params_fit_all_freq_locked, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

# %% fit again all metabolites letting their shifts free

# init: copy initial concentrations found with last fit
fittool.params_init = params_fit_all_freq_locked.copy()

# adapt min/max freq shift to init: make it easier for singlets
fittool.params_min[fit_metabolites, xxx.p_df] = fittool.params_init[fit_metabolites_prefit[0], xxx.p_df] - 5.0
fittool.params_max[fit_metabolites, xxx.p_df] = fittool.params_init[fit_metabolites_prefit[0], xxx.p_df] + 5.0
fittool.params_min[fit_metabolites_prefit, xxx.p_df] = fittool.params_init[fit_metabolites_prefit[0], xxx.p_df] - 10.0
fittool.params_max[fit_metabolites_prefit, xxx.p_df] = fittool.params_init[fit_metabolites_prefit[0], xxx.p_df] + 10.0

# linklock
fittool.params_init.linklock[:] = 1
fittool.params_init.linklock[fit_metabolites[0], :] = [0, 0, 0, -2]
for im in fit_metabolites[1:]:
    fittool.params_init.linklock[im, :] = [0, 0, 0, 2]

# link Cre s/m, Cho s/m, NAA s/m
fittool.params_init.linklock[xxx.m_Cho_CH3, xxx.p_cm] = -3
fittool.params_init.linklock[xxx.m_Cho_CH2, xxx.p_cm] = +3
fittool.params_init.linklock[xxx.m_Cr_CH3, xxx.p_cm] = -4
fittool.params_init.linklock[xxx.m_Cr_CH2, xxx.p_cm] = +4
fittool.params_init.linklock[xxx.m_NAA_CH3, xxx.p_cm] = -5
fittool.params_init.linklock[xxx.m_NAA_CH2, xxx.p_cm] = +5

# go, final fit
fittool.initialize()
[params_fit_all, params_CRBs_abs, params_CRBs_rel, optim_result] = fittool.run()

# %% normalize concentrations

# correct for T2 relaxation assuming some T2s found in literature
params_fit_all_linked_Ts = params_fit_all.correct_T2s(data.te).correct_T1s(data.tr)
params_ref_fit_Ts = params_ref_fit.correct_T2s(data_ref.te).correct_T1s(data_ref.tr)
# normalize concentrations to water
params_fit_all_linked_Ts_abs_water = params_fit_all_linked_Ts.correct_absolute(params_ref_fit_Ts)
# normalize concentrations to creatine
params_fit_all_linked_Ts_abs_cre = params_fit_all_linked_Ts._correct_normalize_concentrations(fit_cre_concentration / params_fit_all_linked_Ts[xxx.m_Cr_CH3, xxx.p_cm])
# ratios to cre
params_fit_all_linked_Ts_ratios = params_fit_all_linked_Ts._correct_normalize_concentrations(1.0 / params_fit_all_linked_Ts[xxx.m_Cr_CH3, xxx.p_cm])

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
                  [sim.params(meta_bs).set_default_human_brain_min(),
                   params_prefit_abs,
                   params_fit_all_linked_Ts_abs_water,
                   params_fit_all_linked_Ts_abs_cre,
                   sim.params(meta_bs).set_default_human_brain_max()],

                  [pZeroes,
                      pZeroes,
                      pZeroes,
                      pZeroes,
                      pZeroes],

                  ['Human lower bounds',
                      'Area integration Abs.',
                      'Fit Ts-cor rel. Water',
                      'Fit Ts-cor rel. Cre',
                      'Human upper bounds'], False, False, False, False, xxx.p_cm, fit_metabolites, 0.1)
plt.ylim([0, 20])

# save fig in SVG format
fname = data.patient_name
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.tight_layout()
plt.pause(1)
plt.savefig('/home/tangir/desktop/' + fname + '_fit.svg', format='svg')

# save variables
with open('/home/tangir/desktop/' + fname + '.pkl', 'wb') as f:
    pickle.dump([params_fit_all, params_fit_all_linked_Ts, params_fit_all_linked_Ts_abs_cre, params_fit_all_linked_Ts_abs_water, params_fit_all_linked_Ts_ratios, params_fit_singlets,params_prefit, params_prefit_abs, params_prefit_ref, params_ref_fit, params_ref_fit_Ts, fit_metabolites, fit_metabolites_prefit, fname, data.birthyear, data.sex], f)
