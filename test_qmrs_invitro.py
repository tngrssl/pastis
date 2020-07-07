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
import mrs.reco as reco
import mrs.log as log

import numpy as np

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = reco.data_db()

# %% retrieve data to process


def my_sc_data(d: reco.MRSData2, p: reco.pipeline):
    """Filter data to get from database."""
    r = ("q_mrs_tests_cre" in d.patient_name)
    return(r)


data_list, _ = rdb.get_datasets(my_sc_data)
data = data_list[4]
data = data.correct_apodization_nd(5.0, True)

# %% prepare simulation machine

# metabolite db
meta_bs = sim.metabolite_basis_set()
meta_bs.initialize()

# sequence
seq = data.sequence
seq.initialize(meta_bs)

# %% set metabolites to fit, water concentration

# some fitting parameters
metabolites_fit = np.sort([
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Water])

# calculate water concentration knowing voxel
vox_vol_mm3 = data.voxel_volume()  # mm
vox_vol_l = vox_vol_mm3 * 0.000001  # L
vox_vol_g = vox_vol_l * 1000.0
vox_vol_mmol = 1000.0 * vox_vol_g / 18.0
water_concentration = vox_vol_mmol / vox_vol_l

# adjust T2s
meta_bs["Cr_CH3"]["T2"] = 300.0
meta_bs["Cr_CH2"]["T2"] = 300.0
meta_bs["Water"]["T2"] = 35.0

# %% area integration stuff

# for non water-suppressed data
prefittool = fit.prefit_tool(data.data_ref, seq)
prefittool.area_integration_peak_ranges = [1]
prefittool.area_integration_peaks = [xxx.m_Water]
prefittool.initialize()
params_ref_area, params_ref_area_pnorm = prefittool.run()

# for water-suppressed data
prefittool = fit.prefit_tool(data, seq)
prefittool.area_integration_peak_ranges = [1, 1]
prefittool.area_integration_peaks = [xxx.m_Cr_CH3, xxx.m_Cr_CH2]
prefittool.initialize()
params_area, params_area_pnorm = prefittool.run()

params_ref_area_pnorm_T2cor = params_ref_area_pnorm.correct_T2s(data.data_ref.te)
params_area_pnorm_T2cor = params_area_pnorm.correct_T2s(data.te)
params_area_pnorm_T2cor_abs = params_area_pnorm_T2cor.correct_absolute(params_ref_area_pnorm_T2cor, water_concentration)
params_area_pnorm_T2cor_abs.print()

# %% fit non water-suppressed data

fittool = fit.fit_tool(data.data_ref, seq)
# min/max fitting bounds
fittool.params_min = fittool.params_min.set_default_water_min()
fittool.params_max = fittool.params_max.set_default_water_max()
# water concentration parameter
fittool.params_min[xxx.m_Water, xxx.p_cm] = 0
fittool.params_max[xxx.m_Water, xxx.p_cm] = 100000.0
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

# run the fit
fittool.initialize()
[params_ref_fit, params_ref_fit_CRBs_abs, params_ref_fit_CRBs_rel, optim_result_ref] = fittool.run()

# %% fit water-suppressed data

fittool = fit.fit_tool(data, seq)

# default fitting bounds from muscle template
fittool.params_min = fittool.params_min.set_default_min()
fittool.params_max = fittool.params_max.set_default_max()

# ranges for concentration
fittool.params_min[:, xxx.p_cm] = 0.0
fittool.params_max[:, xxx.p_cm] = 300.0
fittool.params_max[xxx.m_Water, xxx.p_cm] = 300.0

# linewidth bounds for metabolites
fittool.params_min[xxx.m_All_MBs, xxx.p_dd] = 5.0
fittool.params_max[xxx.m_All_MBs, xxx.p_dd] = 90.0

# initial concentrations
fittool.params_init[:, xxx.p_cm] = 0.0
fittool.params_init[metabolites_fit, xxx.p_cm] = 0.1

# frequency shifts for metabolites
fittool.params_min[xxx.m_All_MBs, xxx.p_df] = -10.0
fittool.params_max[xxx.m_All_MBs, xxx.p_df] = 60.0
fittool.params_init[xxx.m_All_MBs, xxx.p_df] = 50.0

# phase bounds for metabolites and lipids
fittool.params_min[:, xxx.p_dp] = -0.1
fittool.params_max[:, xxx.p_dp] = +0.1

# linklock
fittool.params_init.linklock[:] = 1
fittool.params_init.linklock[metabolites_fit[0], :] = [0, -200, 0, -100]
for im in metabolites_fit[1:]:
    fittool.params_init.linklock[im, :] = [0, 200, 0, 100]

# leave water phase free
fittool.params_init.linklock[xxx.m_Water, :] = [0, 0, 0, 0]

# link Cr CH3 to CH2 ?
fittool.params_init.linklock[xxx.m_Cr_CH3, xxx.p_cm] = -6
fittool.params_init.linklock[xxx.m_Cr_CH2, xxx.p_cm] = +6

# numerical optimization parameters
fittool.display_range_ppm = [0, 6]
fittool.display_frequency = 2

# run the fit
fittool.initialize()
[params_fit_final, params_fit_CRBs_abs, params_fit_CRBs_rel, optim_result] = fittool.run()

# %% normalize and display concentrations

# correct for T2 relaxation
params_fit_final_Ts = params_fit_final.correct_T2s(data.te)
params_ref_fit_Ts = params_ref_fit.correct_T2s(data.data_ref.te)

# normalize concentrations to water
params_fit_final_Ts_abs_water = params_fit_final_Ts.correct_absolute(params_ref_fit_Ts, water_concentration)

# nice display
fig = plt.figure(700)
fig.clf()
fig.canvas.set_window_title("Final results")

# fit plot
ax = plt.subplot(1, 2, 1)
fit.disp_fit(ax, data, params_fit_final, seq, True, True)

# now the bargraphs
ax = plt.subplot(2, 2, 2)
fit.disp_bargraph(ax, [params_fit_final_Ts_abs_water], [params_fit_final_Ts_abs_water * params_fit_CRBs_rel / 100.0], ['Fit T2-cor rel. Water'], False, False, False, False, xxx.p_cm, np.concatenate([metabolites_fit]), 0.2)
plt.ylim([0, 75])

ax = plt.subplot(2, 2, 4)
fit.disp_bargraph(ax, [params_area_pnorm_T2cor_abs], [params_area_pnorm_T2cor_abs *0.0], ['Area Ts-cor rel. Water'], False, False, False, False, xxx.p_cm, np.concatenate([metabolites_fit]), 0.2)
plt.ylim([0, 75])

# adjust subplots
fig.subplots_adjust(left=0.02, bottom=0.1, right=0.98, top=0.95, wspace=0.2, hspace=None)
