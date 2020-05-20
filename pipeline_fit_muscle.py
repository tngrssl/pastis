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
import csv
import os.path

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = reco.data_db()

# %% retrieve data to process
data, data_pipeline = rdb.get_latest_dataset()

# %% set metabolites to fit, water concentration, csv file

# some fitting parameters
metabolites_fit = np.sort([
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Alcar,
    xxx.m_Carni,
    xxx.m_Carno,
    xxx.m_Water])

lipids_fit = np.sort([
    xxx.m_Lip1,
    xxx.m_Lip2,
    xxx.m_Lip3,
    xxx.m_Lip4,
    xxx.m_Lip5,
    xxx.m_Lip6])

water_concentration = 87700.0  # mmol/kg
fit_result_csv_filename = "fit_results.csv"

# %% prepare simulation machine

# metabolite db
meta_bs = sim.metabolite_basis_set()
meta_bs.basis_set_xls_file = "./metabolite_basis_sets/muscle_7T.xls"
meta_bs.non_coupled_only = True
meta_bs.initialize()

# sequence
seq = data.sequence
seq.seqdb_file = "./metabolite_basis_sets/muscle_7T_20200519.pkl"
seq.initialize(meta_bs)

# B0 field switch
gamma1H = 42.576
b0 = seq.f0 / gamma1H
b0_factor = b0 / 7.0

# %% fit non water-suppressed data

fittool = fit.fit_tool(data.data_ref, seq)
# min/max fitting bounds
fittool.params_min = fittool.params_min.set_default_water_min()
fittool.params_max = fittool.params_max.set_default_water_max()
# water concentration parameter
fittool.params_min[xxx.m_Water, xxx.p_cm] = 0
fittool.params_max[xxx.m_Water, xxx.p_cm] = 100000.0
# water linewidth
fittool.params_min[xxx.m_Water, xxx.p_dd] = 10.0 * b0_factor
fittool.params_max[xxx.m_Water, xxx.p_dd] = 500.0 * b0_factor
# water frequency shift
fittool.params_min[xxx.m_Water, xxx.p_df] = -50.0 * b0_factor
fittool.params_max[xxx.m_Water, xxx.p_df] = +50.0 * b0_factor
# water phase
fittool.params_min[xxx.m_Water, xxx.p_dp] = -0.1
fittool.params_max[xxx.m_Water, xxx.p_dp] = +0.1

# initial fit values
fittool.params_init = fittool.params_min.copy()
fittool.params_init[xxx.m_Water, xxx.p_cm] = 1.0
fittool.params_init[xxx.m_Water, xxx.p_dd] = 100.0 * b0_factor
fittool.params_init[xxx.m_Water, xxx.p_df] = 0.0 * b0_factor
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
fittool.params_max[:, xxx.p_cm] = 1000.0
fittool.params_max[xxx.m_Water, xxx.p_cm] = 100000.0

# linewidth bounds for metabolites
fittool.params_min[xxx.m_All_MBs, xxx.p_dd] = 5.0 * b0_factor
fittool.params_max[xxx.m_All_MBs, xxx.p_dd] = 90.0 * b0_factor

# initial concentrations
fittool.params_init[:, xxx.p_cm] = 0.0
fittool.params_init[metabolites_fit, xxx.p_cm] = 0.1

# frequency shifts for metabolites
fittool.params_min[xxx.m_All_MBs, xxx.p_df] = -20.0 * b0_factor
fittool.params_max[xxx.m_All_MBs, xxx.p_df] = 20.0 * b0_factor

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
# fittool.params_init.linklock[xxx.m_Cr_CH3, xxx.p_cm] = -6
# fittool.params_init.linklock[xxx.m_Cr_CH2, xxx.p_cm] = +6

# lipids concentration init
fittool.params_init[lipids_fit, xxx.p_cm] = 0.1

# max lipid cm
fittool.params_max[lipids_fit, xxx.p_cm] = 10000.0

# linewidth bounds for metabolites
fittool.params_min[lipids_fit, xxx.p_dd] = 50.0 * b0_factor
fittool.params_max[lipids_fit, xxx.p_dd] = 225.0 * b0_factor
fittool.params_init[lipids_fit, xxx.p_dd] = 51.0 * b0_factor

# lipids frequency shift init and bound values
fittool.params_min[lipids_fit, xxx.p_df] = -40.0 * b0_factor
fittool.params_max[lipids_fit, xxx.p_df] = +40.0 * b0_factor
fittool.params_init[lipids_fit, xxx.p_df] = 0.0 * b0_factor

fittool.params_init[xxx.m_Lip3, xxx.p_df] = +12
fittool.params_init[xxx.m_Lip4, xxx.p_df] = +3
fittool.params_init[xxx.m_Lip6, xxx.p_df] = -9

# lipids linklock
fittool.params_init.linklock[xxx.m_All_MMs, :] = 1
fittool.params_init.linklock[lipids_fit[0], :] = [0, 200, 0, 100]
for il in lipids_fit[1:]:
    fittool.params_init.linklock[il, :] = [0, 200, 0, 100]

# numerical optimization parameters
fittool.display_range_ppm = [0, 6]
fittool.display_frequency = 2

# run the fit
fittool.initialize()
[params_fit_final, params_fit_CRBs_abs, params_fit_CRBs_rel, optim_result] = fittool.run()

# %% normalize and display concentrations

# correct for T2/T1 relaxation assuming some T1/T2s found in literature
params_fit_final_Ts = params_fit_final.correct_T2s(data.te).correct_T1s(data.tr)
params_ref_fit_Ts = params_ref_fit.correct_T2s(data.data_ref.te).correct_T1s(data.data_ref.tr)

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
fit.disp_bargraph(ax, [params_fit_final_Ts_abs_water], [params_fit_final_Ts_abs_water * params_fit_CRBs_rel / 100.0], ['Fit Ts-cor rel. Water'], False, False, False, False, xxx.p_cm, np.concatenate([metabolites_fit, lipids_fit]), 0.2)

ax = plt.subplot(2, 2, 4)
fit.disp_bargraph(ax, [params_fit_final_Ts_abs_water], [params_fit_final_Ts_abs_water * params_fit_CRBs_rel / 100.0], ['Fit Ts-cor rel. Water'], False, False, False, False, xxx.p_cm, np.concatenate([metabolites_fit, lipids_fit]), 0.2)
plt.ylim([0, 1000])

# adjust subplots
fig.subplots_adjust(left=0.02, bottom=0.1, right=0.98, top=0.95, wspace=0.2, hspace=None)

# %% save results to csv file

# prepare stuff to write in csv file
metabolites_names = np.array(params_fit_final.get_meta_names())
data_legend_caption = data_pipeline._display_legends_list[0]
header_line = ["Dataset"] + list(metabolites_names[np.concatenate([metabolites_fit, lipids_fit])]) + [metabolites_names[xxx.m_Water] + "(REF)"] + [lbl +" (CRB%)" for lbl in list(metabolites_names[np.concatenate([metabolites_fit, lipids_fit])])] + [metabolites_names[xxx.m_Water] + "(REF) (CRB%)"]

data_line = [data_legend_caption] + list(params_fit_final[np.concatenate([metabolites_fit, lipids_fit]), xxx.p_cm]) + [params_ref_fit[xxx.m_Water, xxx.p_cm]] + list(params_fit_CRBs_rel[np.concatenate([metabolites_fit, lipids_fit]), xxx.p_cm]) + [params_ref_fit_CRBs_rel[xxx.m_Water, xxx.p_cm]]

# if csv file does not exist, create it and write header row
if not os.path.isfile(fit_result_csv_filename):
    with open(fit_result_csv_filename, 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        spamwriter.writerow(header_line)
        csvfile.close()

# write data
with open(fit_result_csv_filename, 'a') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    spamwriter.writerow(data_line)
    csvfile.close()

