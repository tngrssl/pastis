#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script to be called once some data was processed using the MRSData2 class and is is currently in memory.

@author: Tangi Roussel
"""

from __future__ import division

import mrs.fit as fit
import mrs.metabase as xxx
import mrs.sim as sim
import numpy as np
import matplotlib.pylab as plt
import warnings
import csv
import os.path

# init
warnings.filterwarnings("ignore", ".*GUI is implemented*")
metabolites_fit = np.sort(metabolites_fit)
lipids_fit = np.sort(lipids_fit)
data = s
data_ref = s_ref

# %% prepare simulation machine

simtool = sim.metabolite_database()
simtool.mega_meta_db_file =  "C:/Users/jsourdon/Documents/Software/Tangi/metabolite_db/20200227.pkl"  
simtool.mrs_seq = sim.mrs_sequence.SPINECHO  # my STEAM is not working for simulation yet
simtool.pg_Fref = data.f0
simtool.te = data.te / 1000.0
simtool.pg_Fs = 1.0 / data.dt
simtool.pg_nPts = data.shape[0]
simtool.pg_amp_factor = 2.0e-5
simtool.sigmaNoise = 0
simtool.initialize()

# B0 field switch
gamma1H = 42.576
b0 = data.f0 / gamma1H
b0_factor = b0 / 7.0

# %% fit non water-suppressed data

fittool = fit.fit_pipeline(data_ref, simtool)
fittool.optim_jacobian = False
# min/max fitting bounds
fittool.params_min = fittool.params_min.set_default_water_min()
fittool.params_max = fittool.params_max.set_default_water_max()
# water concentration parameter
fittool.params_min[xxx.m_Water, xxx.p_cm] = 1.0
fittool.params_max[xxx.m_Water, xxx.p_cm] = 10000.0
# water linewidth
fittool.params_min[xxx.m_Water, xxx.p_dd] = 30.0 * b0_factor
fittool.params_max[xxx.m_Water, xxx.p_dd] = 500.0 * b0_factor
# water frequency shift
fittool.params_min[xxx.m_Water, xxx.p_df] = -50.0 * b0_factor
fittool.params_max[xxx.m_Water, xxx.p_df] = +50.0 * b0_factor
# water phase
fittool.params_min[xxx.m_Water, xxx.p_dp] = -0.1
fittool.params_max[xxx.m_Water, xxx.p_dp] = +0.1

# initial fit values
fittool.params_init = fittool.params_min.copy()
fittool.params_init[xxx.m_Water, xxx.p_cm] = 100.0
fittool.params_init[xxx.m_Water, xxx.p_dd] = 100.0 * b0_factor
fittool.params_init[xxx.m_Water, xxx.p_df] = 0.0 * b0_factor
fittool.params_init[xxx.m_Water, xxx.p_dp] = 0.0

# numerical optimization parameters
fittool.optim_ppm_range = [3, 7]
fittool.optim_xtol = 1e-15
fittool.optim_ftol = 1e-15
fittool.optim_gtol = 1e-15

# run the fit
[params_ref_fit, optim_results] = fittool.run()

# %% fit water-suppressed data

fit_stategy_list = []
fit_stategy_list.append([[xxx.m_Cr_CH3],[]])
fit_stategy_list.append([metabolites_fit,lipids_fit])

for k, this_fit_strategy in enumerate(fit_stategy_list):
    
    this_metabolites_fit = this_fit_strategy[0]
    this_lipids_fit = this_fit_strategy[1]
    
    fittool = fit.fit_pipeline(data, simtool)
    fittool.optim_jacobian = False
    
    # default fitting bounds for brain
    fittool.params_min = fittool.params_min.set_default_human_min()
    fittool.params_max = fittool.params_max.set_default_human_max()
    
    # adapt for muscle: turn on bounds for muscle metabolites
    fittool.params_min[:, xxx.p_cm] = 0.0
    fittool.params_max[:, xxx.p_cm] = 10.0
    
    # linewidth bounds for metabolites
    fittool.params_min[xxx.m_All_MBs, xxx.p_dd] = 5.0 * b0_factor
    fittool.params_max[xxx.m_All_MBs, xxx.p_dd] = 75.0 * b0_factor
    
    # initial
    fittool.params_init[:, xxx.p_cm] = 0.0
    if(k == 0):
        fittool.params_init[xxx.m_All_MBs, xxx.p_dd] = fittool.params_min[xxx.m_All_MBs, xxx.p_dd]
        fittool.params_init[:, xxx.p_dp] = 0.0
        fittool.optim_ppm_range = [2.75, 3.25]
    else:
        fittool.params_init[this_metabolites_fit, xxx.p_dd] = params_fit_global_fshift[this_last_metabolites_fit, xxx.p_dd]
        fittool.params_init[this_metabolites_fit, xxx.p_df] = params_fit_global_fshift[this_last_metabolites_fit, xxx.p_df]
        fittool.optim_ppm_range = [0, 4]
   
    fittool.params_init[this_metabolites_fit, xxx.p_cm] = 0.0
    
    # frequency shifts for metabolites
    fittool.params_min[xxx.m_All_MBs, xxx.p_df] = -10.0 * b0_factor
    fittool.params_max[xxx.m_All_MBs, xxx.p_df] = 10.0 * b0_factor
    
    # phase bounds for metabolites and lipids
    fittool.params_min[:, xxx.p_dp] = -0.1
    fittool.params_max[:, xxx.p_dp] = +0.1
    
    # linklock
    fittool.params_init.linklock[:] = 1
    if(len(this_metabolites_fit)==1):
        fittool.params_init.linklock[this_metabolites_fit, :] = 0
    else:
        fittool.params_init.linklock[this_metabolites_fit[0], :] = [0, -2, -3, -4]
        for im in this_metabolites_fit[1:]:
            fittool.params_init.linklock[im, :] = [0, 2, 3, 4]
    
    # link Cr CH3 to CH2
    # fittool.params_init.linklock[xxx.m_Cr_CH3, xxx.p_cm] = -6
    # fittool.params_init.linklock[xxx.m_Cr_CH2, xxx.p_cm] = +6
    
    if(len(this_lipids_fit) > 0):
        # add lipids
        fittool.params_min.add_macromolecules_min()
        fittool.params_max.add_macromolecules_max()
        
        # lipids concentration init
        fittool.params_init[this_lipids_fit, xxx.p_cm] = 0.1
        
        # linewidth bounds for metabolites
        fittool.params_min[this_lipids_fit, xxx.p_dd] = 50.0 * b0_factor
        fittool.params_max[this_lipids_fit, xxx.p_dd] = 250.0 * b0_factor
        fittool.params_init[this_lipids_fit, xxx.p_dd] = 100.0 * b0_factor
        
        # lipids frequency shift init and bound values
        fittool.params_min[this_lipids_fit, xxx.p_df] = -70.0 * b0_factor
        fittool.params_max[this_lipids_fit, xxx.p_df] = +70.0 * b0_factor
        fittool.params_init[this_lipids_fit, xxx.p_df] = 0.0 * b0_factor
        
        # lipids linklock
        fittool.params_init.linklock[xxx.m_All_MMs, :] = 1
        fittool.params_init.linklock[this_lipids_fit, :] = [0, 200, 100, 4]
        fittool.params_init.linklock[this_lipids_fit[0], :] = [0, -200, -100, 4]
        
        # except Lip4 that needs some special tweak
        fittool.params_min[xxx.m_Lip4, xxx.p_df] = 75.0 * b0_factor
        fittool.params_max[xxx.m_Lip4, xxx.p_df] = 150.0 * b0_factor
        fittool.params_init[xxx.m_Lip4, xxx.p_df] = 100.0 * b0_factor
        fittool.params_init.linklock[xxx.m_Lip4, :] = [0, 200, 0, 4]
    
    # numerical optimization parameters
    fittool.display_range_ppm = [0, 6]
    fittool.optim_xtol = 1e-13
    fittool.optim_ftol = 1e-13
    fittool.optim_gtol = 1e-13
    
    # run the fit
    [params_fit_global_fshift, optim_results] = fittool.run()
    
    # save last strategy
    this_last_metabolites_fit = this_fit_strategy[0]
    this_last_lipids_fit = this_fit_strategy[1]

# %% cool, now let's fit again all metabolites letting them free

# init: copy initial concentrations found with last fit
fittool.params_init = params_fit_global_fshift.copy()

# linklock
fittool.params_init.linklock[metabolites_fit, xxx.p_df] = 0
fittool.params_init.linklock[lipids_fit, xxx.p_df] = 0
# fittool.params_init.linklock[metabolites_fit, xxx.p_dd] = 0
fittool.params_init.linklock[lipids_fit, xxx.p_dd] = 0

# go, final fit
[params_fit_final, optim_results] = fittool.run()

# %% normalize concentrations

# correct for T2/T1 relaxation assuming some T1/T2s found in literature
params_fit_final_Ts = params_fit_final.correct_T2s(data.te).correct_T1s(data.tr)
params_ref_fit_Ts = params_ref_fit.correct_T2s(data_ref.te).correct_T1s(data_ref.tr)
# normalize concentrations to water
params_fit_final_Ts_abs_water = params_fit_final_Ts.correct_absolute(params_ref_fit_Ts, water_concentration)

# nice display
fig = plt.figure(700)
fig.clf()
fig.canvas.set_window_title("Final results")

# fit plot
ax = plt.subplot(1, 2, 1)
sim.disp_fit(ax, data, params_fit_final, simtool, True, True)

# now the bargraph
ax = plt.subplot(2, 2, 2)
sim.disp_bargraph(ax, [params_fit_final_Ts_abs_water], [params_fit_final_Ts_abs_water * 0.0], ['Fit Ts-cor rel. Water'], False, True, False, False, xxx.p_cm, np.concatenate([metabolites_fit, lipids_fit]), 0.1)
plt.ylim([0, 200])

# %% save results to csv file

# prepare stuff to write in csv file
metabolites_names = np.array(simtool.metagroup_names) 
data_legend_caption = p._display_legends_list[0] 
header_line = ["Dataset"] + list(metabolites_names[np.concatenate([metabolites_fit, lipids_fit])]) + [metabolites_names[xxx.m_Water]]
data_line = [data_legend_caption] + list(params_fit_final[np.concatenate([metabolites_fit, lipids_fit]), xxx.p_cm]) + [params_ref_fit[xxx.m_Water, xxx.p_cm]]

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
    
