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

import pdb

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = reco.data_db("/home/tangir/crmbm/data_reco/sc_db.pkl")
fdb = fit.data_db("/home/tangir/crmbm/data_fit/sc_db.pkl")

# %% retrieve data to process


def select_func_twix(d, p):
    r = ("IR" not in d.display_label and
         d.data_ref is not None)
    return(r)


data_list, data_pipeline_list = rdb.get_datasets(select_func_twix)

for data, data_pipeline in zip(data_list, data_pipeline_list):
    data.display_spectrum_1d()

    # %% set metabolites to fit, water concentration, csv file

    # some fitting parameters
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
        xxx.m_Water])

    water_concentration = 55000.0  # mmol/kg

    # %% prepare simulation machine

    # metabolite db
    meta_bs = sim.metabolite_basis_set()
    meta_bs.initialize()

    # sequence
    # take sequence from dataset, usually sLASER
    # seq = data.sequence
    # seq.pulse_rfc_real_shape_enable = False
    # seq.pulse_exc_duration = 0.001
    # seq.pulse_rfc_duration = 0.001
    # seq.spoiler_duration = 0.001

    # make it simple: SE sequence
    seq = sim.mrs_seq_press(data.sequence.te, data.sequence.tr, data.sequence.nuclei, data.sequence.npts, data.sequence.fs, data.sequence.f0)
    seq.initialize(meta_bs)

    # %% area integration

    # for non water-suppressed data
    prefittool = fit.prefit_tool(data.data_ref, seq)
    prefittool.area_integration_peak_ranges = [0.5]
    prefittool.area_integration_peaks = [xxx.m_Water]
    prefittool.initialize()
    params_ref_area, params_ref_area_pnorm = prefittool.run()

    # for water-suppressed data
    prefittool = fit.prefit_tool(data, seq)
    prefittool.area_integration_peak_ranges = [0.1, 0.1, 0.1]
    prefittool.area_integration_peaks = [xxx.m_NAA_CH3, xxx.m_Cr_CH3, xxx.m_Cho_CH3]
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

    # run the fit
    fittool.initialize()
    [params_ref_fit, params_ref_fit_CRBs_abs, params_ref_fit_CRBs_rel, optim_result_ref] = fittool.run()

    # %% fit water-suppressed data

    seq.initialize(meta_bs)

    fittool = fit.fit_tool(data, seq)

    # default fitting bounds from muscle template
    fittool.params_min = fittool.params_min.set_default_min()
    fittool.params_max = fittool.params_max.set_default_max()

    # ranges for concentration
    fittool.params_min[:, xxx.p_cm] = 0.09
    fittool.params_max[:, xxx.p_cm] = 300.0
    fittool.params_max[xxx.m_Water, xxx.p_cm] = 100000.0

    # linewidth bounds for metabolites
    fittool.params_min[xxx.m_All_MBs, xxx.p_dd] = 5.0
    fittool.params_max[xxx.m_All_MBs, xxx.p_dd] = 90.0

    # initial concentrations
    fittool.params_init[:, xxx.p_cm] = 0.0
    fittool.params_init[metabolites_fit, xxx.p_cm] = 0.1

    # frequency shifts for metabolites
    fittool.params_min[xxx.m_All_MBs, xxx.p_df] = -10.0
    fittool.params_max[xxx.m_All_MBs, xxx.p_df] = 10.0

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
    fit.disp_bargraph(ax, [params_fit_final_Ts_abs_water], [params_fit_final_Ts_abs_water * params_fit_CRBs_rel / 100.0], ['Fit Ts-cor rel. Water'], False, False, False, False, xxx.p_cm, metabolites_fit, 0.2)
    plt.ylim([0, 50])

    ax = plt.subplot(2, 2, 4)
    fit.disp_bargraph(ax, [params_area_pnorm_T2cor_abs], [params_area_pnorm_T2cor_abs * 0.0], ['Area Ts-cor rel. Water'], False, False, False, False, xxx.p_cm, metabolites_fit, 0.2)
    plt.ylim([0, 50])

    # adjust subplots
    fig.subplots_adjust(left=0.02, bottom=0.1, right=0.98, top=0.95, wspace=0.2, hspace=None)

    # %% save all this shit

    this_par_dict = {}
    this_par_dict["params_ref_area"] = params_ref_area
    this_par_dict["params_ref_area_pnorm"] = params_ref_area_pnorm
    this_par_dict["params_ref_area_pnorm_T2cor"] = params_ref_area_pnorm_T2cor

    this_par_dict["params_area"] = params_area
    this_par_dict["params_area_pnorm"] = params_area_pnorm
    this_par_dict["params_area_pnorm_T2cor"] = params_area_pnorm_T2cor

    this_par_dict["params_ref_fit"] = params_area
    this_par_dict["params_ref_fit_CRBs_abs"] = params_ref_fit_CRBs_abs
    this_par_dict["params_ref_fit_CRBs_rel"] = params_ref_fit_CRBs_rel
    this_par_dict["optim_result_ref"] = optim_result_ref

    this_par_dict["params_fit_final"] = params_fit_final
    this_par_dict["params_fit_CRBs_abs"] = params_fit_CRBs_abs
    this_par_dict["params_fit_CRBs_rel"] = params_fit_CRBs_rel
    this_par_dict["optim_result"] = optim_result

    this_par_dict["params_fit_final_Ts"] = params_fit_final_Ts
    this_par_dict["params_ref_fit_Ts"] = params_ref_fit_Ts
    this_par_dict["params_fit_final_Ts_abs_water"] = params_fit_final_Ts_abs_water
    this_par_dict["optim_result"] = optim_result

    fdb.save(data, data_pipeline, fittool, this_par_dict)

