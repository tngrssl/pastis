#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of DW-MRS (Yasmin's work).

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.reco as reco
import mrs.sim as sim
import mrs.fit as fit
import mrs.log as log
import mrs.aliases as xxx

import numpy as np
get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

# reco pipeline used here
p = reco.pipeline()
p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                p.job["time-shifting"],
                p.job["channel_combining"],
                p.job["noise_estimation"],
                p.job["zero_filling"],
                p.job["apodizing"],
                p.job["realigning"],
                # p.job["data_rejecting"],
                p.job["averaging"],
                p.job["calibrating"],
                p.job["cropping"],
                p.job["water_removal"],
                # p.job["phasing (suspect)"],
                p.job["displaying"]
                ]

p.settings["display"] = False
p.settings["POI_SNR_range_ppm"] = [4.5, 5.2]
p.job["time-shifting"]["time_shift_us"] = -500  # Yasmin estimated that time-shift (probably originating from the sequence)
p.job["apodizing"]["damping_hz"] = 20

p.save_template("sc_std_dwmrs")

# %% dataset 330_WC_P1_BRAIN
get_ipython().magic("clear")
plt.close("all")

# read data
s = reco.MRSData2("/home/tangir/crmbm/acq_twix/330_WC_P1_BRAIN/meas_MID55_svs_DW_slaser_b_WS_FID46951.dat")

# since we cannot automatically read the dw parameters from the file, we need to input them here:
s.sequence.na = 32
s.sequence.directions = [[0, 0, 1],  # made this up
                         [0, 1, 0],
                         [1, 0, 0],
                         [1, 1, 0],
                         [0, 1, 1],
                         [1, 0, 1]]

s.sequence.bvalues = [600]  # this too
s.sequence.n_b0 = 1

# prepare reconstruction pipeline
p = reco.get_dw_reco_pipeline(s, "sc_std_dwmrs", "330_WC_P1_BRAIN")

# %% run reco pipeline
p.run()

# %% prepare fit strategy

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
dwmrs_fit_strategy = fit.fit_stategy("singlets", metabolites_fit, linklock_arr, None)

# %% run fit

fit_bvalues_list = []
fit_directions_list = []
fit_results_list = []
fit_results_area_list = []

# browse though datasets
for this_dataset in p.dataset:
    # get the data
    this_data = this_dataset["raw"]["data"]

    # display the data
    this_data.display_spectrum_1d()

    # prepare simulation machine

    # metabolite db
    meta_bs = sim.metabolite_basis_set()
    meta_bs.initialize()

    # sequence
    seq = this_data.sequence
    seq.initialize(meta_bs)
    metabolites_fit = dwmrs_fit_strategy.metabolites

    # area integration

    # for water-suppressed data
    prefittool = fit.prefit_tool(this_data, seq)
    prefittool.area_integration_peak_ranges = [0.1, 0.1, 0.1]
    prefittool.area_integration_peaks = [xxx.m_NAA_CH3, xxx.m_Cr_CH3, xxx.m_Cho_CH3]
    prefittool.display_enable = False
    prefittool.initialize()
    _, params_area_pnorm = prefittool.run()

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
    fittool.params_max[:, xxx.p_dd] = 50
    # let the MMs go above
    fittool.params_max[xxx.m_All_MMs, xxx.p_dd] = 400
    # with early minimal damping
    fittool.params_init[:, xxx.p_dd] = fittool.params_min[:, xxx.p_dd] * 1.1

    # initial concentrations
    fittool.params_init[:, xxx.p_cm] = 0.0
    fittool.params_init[metabolites_fit, xxx.p_cm] = 0.1

    # frequency shifts for metabolites
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
    linklock_arr = dwmrs_fit_strategy.linklock
    fittool.params_init.linklock[:] = linklock_arr.copy()

    # numerical optimization parameters
    fittool.display_range_ppm = [0, 6]
    fittool.display_frequency = 2

    fittool.display_enable = False

    # run the fit
    fittool.initialize()
    [params_fit_final, optim_result] = fittool.run()

    # store
    fit_bvalues_list = fit_bvalues_list + [this_data.sequence.bvalues]
    fit_directions_list = fit_directions_list + [this_data.sequence.directions]
    fit_results_list = fit_results_list + [params_fit_final]
    fit_results_area_list = fit_results_area_list + [params_area_pnorm]

# %% calculate ADCs

metabolites_interest = np.array([xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3])
metabolites_names = np.array(params_fit_final.get_meta_names())[metabolites_interest]

# extract concentrations for metabolites of interest
fit_results_np = np.concatenate(fit_results_list, axis=1)
fit_results_np = fit_results_np[metabolites_interest, xxx.p_cm:-1:4]

fit_results_area_np = np.concatenate(fit_results_area_list, axis=1)
fit_results_area_np = fit_results_area_np[metabolites_interest, xxx.p_cm:-1:4]

# data were in this order (see reco.get_dw_reco_pipeline):
# b0
# for directions
#    for bvalues

# separate b0 from the rest
fit_directions_np = np.array(fit_directions_list[1:])
fit_bvalues_np = np.unique(fit_bvalues_list[1:])
#
fit_results_b0 = fit_results_np[:, 0]
fit_results_area_b0 = fit_results_area_np[:, 0]
#
fit_results_np = fit_results_np[:, 1:]
fit_results_area_np = fit_results_area_np[:, 1:]

# build final vector (metabolites, directions, bvalues)
fit_results_3d = np.zeros([len(metabolites_interest), len(fit_directions_np), len(fit_bvalues_np) + 1])
fit_results_area_3d = np.zeros([len(metabolites_interest), len(fit_directions_np), len(fit_bvalues_np) + 1])

# b0 (repeating it for each dir, will make ADC fit easier)
r = 0
fit_results_3d[:, :, 0] = np.tile(fit_results_b0, [len(fit_directions_np), 1]).T
fit_results_area_3d[:, :, 0] = np.tile(fit_results_area_b0, [len(fit_directions_np), 1]).T

# iterating to gather the rest
for d in range(len(fit_directions_np)):
    for b in range(len(fit_bvalues_np)):
        fit_results_3d[:, d, b + 1] = fit_results_np[:, r]
        fit_results_area_3d[:, d, b + 1] = fit_results_area_np[:, r]
        r = r + 1

adc_fit = np.zeros([len(metabolites_interest), len(fit_directions_np)])
adc_area = np.zeros([len(metabolites_interest), len(fit_directions_np)])

# for each metabolite
for m in range(len(metabolites_interest)):
    # for each direction
    for d in range(len(fit_directions_np)):
        # calculate ADC for fit
        x = [0] + list(fit_bvalues_np)
        y_fit = np.log(np.squeeze(fit_results_3d[m, d, :]))
        pf_fit = np.polyfit(x, y_fit, 1)
        adc_fit[m, d] = -pf_fit[0]
        # calculate ADC for area
        x = [0] + list(fit_bvalues_np)
        y_area = np.log(np.squeeze(fit_results_area_3d[m, d, :]))
        pf_area = np.polyfit(x, y_area, 1)
        adc_area[m, d] = -pf_area[0]

        # display
        fig = plt.figure()
        ax = fig.subplots(1, 1)
        ax.set_title("%s | direction %s | ADC_fit = %.2fum2/ms ADC_area = %.2fum2/ms" % (metabolites_names[m], str(fit_directions_np[d]), adc_fit[m, d], adc_area[m, d]))
        ax.plot(x, y_fit, 'rx-', label="fit")
        ax.plot(x, y_area, 'bx-', label="area")
        ax.legend()

    # mean ADC
    print("%s | ADC_fit = %.2fum2/ms ADC_area = %.2fum2/ms" % (metabolites_names[m], np.mean(adc_fit[m, :]), np.mean(adc_area[m, :])))







