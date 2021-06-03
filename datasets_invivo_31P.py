#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of 31P in vitro data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.reco as reco
import mrs.log as log

import numpy as np
import scipy

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.INFO)

# %% 07/05/2021 - nh-t14 - Hany, fid tests, fid variable-TR for T1 estimation
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0
p.settings["display_range_ppm"] = [-25, 25]

p.dataset[0]["legend"] = "0.2ms hard pulse NOE(on)"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0010_fid"]
p.dataset[1]["legend"] = "1ms AHP pulse NOE(on)"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0011_fid"]
p.dataset[2]["legend"] = "9ms AHP pulse (calibrated) NOE(on)"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0012_fid"]
p.dataset[3]["legend"] = "7ms BIR4 pulse (calibrated) NOE(on)"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0013_fid"]
p.dataset[4]["legend"] = "0.2ms hard pulse NOE(off)"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0014_fid"]

p.dataset[5]["legend"] = "TR=15-1s"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0015_fid"]
p.dataset[6]["legend"] = "TR=15-1s"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0016_fid"]
p.dataset[7]["legend"] = "TR=15-1s"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0017_fid"]
p.dataset[8]["legend"] = "TR=15-1s"
p.dataset[8]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0018_fid"]
p.dataset[9]["legend"] = "TR=15-1s"
p.dataset[9]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0019_fid"]
p.dataset[10]["legend"] = "TR=15-1s"
p.dataset[10]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0020_fid"]
p.dataset[11]["legend"] = "TR=15-1s"
p.dataset[11]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0021_fid"]
p.dataset[12]["legend"] = "TR=15-1s"
p.dataset[12]["dcm"]["files"] = ["/home/tangir/crmbm/acq/nh-t14/20210507/01_0022_fid"]

p.job_list = [  p.job["time_shifting"],
                # p.job["phasing"],
                # p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                # p.job["noise_estimation"],
                p.job["apodizing"],
                p.job["cropping"],
                p.job["apodizing"],
                p.job["phasing_suspect"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]
                ]

p.job["time_shifting"]["time_shift_us"] = -1800
p.job["phasing"]["POI_range_ppm"] = [-2, 2]
p.job["phasing_suspect"]["range_ppm"] = [-5, 20]
# p.job["phasing_suspect"]["suspect_method"] = reco.suspect_phasing_method.ACME
p.job["apodizing"]["damping_hz"] = 60.0
p.job["cropping"]["final_npts"] = 512

p.job["calibrating"]["POI_shift_range_ppm"] = [-1, +1]
p.job["calibrating"]["POI_shift_true_ppm"] = 0

#p.job["displaying"]['magnitude_mode'] = True

p.job["analyzing_snr"]["display"] = False
p.job["analyzing_lw"]["display"] = False
p.analyze_enable = False

p.settings["datasets_indexes"] = [0, 1, 2, 3, 4]

p.run()

# %% estimate T1s from previous datasets

# collect all spectra
tr_arr = []
spectra_arr = []
for i in [5, 6, 7, 8, 9, 10, 11, 12]:
    this_data = p.dataset[i]["dcm"]["data"]
    sf_abs = np.abs(this_data.spectrum())
    spectra_arr.append(sf_abs)
    tr_arr.append(this_data.tr)

tr_arr = np.array(tr_arr)
spectra_arr = np.array(spectra_arr)

# for each chemical shift, run T1 fit
fig = plt.figure(1000)
ax = fig.subplots()
t1_ms_arr = []
for i in range(spectra_arr.shape[1]):
    signal_intensity_arr = spectra_arr[:, i]
    fit_res = scipy.optimize.curve_fit(lambda tr, m0, t1: m0 * (1 - np.exp(-tr / t1)), tr_arr, signal_intensity_arr, p0=[1.6e6, 500])
    (m0, t1_ms) = fit_res[0]
    t1_ms_arr.append(t1_ms)

    # debug
    # fit_mod = m0 * (1 - np.exp(-tr_arr / t1_ms))
    # ax.cla()
    # ax.plot(tr_arr, signal_intensity_arr, 'x')
    # ax.plot(tr_arr, fit_mod, '--')
    # plt.pause(0.2)

t1_ms_arr = np.array(t1_ms_arr)

fig = plt.figure(2000)
fig.clf()
ax = fig.subplots(2, 1, sharex='col')
ax[0].plot(this_data.frequency_axis_ppm(), spectra_arr.T)
ax[0].grid()
ax[0].set_xlim([10, -20])
ax[0].set_xlabel("chemical shift (ppm)")
ax[0].set_ylabel("magnitude (u.a)")

ax[1].plot(this_data.frequency_axis_ppm(), t1_ms_arr)
ax[1].set_ylim([0, 10000])
ax[1].grid()
ax[1].set_xlabel("chemical shift (ppm)")
ax[1].set_ylabel("T1 (ms)")

# %% 20/05/2021 - yl-t15 - Yolanda, fid variable-FA for T1 estimation
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0
p.settings["display_range_ppm"] = [-25, 25]

p.dataset[0]["legend"] = "FA=10"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/yl-t15/20210520/01_0013_fid-10"]
p.dataset[1]["legend"] = "FA=30"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/yl-t15/20210520/01_0012_fid-30"]
p.dataset[2]["legend"] = "FA=50"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/yl-t15/20210520/01_0011_fid-50"]
p.dataset[3]["legend"] = "FA=70"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/yl-t15/20210520/01_0010_fid-70"]
p.dataset[4]["legend"] = "FA=80"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/yl-t15/20210520/01_0014_fid-80"]
p.dataset[5]["legend"] = "FA=90"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/yl-t15/20210520/01_0009_fid"]

p.job_list = [  p.job["time_shifting"],
                # p.job["phasing"],
                # p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                # p.job["noise_estimation"],
                p.job["apodizing"],
                p.job["cropping"],
                p.job["apodizing"],
                p.job["phasing_suspect"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]
                ]

p.job["time_shifting"]["time_shift_us"] = -2800
p.job["phasing"]["POI_range_ppm"] = [-2, 2]
p.job["phasing_suspect"]["range_ppm"] = [-5, 20]
# p.job["phasing_suspect"]["suspect_method"] = reco.suspect_phasing_method.ACME
p.job["apodizing"]["damping_hz"] = 60.0
p.job["cropping"]["final_npts"] = 512

p.job["calibrating"]["POI_shift_range_ppm"] = [-1, +1]
p.job["calibrating"]["POI_shift_true_ppm"] = 0

#p.job["displaying"]['magnitude_mode'] = True

p.job["analyzing_snr"]["display"] = False
p.job["analyzing_lw"]["display"] = False
p.analyze_enable = False

p.run()

# %% estimate T1s from previous datasets

# collect all spectra
alpha_arr = [10, 30, 50, 70, 80, 90]
spectra_arr = []
for i in range(6):
    this_data = p.dataset[i]["dcm"]["data"]
    sf_abs = np.abs(this_data.spectrum())
    spectra_arr.append(sf_abs)

alpha_arr = np.array(alpha_arr) * np.pi / 180
spectra_arr = np.array(spectra_arr)
global_tr = this_data.tr

# for each chemical shift, run T1 fit
fig = plt.figure(1000)
ax = fig.subplots()
t1_ms_arr = []
for i in range(spectra_arr.shape[1]):
    signal_intensity_arr = spectra_arr[:, i]
    fit_res = scipy.optimize.curve_fit(lambda alpha, t1: signal_intensity_arr[-1] * np.sin(alpha) * (1 - np.exp(-global_tr / t1)) / (1 - np.cos(alpha) * np.exp(-global_tr / t1)), alpha_arr, signal_intensity_arr)
    (m0, t1_ms) = fit_res[0]
    t1_ms_arr.append(t1_ms)

    # debug
    #fit_mod = m0 * np.sin(alpha_arr) * (1 - np.exp(-global_tr / t1_ms)) / (1 - np.cos(alpha_arr) * np.exp(-global_tr / t1_ms))
    #ax.cla()
    #ax.plot(alpha_arr, signal_intensity_arr, 'x')
    #ax.plot(alpha_arr, fit_mod, '--')
    #plt.pause(0.2)

t1_ms_arr = np.array(t1_ms_arr)

fig = plt.figure(2000)
fig.clf()
ax = fig.subplots(2, 1, sharex='col')
ax[0].plot(this_data.frequency_axis_ppm(), spectra_arr.T)
ax[0].grid()
ax[0].set_xlim([10, -20])
ax[0].set_xlabel("chemical shift (ppm)")
ax[0].set_ylabel("magnitude (u.a)")

ax[1].plot(this_data.frequency_axis_ppm(), t1_ms_arr)
ax[1].set_ylim([0, 10000])
ax[1].grid()
ax[1].set_xlabel("chemical shift (ppm)")
ax[1].set_ylabel("T1 (ms)")
