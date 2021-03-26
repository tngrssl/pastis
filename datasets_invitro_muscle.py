#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vitro data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
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

# %% fat phantom ?
get_ipython().magic("clear")
plt.close("all")

# --- process the water-suppressed (WS) data ---
p = reco.pipeline()

p.dataset[0]["legend"] = ""
p.dataset[0]["raw"]["files"] = ["/crmbm/data_seq/users/JS/2020_1H_MRS/200514_1hmrs_phantom001/meas_MID61_eja_svs_slaser_NOVAPOR_VOX4_FID28498.dat"]

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]
                ]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]]


# let's denoize a little (5Hz exponential apodization)
p.job["apodizing"]["damping_hz"] = 5
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.job["calibrating"]["POI_shift_true_ppm"] = 4.7
p.job["calibrating"]["POI_shift_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.job["calibrating"]["POI_shift_true_ppm"] = 1.4
# p.job["calibrating"]["POI_shift_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1, 1.5]  # signal ppm range
p.job["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.job["analyzing-lw"]["POI_range_ppm"] = [4.5, 5]

# display ppm range
p.job["displaying"]["range_ppm"] = [0, 5]
# run the process pipeline
p.run()

# %% Cre tubes at 3T - NO VAPOR
get_ipython().magic("clear")
plt.close("all")

# --- process the water-suppressed (WS) data ---
p = reco.pipeline()
p.dataset[0]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0001.dcm"]
p.dataset[1]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0002.dcm"]
p.dataset[2]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0003.dcm"]
p.dataset[3]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0004.dcm"]
p.dataset[4]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0005.dcm"]
p.dataset[5]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0006.dcm"]
p.dataset[6]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0007.dcm"]
p.dataset[7]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0008.dcm"]
p.dataset[8]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0009.dcm"]
p.dataset[9]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0010.dcm"]

# found this TEs in the dicom headers
te_list = np.array([23960, 28960, 33960, 38960, 43960, 48960, 53960, 58960, 63960, 68960]) / 1000.0
for i, (d, te) in enumerate(zip(p.dataset, te_list)):
    p.dataset[i]["legend"] = str(te) + "ms"

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]
                ]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]]


# let's denoize a little (5Hz exponential apodization)
p.job["apodizing"]["damping_hz"] = 5
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.job["calibrating"]["POI_shift_true_ppm"] = 4.7
p.job["calibrating"]["POI_shift_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.job["calibrating"]["POI_shift_true_ppm"] = 1.4
# p.job["calibrating"]["POI_shift_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1, 1.5]  # signal ppm range
p.job["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.job["analyzing-lw"]["POI_range_ppm"] = [4.5, 5]

# display ppm range
p.job["displaying"]["range_ppm"] = [0, 5]
# run the process pipeline
p.run()

# %% Cre tubes at 3T - VAPOR
get_ipython().magic("clear")
plt.close("all")

# --- process the water-suppressed (WS) data ---
p = reco.pipeline()

p.dataset[0]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0001.dcm"]
p.dataset[1]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0002.dcm"]
p.dataset[2]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0003.dcm"]
p.dataset[3]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0004.dcm"]
p.dataset[4]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0005.dcm"]
p.dataset[5]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0006.dcm"]
p.dataset[6]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0007.dcm"]
p.dataset[7]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0008.dcm"]
p.dataset[8]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0009.dcm"]
p.dataset[9]["dcm"]["files"] = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0010.dcm"]

# found this TEs in the dicom headers
te_list = np.array([23960, 28960, 33960, 38960, 43960, 48960, 53960, 58960, 63960, 68960]) / 1000.0
for i, (d, te) in enumerate(zip(p.dataset, te_list)):
    p.dataset[i]["legend"] = str(te) + "ms"

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]
                ]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]]


# let's denoize a little (5Hz exponential apodization)
p.job["apodizing"]["damping_hz"] = 5
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.job["calibrating"]["POI_shift_true_ppm"] = 4.7
p.job["calibrating"]["POI_shift_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.job["calibrating"]["POI_shift_true_ppm"] = 1.4
# p.job["calibrating"]["POI_shift_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1, 1.5]  # signal ppm range
p.job["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.job["analyzing-lw"]["POI_range_ppm"] = [4.5, 5]

# display ppm range
p.job["displaying"]["range_ppm"] = [0, 5]
# run the process pipeline
p.run()

# %% T2 fit

# choose peak ppm range
peak_ppm_range = [2, 4]

# for all datasets
peak_intensity_list = []
for d in p.dataset:
    # estimate signal intensity
    _, s, _ = d["dcm"]["data"].analyze_snr_1d(peak_ppm_range)
    peak_intensity_list.append(s)

# logarize
peak_intensity_log_list = np.log(peak_intensity_list)
# regression with polyfit
b, a = np.polynomial.polynomial.polyfit(te_list, peak_intensity_log_list, 1)
# fitted T2 evolution
T2 = -1 / a
fitted_exp = np.exp(b) * np.exp(-te_list / T2)

# scatter plot
fig = plt.figure(400)
fig.clf()
axs = fig.subplots()
axs.scatter(te_list, peak_intensity_list)
axs.plot(te_list, fitted_exp, 'r--')
axs.set_xlabel('TE (ms)')
axs.set_ylabel('%.2f-%.2fppm peak intensity' % (p.job["analyzing-snr"]["POI_SNR_range_ppm"][0], p.job["analyzing-snr"]["POI_SNR_range_ppm"][1]))
axs.set_title("Found T2 = %.2fms" % T2)
axs.grid('on')
