#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vitro data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
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

# %% fat phantom ?
get_ipython().magic("clear")
plt.close("all")

# --- process the water-suppressed (WS) data ---
p = reco.pipeline()
p.data_filepaths = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID61_eja_svs_slaser_NOVAPOR_VOX4_FID28498.dat"]

"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0001.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0002.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0003.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0004.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0005.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0006.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0007.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0008.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0009.dcm",
"C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/01_0018_eja-svs-press-novapor-t2map/original-primary_e09_0010.dcm"

# found this TEs in the dicom headers
#te_list = np.array([23960, 28960, 33960, 38960, 43960, 48960, 53960, 58960, 63960, 68960]) / 1000.0

#p.display_legends = [str(te) + "ms" for te in te_list]

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                # p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]
                ]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]]


# let's denoize a little (5Hz exponential apodization)
p.jobs["apodizing"]["damping_hz"] = 5
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.jobs["calibrating"]["POI_true_ppm"] = 4.7
p.jobs["calibrating"]["POI_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.jobs["calibrating"]["POI_true_ppm"] = 1.4
# p.jobs["calibrating"]["POI_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.jobs["analyzing-snr"]["s_range_ppm"] = [1, 1.5]  # signal ppm range
p.jobs["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.jobs["analyzing-lw"]["range_ppm"] = [4.5, 5]

# display ppm range
p.jobs["displaying"]["range_ppm"] = [0, 5]
# run the process pipeline
datasets = p.run()

# %% Cre tubes at 3T - NO VAPOR
get_ipython().magic("clear")
plt.close("all")

# --- process the water-suppressed (WS) data ---
p = reco.pipeline()
p.data_filepaths = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0001.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0002.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0003.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0004.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0005.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0006.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0007.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0008.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0009.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0038_eja-svs-press-novapor-tube2-t2/original-primary_e09_0010.dcm"]

# found this TEs in the dicom headers
te_list = np.array([23960, 28960, 33960, 38960, 43960, 48960, 53960, 58960, 63960, 68960]) / 1000.0

#p.display_legends = [str(te) + "ms" for te in te_list]

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                # p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]
                ]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]]


# let's denoize a little (5Hz exponential apodization)
p.jobs["apodizing"]["damping_hz"] = 5
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.jobs["calibrating"]["POI_true_ppm"] = 4.7
p.jobs["calibrating"]["POI_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.jobs["calibrating"]["POI_true_ppm"] = 1.4
# p.jobs["calibrating"]["POI_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.jobs["analyzing-snr"]["s_range_ppm"] = [1, 1.5]  # signal ppm range
p.jobs["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.jobs["analyzing-lw"]["range_ppm"] = [4.5, 5]

# display ppm range
p.jobs["displaying"]["range_ppm"] = [0, 5]
# run the process pipeline
datasets = p.run()

# %% Cre tubes at 3T - VAPOR
get_ipython().magic("clear")
plt.close("all")

# --- process the water-suppressed (WS) data ---
p = reco.pipeline()
p.data_filepaths = ["/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0001.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0002.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0003.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0004.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0005.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0006.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0007.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0008.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0009.dcm",
"/crmbm/data_cemerem/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0039_eja-svs-press-vapor-tube2-t2-cr/original-primary_e09_0010.dcm"]

# found this TEs in the dicom headers
te_list = np.array([23960, 28960, 33960, 38960, 43960, 48960, 53960, 58960, 63960, 68960]) / 1000.0

#p.display_legends = [str(te) + "ms" for te in te_list]

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                # p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]
                ]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]]


# let's denoize a little (5Hz exponential apodization)
p.jobs["apodizing"]["damping_hz"] = 5
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.jobs["calibrating"]["POI_true_ppm"] = 4.7
p.jobs["calibrating"]["POI_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.jobs["calibrating"]["POI_true_ppm"] = 1.4
# p.jobs["calibrating"]["POI_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.jobs["analyzing-snr"]["s_range_ppm"] = [1, 1.5]  # signal ppm range
p.jobs["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.jobs["analyzing-lw"]["range_ppm"] = [4.5, 5]

# display ppm range
p.jobs["displaying"]["range_ppm"] = [0, 5]
# run the process pipeline
datasets = p.run()

# %% T2 fit

# choose peak ppm range
peak_ppm_range = [2, 4]

# for all datasets
peak_intensity_list = []
for d in datasets:
    # estimate signal intensity
    _, s, _ = d.analyze_snr_1d(peak_ppm_range)
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
axs.set_ylabel('%.2f-%.2fppm peak intensity' % (p.jobs["analyzing-snr"]["s_range_ppm"][0], p.jobs["analyzing-snr"]["s_range_ppm"][1]))
axs.set_title("Found T2 = %.2fms" % T2)
axs.grid('on')
