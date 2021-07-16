#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vivo data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.reco as reco
import mrs.log as log
get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.INFO)

# %% "brain" reconstruction template
template_name = "heart"
p = reco.pipeline()

# in general, use water peak for processing (phasing, realigning, etc.)
p.settings["POI_range_ppm"] = [4, 5]
# calibrate spectrum using the water peak
p.settings["POI_shift_range_ppm"] = [4, 5]
p.settings["POI_shift_true_ppm"] = 4.7
# measure SNR on water peak
p.settings["POI_SNR_range_ppm"] = [4, 5]
# measure linewidth on water
p.settings["POI_LW_range_ppm"] = [4, 5]

p.job_list = [  # p.job["displaying_anatomy"],
                p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["noise_estimation"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["apodizing"],
                p.job["realigning"],
                p.job["data_rejecting"],
                p.job["averaging"],
                p.job["calibrating"],
                # p.job["water_removal"],
                p.job["cropping"],
                p.job["displaying"]
                ]

# measure SNR on water during data rejecting
p.job["data_rejecting"]["POI_SNR_range_ppm"] = [4, 5]

p.save_template(template_name)

# %% 20210712 - noWS data missing
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("heart")

p.dataset[0]["legend"] = "test heart"
p.dataset[0]["raw"]["files"] = ["/home/tangir/desktop/BUCHINGER-WILHELMI FASTING/test MR protocol/meas_MID00617_FID475565_eja_svs_press.dat",
                                "/home/tangir/desktop/BUCHINGER-WILHELMI FASTING/test MR protocol/meas_MID00617_FID475565_eja_svs_press.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/desktop/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/eja_svs_press_31/IM-0008-0001.dcm",
                                "/home/tangir/desktop/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/eja_svs_press_31/IM-0008-0001.dcm"]

p.run()
p.check_analyze_results(True)

