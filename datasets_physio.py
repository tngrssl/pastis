#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of physio data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.reco as reco
import mrs.log as log

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

# %% 27/08/2019 - 308-rs-p1-moelle - Ocha - short TR resp test
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "sLASER 20/1 SC short TR"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID177_steam_shortTE_SNR+_FID38922.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.analyze_enable = False
p.run()

# %% 05/09/2019 - 311-sl-p1-moelle - Simon, test respiration
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "STEAM noWS no trig"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID90_steam_shortTE_SNR+_FID39702.dat"]

p.dataset[1]["legend"] = "STEAM noWS resp trig 20%"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID91_steam_shortTE_SNR+_FID39703.dat"]

p.dataset[2]["legend"] = "STEAM noWS no trig"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID115_steam_shortTE_SNR+_FID39727.dat"]

p.dataset[3]["legend"] = "STEAM noWS resp trig 20%"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID116_steam_shortTE_SNR+_FID39728.dat"]

p.dataset[4]["legend"] = "STEAM noWS resp trig 10%"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID117_steam_shortTE_SNR+_FID39729.dat"]

p.dataset[5]["legend"] = "STEAM noWS resp trig 30%"
p.dataset[5]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID118_steam_shortTE_SNR+_FID39730.dat"]

p.job_list = [  p.job["phasing"],
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
                p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["zero_filling"]["npts"] = 4096 * 8
p.job["apodizing"]["damping_hz"] = 5

p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [4.5, 4.8]
p.job["analyzing-lw"]["POI_range_ppm"] = [4.5, 4.8]

p.run()

# %% 23/09/2019 - 313-ft-p1-moelle - Fransiska, test respiration
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "STEAM TR~1070ms NA=90"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/313-ft-p1-moelle/meas_MID60_steam_shortTE_SNR+_FID41492.dat"]
p.dataset[0]["physio-file"] = "/home/tangir/crmbm/acq_physio/313_FT_P1_MOELLE_1.resp"

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.job["zero_filling"]["npts"] = 4096 * 2
p.job["physio_analysis"]["POI_range_ppm"] = [4, 5]
p.job["physio_analysis"]["delta_time_ms"] = 40000.0
p.analyze_enable = False
p.run()

# %% 25/09/2019 - 314-yt-p1-moelle - Yolanda, test physio
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "sLASER 20/1 SC"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID78_steam_shortTE_SNR+_FID41676.dat"]
p.dataset[0]["physio-file"] = "/home/tangir/crmbm/acq_physio/314_YT_P1_MOELLE_1.resp"

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.job["zero_filling"]["npts"] = 4096 * 2
p.job["physio_analysis"]["POI_range_ppm"] = [4, 5]
p.job["physio_analysis"]["delta_time_ms"] = 40000.0
p.analyze_enable = False
p.run()

# %% 03/10/2019 - 316-ap-p1-moelle - Anissa, test resp
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "STEAM TR=1010ms NA=256"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID44_steam_shortTE_SNR+_FID42203.dat"]
p.dataset[0]["physio-file"] = "/home/tangir/crmbm/acq_physio/316_AP_P1_MOELLE.resp"

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.job["zero_filling"]["npts"] = 4096 * 2
p.job["physio_analysis"]["POI_range_ppm"] = [4, 5]
p.job["physio_analysis"]["delta_time_ms"] = 10000.0
p.analyze_enable = False
p.run()

# %% 05/11/2019 - 328-af-p1-moelle - Anne, test apnea/breath hold
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "STEAM TR=1010ms NA=300"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/328-af-p1-moelle/meas_MID69_steam_shortTE_SNR+_FID45776.dat"]
p.dataset[0]["physio-file"] = "/home/tangir/crmbm/acq_physio/328_AF_P1_MOELLE.resp"

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["apodizing"],
                # p.job["cropping"],
                # p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.job["zero_filling"]["npts"] = 4096 * 2
p.job["physio_analysis"]["POI_range_ppm"] = [4, 5]
p.job["physio_analysis"]["delta_time_ms"] = 40000.0
p.analyze_enable = False
p.run()
