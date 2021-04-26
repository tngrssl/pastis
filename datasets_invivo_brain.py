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

# display stuff?
display_stuff = True

# stop pipeline if the reco went bad
raise_error_on_bad_reco = True

# template to use here (see definitions below)
reco_template = "brain"

# %% "brain" reconstruction template
template_name = "brain"

p = reco.pipeline()
p.settings["storage_file"] = "/home/tangir/crmbm/acq_db/%s.pkl" % template_name
p.settings["POI_range_ppm"] = [4.5, 5.2]
p.settings["POI_shift_range_ppm"] = [1.8, 2.2]
p.settings["POI_shift_true_ppm"] = 2.008
p.settings["POI_LW_range_ppm"] = [1.8, 2.2]
p.settings["POI_LW_range_ppm"] = [4.5, 5.2]
p.settings["allowed_apodization"] = 1.0
p.settings["display"] = display_stuff

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

p.job["cropping"]["final_npts"] = 2048
p.job["displaying"]["apodization_factor"] = 5.0

# SNR like LCModel...
p.job["analyzing_snr"]["half_factor"] = True
p.job["ref_data_analyzing_snr"]["half_factor"] = True

p.save_template(template_name)

# %% 20/06/2019 - 296_ym_p1_brainmoelle - Yasmin :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline(reco_template)

p.dataset[0]["legend"] = "brain - sLASER R:N=5:5"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID81_slaser_R_N=5_5+_shortTE_SNR++_FID33880.dat",
                                "/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID80_slaser_R_N=5_5+_shortTE_SNR++_FID33879.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0014_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0013_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "brain - sLASER R:N=25:1"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID79_slaser_R_N=20+_1_longTE_SNR++++_FID33878.dat",
                                "/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat"]
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0012_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0009_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[2]["legend"] = "brain - STEAM"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID71_steam_shortTE_SNR+_FID33870.dat",
                                "/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID69_steam_shortTE_SNR+_FID33868.dat"]
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0008_steam-shortte-snr/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0006_steam-shortte-snr/original-primary_e09_0001.dcm"]


p.settings["datasets_indexes"] = [0, 1]
p.run()
p.check_analyze_results()
p.save_datasets()

# %% 25/06/2019 - 296_ym_p1_brainmoelle - FID modulus tests
get_ipython().magic("clear")

p = reco.pipeline(reco_template)

p.dataset[0]["legend"] = "brain - sLASER no VAPOR + conventionnal process"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat"]

p.job_list.insert(11, p.job["water_removal"])
p.job_list.remove(p.job["data_rejecting"])
p.settings["POI_range_ppm"] = [4.5, 5.2]
p.run()

p = reco.pipeline(reco_template)

p.dataset[0]["legend"] = "brain - sLASER no VAPOR + FID process"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat"]

p.job_list.insert(2, p.job["FID modulus"])
p.job_list.insert(12, p.job["water_removal"])
p.job_list.remove(p.job["data_rejecting"])
p.settings["POI_range_ppm"] = [4.5, 5.2]
p.run()

# %% 27/08/2019 - 308-rs-p1-moelle - Ocha
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline(reco_template)

p.dataset[0]["legend"] = "brain - sLASER 20:1 resp trig"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID210_slaser_R_N=20+_1_longTE_SNR++++_FID38955.dat",
                                "/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID211_slaser_R_N=20+_1_longTE_SNR++++_FID38956.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0024_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0025_slaser-r-n/original-primary_e09_0001.dcm"]

p.run()
p.check_analyze_results()
p.save_datasets()

# %% 23/01/2019 - 347-re-p1-moelle - Renaud
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline(reco_template)

p.dataset[0]["legend"] = "sLASER 10:2 IR TE=40ms"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID216_slaser_R_N=10_2_longTE_SNR+++_FID50575.dat",
                                "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID218_slaser_R_N=10_2_longTE_SNR+++_FID50577.dat"]

p.dataset[1]["legend"] = "sLASER 10:2 IR TE=30ms"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID221_slaser_R_N=10_2_longTE_SNR+++_FID50580.dat",
                                "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID222_slaser_R_N=10_2_longTE_SNR+++_FID50581.dat"]

p.dataset[2]["legend"] = "STEAM IR TE=3ms"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID224_steam_shortTE_SNR+_FID50583.dat",
                                "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID226_steam_shortTE_SNR+_FID50585.dat"]

p.dataset[3]["legend"] = "sLASER 10:2 TE=30ms"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID228_slaser_R_N=10_2_longTE_SNR+++_FID50587.dat",
                                "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID229_slaser_R_N=10_2_longTE_SNR+++_FID50588.dat"]
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/347-re-p1-moelle/20200123/01_0022_slaser-r-n",
                                "/home/tangir/crmbm/acq/347-re-p1-moelle/20200123/01_0023_slaser-r-n"]

p.settings["datasets_indexes"] = 3
p.run()
p.check_analyze_results()
p.save_datasets()
