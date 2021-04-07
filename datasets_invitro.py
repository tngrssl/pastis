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

# %% 27/02/2019 - spinal cord phantom - nice WS STEAM, testing TWIX reco
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "SC phantom"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/braino_009/meas_MID49_svs_st_vapor_643_FID28469.dat",
                                "/home/tangir/crmbm/acq_twix/braino_009/meas_MID50_svs_st_vapor_643_FID28470.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0024_svs-st-vapor-643",
                                "/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0025_svs-st-vapor-643"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - OVS optim
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "original - OVS duration 3000µs"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0026_svs-st-vapor-643"]

p.dataset[1]["legend"] = "OVS duration 1000µs"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0027_svs-st-vapor-643"]

p.dataset[2]["legend"] = "OVS duration 5000µs"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0028_svs-st-vapor-643"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.job["phasing"]["offset"] = np.pi
p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - spoiler optim
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "original - spoiler 20mT/m for 500µs"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0029_svs-st-vapor-643"]

p.dataset[1]["legend"] = "spoiler 20mT/m for 250µs"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0030_svs-st-vapor-643"]

p.dataset[2]["legend"] = "spoiler 20mT/m for 100µs"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0031_svs-st-vapor-643"]

p.dataset[3]["legend"] = "spoiler 10mT/m for 500µs"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0032_svs-st-vapor-643"]

p.dataset[4]["legend"] = "spoiler 5mT/m for 500µs"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0033_svs-st-vapor-643"]

p.dataset[5]["legend"] = "spoiler 0mT/m for 500µs"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0034_svs-st-vapor-643"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.job["phasing"]["offset"] = np.pi
p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - various optim
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "original"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0035_svs-st-vapor-643"]

p.dataset[1]["legend"] = "asym->sym pulses"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0036_svs-st-vapor-643"]

p.dataset[2]["legend"] = "2048->4096 pts"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0037_svs-st-vapor-643"]

p.dataset[3]["legend"] = "0->200µs acq win shift"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0038_svs-st-vapor-643"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - TM optim
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "original - TM=43ms"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0039_svs-st-vapor-643"]

p.dataset[1]["legend"] = "TM=50ms"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0040_svs-st-vapor-643"]

p.dataset[2]["legend"] = "TM=60ms"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0041_svs-st-vapor-643"]

p.dataset[3]["legend"] = "TM=70ms"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0042_svs-st-vapor-643"]

p.dataset[4]["legend"] = "TM=80ms"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0043_svs-st-vapor-643"]

p.dataset[5]["legend"] = "TM=90ms"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0044_svs-st-vapor-643"]

p.dataset[6]["legend"] = "TM=100ms"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0045_svs-st-vapor-643"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - compare svs-st-vapor-643, eja-svs-steam and eja_svs_slaser
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "svs-st-vapor-643"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0039_svs-st-vapor-643"]

p.dataset[1]["legend"] = "eja-svs-steam"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0051_eja-svs-steam"]

p.dataset[2]["legend"] = "eja-svs-slaser"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0079_eja-svs-slaser"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 12/03/2019 - spinal cord phantom - TR optim
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "TR=1500ms"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0026_svs-st-vapor-643"]

p.dataset[1]["legend"] = "TR=2000ms"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0027_svs-st-vapor-643"]

p.dataset[2]["legend"] = "TR=2500ms"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0028_svs-st-vapor-643"]

p.dataset[3]["legend"] = "TR=3000ms"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0029_svs-st-vapor-643"]

p.dataset[4]["legend"] = "TR=3500ms"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0030_svs-st-vapor-643"]

p.dataset[5]["legend"] = "TR=4000ms"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0031_svs-st-vapor-643"]

p.dataset[6]["legend"] = "TR=4500ms"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0032_svs-st-vapor-643"]

p.dataset[7]["legend"] = "TR=5000ms"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0033_svs-st-vapor-643"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "initial STEAM"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0016_svs-st-vapor-643-optim-trig"]

p.dataset[1]["legend"] = "sym. pulses"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0019_svs-st-vapor-643-optim-trig"]

p.dataset[2]["legend"] = "increased settling time"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0022_svs-st-vapor-643-optim-trig"]

p.dataset[3]["legend"] = "increasing ramp time"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0025_svs-st-vapor-643-optim-trig"]

p.dataset[4]["legend"] = "short strong spoilers"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0028_svs-st-vapor-643-optim-trig"]

p.dataset[5]["legend"] = "mini spoiler"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0031_svs-st-vapor-643-optim-trig"]

p.dataset[6]["legend"] = "OVS on"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0034_svs-st-vapor-643-optim-trig"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - VOI traces
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0018_svs-st-vapor-643-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0021_svs-st-vapor-643-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0024_svs-st-vapor-643-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0027_svs-st-vapor-643-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0030_svs-st-vapor-643-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0033_svs-st-vapor-643-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0036_svs-st-vapor-643-optim-trig"]

legends_list = ["initial STEAM",
                "sym. pulses",
                "increased settling time",
                "increasing ramp time",
                "short strong spoilers",
                "mini spoiler",
                "OVS on"]

analyze_selectivity_range_list = [[800, 3550],
                                  [-10600, -7800],
                                  [-3650, 1850]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - sLASER
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "initial sLASER"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0040_eja-svs-slaser-optim-trig"]

p.dataset[1]["legend"] = "OVS on"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0044_eja-svs-slaser-optim-trig"]

p.dataset[2]["legend"] = "short strong spoilers"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0047_eja-svs-slaser-optim-trig"]

p.dataset[3]["legend"] = "reducing ramp time"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0051_eja-svs-slaser-optim-trig"]

p.dataset[4]["legend"] = "increasing settling time"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0054_eja-svs-slaser-optim-trig"]

p.dataset[5]["legend"] = "reducing R 20->10"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0057_eja-svs-slaser-optim-trig"]

p.dataset[6]["legend"] = "increasing R 20->25"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0060_eja-svs-slaser-optim-trig"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM versus sLASER
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "initial STEAM"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0016_svs-st-vapor-643-optim-trig"]

p.dataset[1]["legend"] = "initial sLASER"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0040_eja-svs-slaser-optim-trig"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002 - sLASER optim invitro in invivo B1/B0 conditions - VOI 1/2 profiles
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0043_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0050_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0053_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0056_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0059_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0062_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0063_eja-svs-slaser-optim-trig"]

legends_list = ["initial sLASER",
                "OVS on",
                "short strong spoilers",
                "reducing ramp time",
                "increasing settling time",
                "reducing R 20->10",
                "increasing R 20->25"]

analyze_selectivity_range_list = [[800, 3550],
                                  [-10600, -7800],
                                  [-3650, 1850]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 20/03/2019 - fatty_braino_002 - sLASER optim invitro in invivo B1/B0 conditions, R factor - VOI 2/2 profiles
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0063_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0064_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0065_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0066_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0067_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0068_eja-svs-slaser-optim-trig",
                    "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0069_eja-svs-slaser-optim-trig"]

legends_list = ["increasing R=5",
                "increasing R=10",
                "increasing R=15",
                "increasing R=20",
                "increasing R=25",
                "increasing R=30",
                "increasing R=40"]

analyze_selectivity_range_list = [[800, 3550],
                                  [-10600, -7800],
                                  [-3650, 1850]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %%20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - TWIX vs DCM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "ref 0th order phasing"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/fatty_braino_002/meas_MID53_svs_st_vapor_643_optim_trig_FID29712.dat",
                                "/home/tangir/crmbm/acq_twix/fatty_braino_002/meas_MID54_svs_st_vapor_643_optim_trig_FID29713.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0016_svs-st-vapor-643-optim-trig",
                                "/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0017_svs-st-vapor-643-optim-trig"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - sym vs asym STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "asym STEAM"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0026_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0027_svs-st-vapor-643-optim"]

p.dataset[1]["legend"] = "sym STEAM"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0024_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0025_svs-st-vapor-643-optim"]

p.job_list = [  p.job["phasing"],
                # p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - sym vs asym STEAM - VOI profiles
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0021_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0022_svs-st-vapor-643-optim"]

legends_list = ["sym STEAM",
                "asym sSTEAM"]

analyze_selectivity_range_list = [[-2000, 1600],
                                  [-12500, -9000],
                                  [-3800, 2600]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - settling time in STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "STEAM 500us settling"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0028_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0029_svs-st-vapor-643-optim"]

p.dataset[1]["legend"] = "STEAM 300us settling"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0030_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0031_svs-st-vapor-643-optim"]

p.dataset[2]["legend"] = "STEAM 10us settling"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0032_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0033_svs-st-vapor-643-optim"]

p.job_list = [  p.job["phasing"],
                # p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - settling time in STEAM - VOI profiles
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0034_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0035_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0036_svs-st-vapor-643-optim"]

legends_list = ["STEAM 10us settling",
                "STEAM 300us settling",
                "STEAM 500us settling"]

analyze_selectivity_range_list = [[-2000, 1600],
                                  [-12500, -9000],
                                  [-3800, 2600]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - OVS in STEAM
get_ipython().magic("clear")
plt.close("all")

p.dataset[0]["legend"] = "STEAM no OVS"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0006_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0007_svs-st-vapor-643-optim"]

p.dataset[1]["legend"] = "STEAM OVS 45deg (SAR limit)"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0008_svs-st-vapor-643-optim/",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0009_svs-st-vapor-643-optim"]

p.job_list = [  p.job["phasing"],
                # p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - OVS in STEAM - VOI profiles (1/2)
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0012_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0013_svs-st-vapor-643-optim"]

legends_list = ["no OVS",
                "OVS 45deg"]

analyze_selectivity_range_list = [[-2000, 1600],
                                  [-12500, -9000],
                                  [-3800, 2600]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - ramp time in STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "200us ramp"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0014_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0015_svs-st-vapor-643-optim"]

p.dataset[1]["legend"] = "100us ramp"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0016_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0017_svs-st-vapor-643-optim"]

p.dataset[2]["legend"] = "300us ramp"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0018_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0019_svs-st-vapor-643-optim"]

p.job_list = [  p.job["phasing"],
                # p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - OVS in STEAM - VOI profiles (2/2)
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0020_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0021_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0022_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0023_svs-st-vapor-643-optim"]

legends_list = ["200 us ramp",
                "100 us ramp",
                "300 us ramp",
                "100 us ramp FOV 500"]

analyze_selectivity_range_list = [[-2000, 1600],
                                  [-12500, -9000],
                                  [-3800, 2600]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - spoilers in STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "20mT/m - 500ms"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0024_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0025_svs-st-vapor-643-optim"]

p.dataset[1]["legend"] = "35mT/m - 100ms"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0027_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0026_svs-st-vapor-643-optim"]

p.dataset[2]["legend"] = "20mT/m - 100ms"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0028_svs-st-vapor-643-optim/",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0029_svs-st-vapor-643-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - spoilers in STEAM - VOI profiles
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0030_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0031_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0032_svs-st-vapor-643-optim"]

legends_list = ["20mT/m - 500ms",
                "35mT/m - 100ms",
                "20mT/m - 100ms"]

analyze_selectivity_range_list = [[-2000, 1600],
                                  [-12500, -9000],
                                  [-3800, 2600]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - VAPOR BW in STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "135Hz"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0033_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0034_svs-st-vapor-643-optim"]

p.dataset[1]["legend"] = "200Hz"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0036_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0037_svs-st-vapor-643-optim"]

p.dataset[2]["legend"] = "250Hz"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0038_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0039_svs-st-vapor-643-optim"]

p.dataset[3]["legend"] = "optim"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0043_svs-st-vapor-643-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0044_svs-st-vapor-643-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - testing 30/30/60 OVS in STEAM - VOI profiles
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0045_svs-st-vapor-643-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0046_svs-st-vapor-643-optim"]

legends_list = ["30/30/60 OVS",
                "no OVS"]

analyze_selectivity_range_list = [[-2000, 1600],
                                  [-12500, -9000],
                                  [-3800, 2600]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - R/pulse length in sLASER - VOI profiles
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0051_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0052_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0053_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0054_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0055_eja-svs-slaser-optim"]

legends_list = ["R=20",
                "R=10",
                "OVS 90",
                "+ 30mT/m spoiler",
                "+ 40mT/m spoiler"]

analyze_selectivity_range_list = [[-2000, 1600],
                                  [-12500, -9000],
                                  [-3800, 2600]]

reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - R/N pulses in sLASER - VOI profiles - MAP
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0056_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0057_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0058_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0059_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0060_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0061_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0062_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0063_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0064_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0065_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0066_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0067_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0068_eja-svs-slaser-optim",
                    "/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0069_eja-svs-slaser-optim"]

legends_list = ["N=1,R=5",
                "N=1,R=10",
                "N=1,R=15",
                "N=1,R=20",
                "N=1,R=25",
                "N=1,R=30",
                "N=2,R=5",
                "N=2,R=10",
                "N=2,R=15",
                "N=3,R=5",
                "N=3,R=10",
                "N=4,R=5",
                "N=4,R=10",
                "N=5,R=5"]

analyze_selectivity_range_list = [[800, 3550],
                                  [-10600, -7800],
                                  [-3650, 1450]]

sel_results = reco.reco_spatial_select_profile(dcm_folders_list,
                                 legends_list,
                                 analyze_selectivity_range_list)

# R and N parameters
n_r_list = []
minTE_list = []
for d in dcm_folders_list:
    s = reco.MRSData2(d)
    this_r = s.sequence.pulse_rfc_r
    this_n = s.sequence.pulse_rfc_n
    this_te = s.te
    n_r_list.append([this_r, this_n])
    minTE_list.append(this_te)

n_r_list = np.array(n_r_list)
r_list_uniq = np.unique(n_r_list[:, 0])
n_list_uniq = np.unique(n_r_list[:, 1])

# reshaping to [5,7,3,2]: 5 possible Ns, 7 possible Rs, 3 axes, 2 conditions IN/OUT
sel_results_2d = np.full([len(n_list_uniq), len(r_list_uniq), 3, 2], np.nan)
min_TE2d = np.full([len(n_list_uniq), len(r_list_uniq)], np.nan)
for this_n_r, this_sel_results, this_minTE in zip(n_r_list, sel_results[:], minTE_list):
    this_r = this_n_r[0]
    this_n = this_n_r[1]
    this_r_ind = np.where(r_list_uniq == this_r)
    this_n_ind = np.where(n_list_uniq == this_n)
    sel_results_2d[this_n_ind[0][0], this_r_ind[0][0], :, :] = this_sel_results
    min_TE2d[this_n_ind, this_r_ind] = this_minTE

# and plotting
fig, ax = plt.subplots(2, 3)

ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 4)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 5)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 6)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

# ratio
sel_results_2d_ratio = np.full([5, 7, 3], np.nan)
sel_results_2d_ratio = sel_results_2d[:, :, :, 0] / sel_results_2d[:, :, :, 1]

fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 2]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

# and the min TE heat map
fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(min_TE2d)
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('Mininum TE')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('TE (ms)')
plt.show()

# %% 07/05/2019 - fatty_braino_005, bad B1/B0 conditions - R/N pulses in sLASER TE fixed at 50ms - MAP
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "N=1,R=5"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0005_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0006_eja-svs-slaser-optim"]

p.dataset[1]["legend"] = "N=1,R=10"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0007_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0008_eja-svs-slaser-optim"]

p.dataset[2]["legend"] = "N=1,R=15"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0009_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0010_eja-svs-slaser-optim"]

p.dataset[3]["legend"] = "N=1,R=20"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0011_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0012_eja-svs-slaser-optim"]

p.dataset[4]["legend"] = "N=1,R=25"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0013_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0014_eja-svs-slaser-optim"]

p.dataset[5]["legend"] = "N=1,R=30"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0015_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0016_eja-svs-slaser-optim"]

p.dataset[6]["legend"] = "N=2,R=5"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0017_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0018_eja-svs-slaser-optim"]

p.dataset[7]["legend"] = "N=2,R=10"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0019_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0020_eja-svs-slaser-optim"]

p.dataset[8]["legend"] = "N=2,R=15"
p.dataset[8]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0021_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0022_eja-svs-slaser-optim"]

p.dataset[9]["legend"] = "N=3,R=5"
p.dataset[9]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0023_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0024_eja-svs-slaser-optim"]

p.dataset[10]["legend"] = "N=3,R=10"
p.dataset[10]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0025_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0026_eja-svs-slaser-optim"]

p.dataset[11]["legend"] = "N=4,R=5"
p.dataset[11]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0027_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0028_eja-svs-slaser-optim"]

p.dataset[12]["legend"] = "N=4,R=10"
p.dataset[12]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0029_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0030_eja-svs-slaser-optim"]

p.dataset[13]["legend"] = "N=5,R=5"
p.dataset[13]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0031_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0032_eja-svs-slaser-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.run()


# R and N parameters
n_r_list = []
for d in p.dataset:
    this_r = d["dcm"]["data"].sequence.pulse_rfc_r
    this_n = d["dcm"]["data"].sequence.pulse_rfc_n
    n_r_list.append([this_r, this_n])

n_r_list = np.array(n_r_list)
r_list_uniq = np.unique(n_r_list[:, 0])
n_list_uniq = np.unique(n_r_list[:, 1])

# this is shaped [14]: 14 scans
_, snr_results, _, _ = p.get_final_analyze_results()

# reshaping to [5,7]: 5 possible Ns, 7 possible Rs
snr_results_2d = np.full([len(n_list_uniq), len(r_list_uniq)], np.nan)
for this_n_r, this_snr in zip(n_r_list, snr_results):
    this_r = this_n_r[0]
    this_n = this_n_r[1]
    this_r_ind = np.where(r_list_uniq == this_r)
    this_n_ind = np.where(n_list_uniq == this_n)
    snr_results_2d[this_n_ind[0][0], this_r_ind[0][0]] = this_snr

# and plotting
fig, ax = plt.subplots()

ax = plt.subplot(2, 2, 1)
im = ax.imshow(snr_results_2d)
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('SNR on NAA peak')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('SNR')
plt.show()

# %% 16/05/2019 - fatty_braino_007, bad B1/B0 conditions - R/N pulses in sLASER TE minimal - MAP
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "N=1,R=5"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0006_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0007_eja-svs-slaser-optim"]

p.dataset[1]["legend"] = "N=1,R=10"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0008_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0009_eja-svs-slaser-optim"]

p.dataset[2]["legend"] = "N=1,R=15"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0010_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0011_eja-svs-slaser-optim"]

p.dataset[3]["legend"] = "N=1,R=20"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0012_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0013_eja-svs-slaser-optim"]

p.dataset[4]["legend"] = "N=1,R=25"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0014_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0015_eja-svs-slaser-optim"]

p.dataset[5]["legend"] = "N=1,R=30"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0016_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0017_eja-svs-slaser-optim"]

p.dataset[6]["legend"] = "N=2,R=5"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0018_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0019_eja-svs-slaser-optim"]

p.dataset[7]["legend"] = "N=2,R=10"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0020_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0021_eja-svs-slaser-optim"]

p.dataset[8]["legend"] = "N=2,R=15"
p.dataset[8]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0022_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0023_eja-svs-slaser-optim"]

p.dataset[9]["legend"] = "N=3,R=5"
p.dataset[9]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0024_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0025_eja-svs-slaser-optim"]

p.dataset[10]["legend"] = "N=3,R=10"
p.dataset[10]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0026_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0027_eja-svs-slaser-optim"]

p.dataset[11]["legend"] = "N=4,R=5"
p.dataset[11]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0028_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0029_eja-svs-slaser-optim"]

p.dataset[12]["legend"] = "N=4,R=10"
p.dataset[12]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0030_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0031_eja-svs-slaser-optim"]

p.dataset[13]["legend"] = "N=5,R=5"
p.dataset[13]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0032_eja-svs-slaser-optim",
                                 "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0033_eja-svs-slaser-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.run()

# R and N parameters
n_r_list = []
for d in p.dataset:
    this_r = d["dcm"]["data"].sequence.pulse_rfc_r
    this_n = d["dcm"]["data"].sequence.pulse_rfc_n
    n_r_list.append([this_r, this_n])

n_r_list = np.array(n_r_list)
r_list_uniq = np.unique(n_r_list[:, 0])
n_list_uniq = np.unique(n_r_list[:, 1])

# this is shaped [14]: 14 scans
_, snr_results, _, _ = p.get_final_analyze_results()

# reshaping to [5,7]: 5 possible Ns, 7 possible Rs
snr_results_2d = np.full([len(n_list_uniq), len(r_list_uniq)], np.nan)
for this_n_r, this_snr in zip(n_r_list, snr_results):
    this_r = this_n_r[0]
    this_n = this_n_r[1]
    this_r_ind = np.where(r_list_uniq == this_r)
    this_n_ind = np.where(n_list_uniq == this_n)
    snr_results_2d[this_n_ind[0][0], this_r_ind[0][0]] = this_snr

# and plotting
fig, ax = plt.subplots()

ax = plt.subplot(2, 2, 1)
im = ax.imshow(snr_results_2d)

ax.set_xticks(range(len(n_list_uniq)))
ax.set_xticklabels(n_list_uniq)
ax.set_xlabel('R')

ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')

ax.grid('on')
ax.set_title('SNR on NAA peak')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('SNR')
plt.show()

# %%16/05/2019 - fatty_braino_007, bad B1/B0 conditions - STEAM - testing TWIX reco
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "0th"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID30_eja_svs_slaser_optim_FID32219.dat",
                                "/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID31_eja_svs_slaser_optim_FID32220.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0032_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0033_eja-svs-slaser-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["order"] = 0
p.analyze_enable = True
p.run()

# %% 21/05/2019 - fatty_braino_007, bad B1/B0 conditions - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "steam"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0008_svs-st-vapor-643",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0009_svs-st-vapor-643"]

p.dataset[1]["legend"] = "N=1,R=5"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0015_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0017_eja-svs-slaser-optim"]

p.dataset[2]["legend"] = "N=1,R=20"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0018_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0019_eja-svs-slaser-optim"]

p.dataset[3]["legend"] = "N=2,R=10"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0020_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0021_eja-svs-slaser-optim"]

p.dataset[4]["legend"] = "N=5,R=5"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0022_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0023_eja-svs-slaser-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.analyze_enable = False
p.run()

# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "steam"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0026_svs-st-vapor-643",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0027_svs-st-vapor-643"]

p.dataset[1]["legend"] = "N=1,R=5"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0030_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0031_eja-svs-slaser-optim"]

p.dataset[2]["legend"] = "N=1,R=20"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0032_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0033_eja-svs-slaser-optim"]

p.dataset[3]["legend"] = "N=2,R=10"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0034_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0035_eja-svs-slaser-optim"]

p.dataset[4]["legend"] = "N=5,R=5"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0036_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0037_eja-svs-slaser-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.analyze_enable = False
p.run()

# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions, stronger spoilers - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "steam"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0026_svs-st-vapor-643",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0027_svs-st-vapor-643"]

p.dataset[1]["legend"] = "N=1,R=5"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0041_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0042_eja-svs-slaser-optim"]

p.dataset[2]["legend"] = "N=1,R=20"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0043_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0044_eja-svs-slaser-optim"]

p.dataset[3]["legend"] = "N=2,R=10"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0045_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0046_eja-svs-slaser-optim"]

p.dataset[4]["legend"] = "N=5,R=5"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0047_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0048_eja-svs-slaser-optim"]

p.dataset[5]["legend"] = "N=1,R=20, big spoil"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0049_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0050_eja-svs-slaser-optim"]

p.dataset[6]["legend"] = "N=5,R=5, big spoil"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0051_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0052_eja-svs-slaser-optim"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.analyze_enable = False
p.run()

# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions, even more stronger spoilers - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "N=1,R=5"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0053_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0054_eja-svs-slaser-optim"]

p.dataset[1]["legend"] = "N=1,R=20"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0055_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0056_eja-svs-slaser-optim"]

p.dataset[2]["legend"] = "N=2,R=10"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0057_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0058_eja-svs-slaser-optim"]

p.dataset[3]["legend"] = "N=5,R=5"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0059_eja-svs-slaser-optim",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0060_eja-svs-slaser-optim"]

p.dataset[4]["legend"] = "steam"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0061_svs-st-vapor-643",
                                "/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0062_svs-st-vapor-643"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.analyze_enable = False
p.run()

# %% 03/06/2019 - fatty_braino_head_coil_002, test non-WS steam
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "steam"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0010_steam-shortte-snr",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0011_steam-shortte-snr"]

p.dataset[1]["legend"] = "steam nonWS"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0012_steam-shortte-snr",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0012_steam-shortte-snr"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                # p.job["cropping"],
                p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.job["phasing"]["order"] = 0
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.analyze_enable = False
p.run()

# %% 05/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - VOI profiles - MAP
get_ipython().magic("clear")
plt.close("all")

dcm_folders_list = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0040_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0037_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0034_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0031_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0027_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0024_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0021_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0018_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0051_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0068_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0077_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0052_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0071_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0053_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0074_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0056_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0059_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0062_slaser-r-n",
                    "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0065_slaser-r-n"]

legends_list = ["N=1,R=5",
                "N=1,R=10",
                "N=1,R=15",
                "N=1,R=20",
                "N=1,R=25",
                "N=1,R=30",
                "N=1,R=35",
                "N=1,R=40",
                "N=2,R=5",
                "N=2,R=10",
                "N=2,R=15",
                "N=3,R=5",
                "N=3,R=10",
                "N=4,R=5",
                "N=4,R=10",
                "N=5,R=5",
                "N=6,R=5",
                "N=7,R=5",
                "N=8,R=5"]

analyze_selectivity_range_list = [[-6000, -2000],
                                  [-9000, -5000],
                                  [-5000, -2000]]

# this is shaped [20,3,2]: 20 scans, 3 axes, 2 conditions IN/OUT
sel_results = reco.reco_spatial_select_profile(dcm_folders_list,
                                                       legends_list,
                                                       analyze_selectivity_range_list)

# R and N parameters
n_r_list = []
minTE_list = []
for d in dcm_folders_list:
    s = reco.MRSData2(d)
    this_r = s.sequence.pulse_rfc_r
    this_n = s.sequence.pulse_rfc_n
    this_te = s.te
    n_r_list.append([this_r, this_n])
    minTE_list.append(this_te)

n_r_list = np.array(n_r_list)
r_list_uniq = np.unique(n_r_list[:, 0])
n_list_uniq = np.unique(n_r_list[:, 1])

# reshaping to [5,8,3,2]: 5 possible Ns, 8 possible Rs, 3 axes, 2 conditions IN/OUT
sel_results_2d = np.full([len(n_list_uniq), len(r_list_uniq), 3, 2], np.nan)
min_TE2d = np.full([len(n_list_uniq), len(r_list_uniq)], np.nan)
for this_n_r, this_sel_results, this_minTE in zip(n_r_list, sel_results[:], minTE_list):
    this_r = this_n_r[0]
    this_n = this_n_r[1]
    this_r_ind = np.where(r_list_uniq == this_r)
    this_n_ind = np.where(n_list_uniq == this_n)
    sel_results_2d[this_n_ind[0][0], this_r_ind[0][0], :, :] = this_sel_results
    min_TE2d[this_n_ind, this_r_ind] = this_minTE

# and plotting
fig, ax = plt.subplots(2, 3)

ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 4)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 5)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 6)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

# ratio
sel_results_2d_ratio = np.full([8, 8, 3], np.nan)
sel_results_2d_ratio = sel_results_2d[:, :, :, 0] / sel_results_2d[:, :, :, 1]

fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 0]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 1]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 2]))
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

# and the min TE heat map
fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(min_TE2d)
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('Mininum TE')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('TE (ms)')
plt.show()

# %% 05/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - sLASER spectra - MAP
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "N=1,R=5"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0038_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0039_slaser-r-n"]

p.dataset[1]["legend"] = "N=1,R=10"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0035_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0036_slaser-r-n"]

p.dataset[2]["legend"] = "N=1,R=15"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0032_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0033_slaser-r-n"]

p.dataset[3]["legend"] = "N=1,R=20"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0028_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0029_slaser-r-n"]

p.dataset[4]["legend"] = "N=1,R=25"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0025_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0026_slaser-r-n"]

p.dataset[5]["legend"] = "N=1,R=30"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0022_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0023_slaser-r-n"]

p.dataset[6]["legend"] = "N=1,R=35"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0019_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0020_slaser-r-n"]

p.dataset[7]["legend"] = "N=1,R=40"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0016_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0017_slaser-r-n"]

p.dataset[8]["legend"] = "N=2,R=5"
p.dataset[8]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0041_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0042_slaser-r-n"]

p.dataset[9]["legend"] = "N=2,R=10"
p.dataset[9]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0066_slaser-r-n",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0067_slaser-r-n"]

p.dataset[10]["legend"] = "N=2,R=15"
p.dataset[10]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0075_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0076_slaser-r-n"]

p.dataset[11]["legend"] = "N=3,R=5"
p.dataset[11]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0045_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0046_slaser-r-n"]

p.dataset[12]["legend"] = "N=3,R=10"
p.dataset[12]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0069_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0070_slaser-r-n"]

p.dataset[13]["legend"] = "N=4,R=5"
p.dataset[13]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0048_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0049_slaser-r-n"]

p.dataset[14]["legend"] = "N=4,R=10"
p.dataset[14]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0072_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0073_slaser-r-n"]

p.dataset[15]["legend"] = "N=5,R=5"
p.dataset[15]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0054_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0055_slaser-r-n"]

p.dataset[16]["legend"] = "N=5,R=6"
p.dataset[16]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0057_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0058_slaser-r-n"]

p.dataset[17]["legend"] = "N=5,R=7"
p.dataset[17]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0060_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0061_slaser-r-n"]

p.dataset[18]["legend"] = "N=5,R=8"
p.dataset[18]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0063_slaser-r-n",
                                 "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0064_slaser-r-n"]

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
                p.job["calibrating"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.settings["display_offset"] = 30000
p.run()

# R and N parameters
n_r_list = []
for d in p.dataset:
    this_r = d["dcm"]["data"].sequence.pulse_rfc_r
    this_n = d["dcm"]["data"].sequence.pulse_rfc_n
    this_te = d["dcm"]["data"].te
    n_r_list.append([this_r, this_n])

n_r_list = np.array(n_r_list)
r_list_uniq = np.unique(n_r_list[:, 0])
n_list_uniq = np.unique(n_r_list[:, 1])

# this is shaped [19]: 19 scans
_, snr_results, _, _ = p.get_final_analyze_results()

# reshaping to [8,8]: 8 possible Ns, 8 possible Rs
snr_results_2d = np.full([len(n_list_uniq), len(r_list_uniq)], np.nan)
for this_n_r, this_snr_results in zip(n_r_list, snr_results):
    this_r = this_n_r[0]
    this_n = this_n_r[1]
    this_r_ind = np.where(r_list_uniq == this_r)
    this_n_ind = np.where(n_list_uniq == this_n)
    snr_results_2d[this_n_ind[0][0], this_r_ind[0][0]] = this_snr_results

# and plotting
fig, ax = plt.subplots()

ax = plt.subplot(2, 2, 1)
im = ax.imshow(snr_results_2d)
ax.set_xticks(range(len(r_list_uniq)))
ax.set_xticklabels(r_list_uniq)
ax.set_xlabel('R')
ax.set_yticks(range(len(n_list_uniq)))
ax.set_yticklabels(n_list_uniq)
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('SNR on NAA peak')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('SNR')
plt.show()

# %% 05/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - STEAM WS spectra
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "steam WS"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0010_steam-shortte-snr",
                                "/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0012_steam-shortte-snr"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                # p.job["cropping"],
                p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]


p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.job["phasing"]["order"] = 0
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.9, 2.1]
p.job["displaying"]["range_ppm"] = [1, 6]

p.analyze_enable = True
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing non WS STEAM spectra
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "steam ws 90"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr"]

p.dataset[1]["legend"] = "steam ws 30deg"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0012_steam-shortte-snr"]

p.dataset[2]["legend"] = "steam ws 1deg"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr"]

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                # p.job["cropping"],
                p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]]

p.analyze_enable = False

p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [1.5, 2]
p.job["displaying"]["range_ppm"] = [1, 6]
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing DCM reco with right ref scans...
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "steam dcm not phased"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr",
                                "/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr"]

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False

p.job["phasing"]["order"] = 0
p.job["displaying"]["range_ppm"] = [1, 6]
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing TWIX reco with right ref scans...
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "phased"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr",
                                "/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID90_steam_shortTE_SNR+_FID33706.dat",
                                "/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID92_steam_shortTE_SNR+_FID33708.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.job["apodizing"]["damping_hz"] = 20
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing TWIX reco and non WS STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "VAPOR FA=30deg"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID91_steam_shortTE_SNR+_FID33707.dat"]

p.dataset[1]["legend"] = "VAPOR FA=1deg"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID92_steam_shortTE_SNR+_FID33708.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                # p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                p.job["cropping"],
                p.job["water_removal"],
                # p.job["calibrating"],
                p.job["displaying"]]

p.job["phasing"]["POI_range_ppm"] = [4, 5]
p.job["phasing"]["order"] = 1
p.job["apodizing"]["damping_hz"] = 20

p.analyze_enable = False
p.run()

# %%16/05/2019 - fatty_braino_007, bad B1/B0 conditions - STEAM - testing TWIX reco
get_ipython().magic("clear")
plt.close("all")

for chan_to_turn_off in range(0, 9):

    print("chan_to_turn_off=" + str(chan_to_turn_off))
    p = reco.pipeline()
    p.dataset[0]["legend"] = "off channel = " + str(chan_to_turn_off)
    p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID30_eja_svs_slaser_optim_FID32219.dat",
                                    "/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID31_eja_svs_slaser_optim_FID32220.dat"]

    p.job_list = [  p.job["phasing"],
                    p.job["scaling"],
                    # p.job["FID modulus"],
                    p.job["channel_combining"],
                    # p.job["concatenate"],
                    # p.job["zero_filling"],
                    # p.job["physio_analysis"],
                    # p.job["data_rejecting"],
                    p.job["realigning"],
                    p.job["averaging"],
                    p.job["noise_estimation"],
                    p.job["apodizing"],
                    # p.job["cropping"],
                    # p.job["water_removal"],
                    p.job["calibrating"],
                    p.job["displaying"]]

    p.job["phasing"]["order"] = 1
    p.job["channel_combining"]["phasing"] = False

    combine_weights = np.full([8, ], True)
    if(chan_to_turn_off < 8):
        combine_weights[chan_to_turn_off] = False
    p.job["channel_combining"]["weights"] = combine_weights

    p.run()

# %%27/11/2019 - phantom_test_inv, sLASER - inv test
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "WS"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID145_slaser_R_N=20+_1_longTE_SNR++++_FID47558.dat","/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID147_slaser_R_N=20+_1_longTE_SNR++++_FID47560.dat"]

p.dataset[1]["legend"] = "WS TI=770ms"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID149_slaser_R_N=20+_1_longTE_SNR++++_FID47562.dat","/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID147_slaser_R_N=20+_1_longTE_SNR++++_FID47560.dat"]

p.dataset[2]["legend"] = "WS double inv TI=770ms"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID151_slaser_R_N=20+_1_longTE_SNR++++_FID47564.dat","/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID147_slaser_R_N=20+_1_longTE_SNR++++_FID47560.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %%20/12/2019 - test_IR_sLASER - inv tests
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.dataset[0]["legend"] = "15"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "25"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0002.dcm"]

p.dataset[2]["legend"] = "35"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0003.dcm"]

p.dataset[3]["legend"] = "45"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0004.dcm"]

p.dataset[4]["legend"] = "55"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0005.dcm"]

p.dataset[5]["legend"] = "65"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0006.dcm"]

p.dataset[6]["legend"] = "75"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0007.dcm"]

p.dataset[7]["legend"] = "85"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0008.dcm"]

p.dataset[8]["legend"] = "95"
p.dataset[8]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0009.dcm"]

p.dataset[9]["legend"] = "105"
p.dataset[9]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0010.dcm"]

p.dataset[10]["legend"] = "115"
p.dataset[10]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0011.dcm"]

p.dataset[11]["legend"] = "125"
p.dataset[11]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0012.dcm"]

p.dataset[12]["legend"] = "135"
p.dataset[12]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0013.dcm"]

p.dataset[13]["legend"] = "145"
p.dataset[13]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0014.dcm"]

p.dataset[14]["legend"] = "155"
p.dataset[14]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0015.dcm"]

p.dataset[15]["legend"] = "165"
p.dataset[15]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0016.dcm"]

p.dataset[16]["legend"] = "175"
p.dataset[16]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0017.dcm"]

p.dataset[17]["legend"] = "185"
p.dataset[17]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0018.dcm"]

p.dataset[18]["legend"] = "195"
p.dataset[18]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0019.dcm"]

p.dataset[19]["legend"] = "205"
p.dataset[19]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0020.dcm"]

p.dataset[20]["legend"] = "215"
p.dataset[20]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0021.dcm"]

p.dataset[21]["legend"] = "225"
p.dataset[21]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0022.dcm"]

p.dataset[22]["legend"] = "235"
p.dataset[22]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0023.dcm"]

p.dataset[23]["legend"] = "245"
p.dataset[23]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0024.dcm"]

p.dataset[24]["legend"] = "255"
p.dataset[24]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0025.dcm"]

p.dataset[25]["legend"] = "265"
p.dataset[25]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0026.dcm"]

p.dataset[26]["legend"] = "275"
p.dataset[26]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0027.dcm"]

p.dataset[27]["legend"] = "285"
p.dataset[27]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0028.dcm"]

p.dataset[28]["legend"] = "295"
p.dataset[28]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0029.dcm"]

p.dataset[29]["legend"] = "305"
p.dataset[29]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0030.dcm"]

p.dataset[30]["legend"] = "315"
p.dataset[30]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0031.dcm"]

p.dataset[31]["legend"] = "325"
p.dataset[31]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0032.dcm"]

p.dataset[32]["legend"] = "335"
p.dataset[32]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0033.dcm"]

p.dataset[33]["legend"] = "345"
p.dataset[33]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0034.dcm"]

p.dataset[34]["legend"] = "355"
p.dataset[34]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0035.dcm"]

p.dataset[35]["legend"] = "365"
p.dataset[35]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0036.dcm"]

p.dataset[36]["legend"] = "375"
p.dataset[36]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0037.dcm"]

p.dataset[37]["legend"] = "385"
p.dataset[37]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0038.dcm"]

p.dataset[38]["legend"] = "395"
p.dataset[38]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0039.dcm"]

p.dataset[39]["legend"] = "405"
p.dataset[39]["dcm"]["files"] = ["/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0040.dcm"]

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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

p.analyze_enable = False
p.run()

# %%20/12/2019 - test_IR_sLASER - head coil - WS/noWS/FID modulus tests
plt.close('all')

p = reco.pipeline()
p.dataset[0]["legend"] = "WS"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID48889.dat",
                                "/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID70_slaser_R_N=20+_1_longTE_SNR++++_FID48891.dat"]

p.dataset[1]["legend"] = "noWS"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID69_slaser_R_N=20+_1_longTE_SNR++++_FID48890.dat",
                                "/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID70_slaser_R_N=20+_1_longTE_SNR++++_FID48891.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                p.job["cropping"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]]

p.settings["datasets_indexes"] = [0]
p.run()

p = reco.pipeline()
p.dataset[0]["legend"] = "WS"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID48889.dat",
                                "/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID70_slaser_R_N=20+_1_longTE_SNR++++_FID48891.dat"]

p.dataset[1]["legend"] = "noWS"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID69_slaser_R_N=20+_1_longTE_SNR++++_FID48890.dat",
                                "/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID70_slaser_R_N=20+_1_longTE_SNR++++_FID48891.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                p.job["cropping"],
                p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]]

p.settings["datasets_indexes"] = [1]
p.run()

# %% 22/01/2020 - test_IR_sLASER - head coil - WS tests
plt.close('all')

p = reco.pipeline()

p.dataset[0]["legend"] = "110"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID72_slaser_R_N=20+_1_longTE_SNR++++_FID49978.dat"]

p.dataset[1]["legend"] = "45"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID49979.dat"]

p.dataset[2]["legend"] = "30"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID74_slaser_R_N=20+_1_longTE_SNR++++_FID49980.dat"]

p.dataset[3]["legend"] = "10"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID75_slaser_R_N=20+_1_longTE_SNR++++_FID49981.dat"]

p.dataset[4]["legend"] = "0"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID76_slaser_R_N=20+_1_longTE_SNR++++_FID49982.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                p.job["cropping"],
                p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [2.8, 3.2]
p.job["analyzing-lw"]["POI_range_ppm"] = [2.8, 3.2]
p.job["calibrating"]["POI_shift_range_ppm"] = [4, 5]
p.job["calibrating"]["POI_shift_true_ppm"] = 4.7

p.run()

# %% 09/03/2020 - Comparing sLASER w/o offset exc
plt.close('all')

p = reco.pipeline()
p.dataset[0]["legend"] = "offset = 0"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID100_slaser_R_N=20+_1_longTE_SNR++++_FID54192.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID101_slaser_R_N=20+_1_longTE_SNR++++_FID54193.dat"]

p.dataset[1]["legend"] = "offset = -1.7ppm"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID54194.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID103_slaser_R_N=20+_1_longTE_SNR++++_FID54195.dat"]

p.dataset[2]["legend"] = "offset = -1.7ppm for WS, 0 for REF"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID54194.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID104_slaser_R_N=20+_1_longTE_SNR++++_FID54196.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
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


p.job["phasing"]["order"] = 1
p.job["zero_filling"]["npts"] = 4096 * 6
p.job["apodizing"]["damping_hz"] = 1

p.settings["display_offset"] = 0.0
p.analyze_enable = False
p.run()

# %% 11/03/2020 - Jojo's creatine phantom : Comparing different concentrations at 7T
plt.close('all')

p = reco.pipeline()
p.dataset[0]["legend"] = "sLASER on 10mM Cre"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID30_slaser_[10]_WSAT_FID53582.dat",
                                "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID31_slaser_[10]_NOWSAT_FID53583.dat"]

p.dataset[1]["legend"] = "STEAM on 10mM Cre"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID35_steam_[10]_WSAT_FID53587.dat",
                                "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID36_steam_[10]_NOWSAT_FID53588.dat"]

p.dataset[2]["legend"] = "STEAM on 25mM Cre"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID54_steam_[25]_WSAT_FID53606.dat",
                                "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID55_steam_[25]_NOWSAT_FID53607.dat"]

p.dataset[3]["legend"] = "sLASER on 25mM Cre"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID59_slaser_[25]_WSAT_FID53611.dat",
                                "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID60_slaser_[25]_NOWSAT_FID53612.dat"]

p.dataset[4]["legend"] = "STEAM on 50mM Cre"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID75_steam_[50]_WSAT_FID53627.dat",
                                "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID81_slaser_[50]_NOWSAT_FID53633.dat"]

p.dataset[5]["legend"] = "sLASER on 50mM Cre"
p.dataset[5]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID80_slaser_[50]_WSAT_FID53632.dat",
                                "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID82_steam_[50]_NOWSAT_FID53634.dat"]

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                p.job["apodizing"],
                p.job["cropping"],
                p.job["water_removal"],
                p.job["calibrating"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.job["phasing"]["order"] = 1
p.job["phasing"]["offset"] = 0.0  # np.pi
p.job["channel_combining"]["phasing"] = False

p.job["calibrating"]["POI_shift_range_ppm"] = [2.5, 3.5]
p.job["calibrating"]["POI_shift_true_ppm"] = 3.03

p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [2.5, 3.5]
p.job["analyzing-lw"]["POI_range_ppm"] = [2.5, 3.5]

p.settings["datasets_indexes"] = [2]
data_list = p.run()

# %% 26/05/2020 - Ethanol tests
plt.close('all')

p = reco.pipeline()

p.dataset[0]["legend"] = "Ethanol - STEAM"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID95_steam_shortTE_SNR+_FID54187.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID96_steam_shortTE_SNR+_FID54188.dat"]

p.dataset[1]["legend"] = "Ethanol - sLASER 50/1 TE=50ms"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID54194.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID103_slaser_R_N=20+_1_longTE_SNR++++_FID54195.dat"]

p.dataset[2]["legend"] = "Ethanol - sLASER 50/1 TE=90ms"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID105_slaser_R_N=20+_1_longTE_SNR++++_FID54197.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID106_slaser_R_N=20+_1_longTE_SNR++++_FID54198.dat"]

p.dataset[3]["legend"] = "Ethanol - sLASER 50/1 TE=130ms"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID107_slaser_R_N=20+_1_longTE_SNR++++_FID54199.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID108_slaser_R_N=20+_1_longTE_SNR++++_FID54200.dat"]

p.dataset[4]["legend"] = "Ethanol - sLASER 20/1 TE=130ms"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID109_slaser_R_N=20+_1_longTE_SNR++++_FID54201.dat",
                                "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID110_slaser_R_N=20+_1_longTE_SNR++++_FID54202.dat"]

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
                p.job["cropping"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["phasing (suspect)"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.job["phasing"]["order"] = 1
p.job["zero_filling"]["npts"] = 16384 + 1024
p.job["cropping"]["final_npts"] = 16384 + 1024
p.job["apodizing"]["damping_hz"] = 15
p.job["calibrating"]["POI_shift_range_ppm"] = [4.5, 5]
p.job["calibrating"]["POI_shift_true_ppm"] = 4.7

p.analyze_enable = False
p.settings["datasets_indexes"] = [3, 4]
p.run()

# %% 03/06/2020 - Ethanol night tests
plt.close('all')

p = reco.pipeline()

p.dataset[0]["legend"] = "Ethanol2 - 20/1 R=70"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID200_slaser_R_N=20+_1_longTE_SNR++++_FID56493.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID201_slaser_R_N=20+_1_longTE_SNR++++_FID56494.dat"]

p.dataset[1]["legend"] = "Ethanol2 - 20/1 R=20"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID203_slaser_R_N=20+_1_longTE_SNR++++_FID56496.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID204_slaser_R_N=20+_1_longTE_SNR++++_FID56497.dat"]

p.dataset[2]["legend"] = "Ethanol2 - 20/1 R=70 5ms"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID205_slaser_R_N=20+_1_longTE_SNR++++_FID56498.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID206_slaser_R_N=20+_1_longTE_SNR++++_FID56499.dat"]

p.dataset[3]["legend"] = "Ethanol2 - 20/1 TE=100 R=70"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID208_slaser_R_N=20+_1_longTE_SNR++++_FID56501.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID209_slaser_R_N=20+_1_longTE_SNR++++_FID56502.dat"]

p.dataset[4]["legend"] = "Ethanol2 - 20/1 TE=100 R=20"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID210_slaser_R_N=20+_1_longTE_SNR++++_FID56503.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID211_slaser_R_N=20+_1_longTE_SNR++++_FID56504.dat"]

p.dataset[5]["legend"] = "Ethanol2 - 20/1 TE=100 R=70 5ms"
p.dataset[5]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID212_slaser_R_N=20+_1_longTE_SNR++++_FID56505.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID213_slaser_R_N=20+_1_longTE_SNR++++_FID56506.dat"]

p.dataset[6]["legend"] = "Ethanol2 - PRESS TE=50ms"
p.dataset[6]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID218_svs_se_30_FID56511.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID219_svs_se_30_FID56512.dat"]

p.dataset[7]["legend"] = "Ethanol2 - PRESS TE=30ms"
p.dataset[7]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID221_svs_se_30_FID56514.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID222_svs_se_30_FID56515.dat"]

p.dataset[8]["legend"] = "Ethanol2 - ejaPRESS TE=30ms"
p.dataset[8]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID228_eja_svs_press_FID56521.dat",
                                "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID230_eja_svs_press_FID56523.dat"]

p.dataset[9]["legend"] = "Glu - 20/1 R=70"
p.dataset[9]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/glu_slaser_head/meas_MID263_slaser_R_N=20+_1_longTE_SNR++++_FID56556.dat",
                                "/home/tangir/crmbm/acq_twix/glu_slaser_head/meas_MID264_slaser_R_N=20+_1_longTE_SNR++++_FID56557.dat"]

p.dataset[10]["legend"] = "Glycerol - 20/1 R=70"
p.dataset[10]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID204_slaser_R_N=20+_1_longTE_SNR++++_FID56867.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID205_slaser_R_N=20+_1_longTE_SNR++++_FID56868.dat"]

p.dataset[11]["legend"] = "Glycerol - 20/1 R=20"
p.dataset[11]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID206_slaser_R_N=20+_1_longTE_SNR++++_FID56869.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID207_slaser_R_N=20+_1_longTE_SNR++++_FID56870.dat"]

p.dataset[12]["legend"] = "Glycerol - 20/1 R=70 5ms"
p.dataset[12]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID208_slaser_R_N=20+_1_longTE_SNR++++_FID56871.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID209_slaser_R_N=20+_1_longTE_SNR++++_FID56872.dat"]

p.dataset[13]["legend"] = "Ethanol - 20/1 R=70"
p.dataset[13]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID226_slaser_R_N=20+_1_longTE_SNR++++_FID56889.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID227_slaser_R_N=20+_1_longTE_SNR++++_FID56890.dat"]

p.dataset[14]["legend"] = "Ethanol - 20/1 R=20"
p.dataset[14]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID228_slaser_R_N=20+_1_longTE_SNR++++_FID56891.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID229_slaser_R_N=20+_1_longTE_SNR++++_FID56892.dat"]

p.dataset[15]["legend"] = "Ethanol - 20/1 R=70 5ms"
p.dataset[15]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID230_slaser_R_N=20+_1_longTE_SNR++++_FID56893.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID231_slaser_R_N=20+_1_longTE_SNR++++_FID56894.dat"]

p.dataset[16]["legend"] = "Eth3 - 10/2"
p.dataset[16]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID252_slaser_R_N=10_2_longTE_SNR+++_FID56915.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID253_slaser_R_N=10_2_longTE_SNR+++_FID56916.dat"]

p.dataset[17]["legend"] = "Eth3 - 20/1 6ms"
p.dataset[17]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID255_slaser_R_N=20+_1_longTE_SNR++++_FID56918.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID256_slaser_R_N=20+_1_longTE_SNR++++_FID56919.dat"]

p.dataset[18]["legend"] = "Eth3 - 20/1 7ms"
p.dataset[18]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID257_slaser_R_N=20+_1_longTE_SNR++++_FID56920.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID258_slaser_R_N=20+_1_longTE_SNR++++_FID56921.dat"]

p.dataset[19]["legend"] = "Eth3 - 20/1 8ms"
p.dataset[19]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID259_slaser_R_N=20+_1_longTE_SNR++++_FID56922.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID260_slaser_R_N=20+_1_longTE_SNR++++_FID56923.dat"]

p.dataset[20]["legend"] = "Eth3 - 20/1 9ms"
p.dataset[20]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID261_slaser_R_N=20+_1_longTE_SNR++++_FID56924.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID262_slaser_R_N=20+_1_longTE_SNR++++_FID56925.dat"]

p.dataset[21]["legend"] = "Eth3 - 20/1 10ms"
p.dataset[21]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID263_slaser_R_N=20+_1_longTE_SNR++++_FID56926.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID264_slaser_R_N=20+_1_longTE_SNR++++_FID56927.dat"]

p.dataset[22]["legend"] = "Gly3 - 10/2"
p.dataset[22]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID276_slaser_R_N=10_2_longTE_SNR+++_FID56939.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID277_slaser_R_N=10_2_longTE_SNR+++_FID56940.dat"]

p.dataset[23]["legend"] = "Gly3 - 20/1 6ms"
p.dataset[23]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID280_slaser_R_N=20+_1_longTE_SNR++++_FID56943.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID282_slaser_R_N=20+_1_longTE_SNR++++_FID56945.dat"]

p.dataset[24]["legend"] = "Gly3 - 20/1 7ms"
p.dataset[24]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID283_slaser_R_N=20+_1_longTE_SNR++++_FID56946.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID284_slaser_R_N=20+_1_longTE_SNR++++_FID56947.dat"]

p.dataset[25]["legend"] = "Gly3 - 20/1 8ms"
p.dataset[25]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID285_slaser_R_N=20+_1_longTE_SNR++++_FID56948.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID286_slaser_R_N=20+_1_longTE_SNR++++_FID56949.dat"]

p.dataset[26]["legend"] = "Gly3 - 20/1 9ms"
p.dataset[26]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID287_slaser_R_N=20+_1_longTE_SNR++++_FID56950.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID288_slaser_R_N=20+_1_longTE_SNR++++_FID56951.dat"]

p.dataset[27]["legend"] = "Gly3 - 20/1 10ms"
p.dataset[27]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID289_slaser_R_N=20+_1_longTE_SNR++++_FID56952.dat",
                                 "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID290_slaser_R_N=20+_1_longTE_SNR++++_FID56953.dat"]

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
                p.job["cropping"],
                # p.job["water_removal"],
                p.job["calibrating"],
                # p.job["phasing (suspect)"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.job["phasing"]["display"] = False
p.job["phasing"]["order"] = 0
p.job["channel_combining"]["phasing"] = False
p.job["zero_filling"]["npts"] = 16384 + 1024
p.job["cropping"]["final_npts"] = 16384 + 1024
p.job["apodizing"]["damping_hz"] = 5
p.job["calibrating"]["POI_shift_range_ppm"] = [4.5, 5]
p.job["calibrating"]["POI_shift_true_ppm"] = 4.7
p.job["displaying"]["range_ppm"] = [0, 5]

p.analyze_enable = False
p.settings["datasets_indexes"] = np.arange(16,28)
p.run()

# %% 03/06/2020 - Creatine tubes in water/agar
plt.close('all')

p = reco.pipeline()

p.dataset[0]["legend"] = "Cre - 25M - TE=50ms"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID107_slaser_R_N=20+_1_longTE_SNR++++_FID58514.dat",
                                "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID108_slaser_R_N=20+_1_longTE_SNR++++_FID58515.dat"]

p.dataset[1]["legend"] = "Cre - 25M - TE=90ms"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID109_slaser_R_N=20+_1_longTE_SNR++++_FID58516.dat",
                                "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID110_slaser_R_N=20+_1_longTE_SNR++++_FID58517.dat"]

p.dataset[2]["legend"] = "Cre - 25M - TE=130ms"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID111_slaser_R_N=20+_1_longTE_SNR++++_FID58518.dat",
                                "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID112_slaser_R_N=20+_1_longTE_SNR++++_FID58519.dat"]

p.dataset[3]["legend"] = "Cre - 10M - TE=50ms"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID248_slaser_R_N=20+_1_longTE_SNR++++_FID58655.dat",
                                "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID251_slaser_R_N=20+_1_longTE_SNR++++_FID58658.dat"]

p.dataset[4]["legend"] = "Cre - 10M - TE=90ms"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID252_slaser_R_N=20+_1_longTE_SNR++++_FID58659.dat",
                                "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID253_slaser_R_N=20+_1_longTE_SNR++++_FID58660.dat"]

p.dataset[5]["legend"] = "Cre - 10M - TE=130ms"
p.dataset[5]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID254_slaser_R_N=20+_1_longTE_SNR++++_FID58661.dat",
                                "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID255_slaser_R_N=20+_1_longTE_SNR++++_FID58662.dat"]

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
                p.job["cropping"],
                # p.job["water_removal"],
                p.job["calibrating"],
                p.job["phasing (suspect)"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.job["zero_filling"]["npts"] = 16384 + 1024
p.job["cropping"]["final_npts"] = 16384 + 1024
p.job["apodizing"]["damping_hz"] = 5
p.job["calibrating"]["POI_shift_range_ppm"] = [4.5, 5]
p.job["calibrating"]["POI_shift_true_ppm"] = 4.7
p.job["displaying"]["range_ppm"] = [0, 5]

p.analyze_enable = False
p.run()
