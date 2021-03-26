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
get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

# %% 08/08/2019 - phosphates_mix_001 - levure/lessive/javel/serum
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "un peu de lessive"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_001/phosphates-mix/20190808/01_0006_fid/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "beaucoup"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_001/phosphates-mix/20190808/01_0007_fid/original-primary_e09_0001.dcm"]

p.dataset[2]["legend"] = "beaucoup beaucoup"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_001/phosphates-mix/20190808/01_0008_fid/original-primary_e09_0001.dcm"]

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

p.job["phasing"]["POI_range_ppm"] = [-5, 5]
p.job["apodizing"]["damping_hz"] = 60.0
p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]

p.analyze_enable = False
p.run()

# %% 08/08/2019 - phosphates_mix_002 - levure/lessive/javel/serum
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "non-localized"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0004_fid/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "FA test"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0001.dcm"]

p.dataset[2]["legend"] = "FA test"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0002.dcm"]

p.dataset[3]["legend"] = "FA test"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0003.dcm"]

p.dataset[4]["legend"] = "FA test"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0004.dcm"]

p.dataset[5]["legend"] = "FA test"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0005.dcm"]

p.dataset[6]["legend"] = "FA test"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0006.dcm"]

p.dataset[7]["legend"] = "FA test"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0007.dcm"]

p.dataset[8]["legend"] = "FA test"
p.dataset[8]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0008.dcm"]

p.dataset[9]["legend"] = "FA test"
p.dataset[9]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0009.dcm"]

p.dataset[10]["legend"] = "FA test"
p.dataset[10]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0010.dcm"]

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

p.job["phasing"]["POI_range_ppm"] = [-5, 5]
p.job["apodizing"]["damping_hz"] = 60
p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]

p.analyze_enable = False
p.run()

# %% 08/08/2019 - phosphates_mix_002 - levure/lessive/javel/serum, tests STEAM, sLASER, ISIS
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "sLASER"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0010_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "ISIS"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0012_bow-isis-15/original-primary_e09_0001.dcm"]

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

p.job["phasing"]["POI_range_ppm"] = [-5, 5]
p.job["apodizing"]["damping_hz"] = 110
p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]

p.analyze_enable = False
p.run()

# %% 08/08/2019 - phosphates_mix_002 - test raw data ISIS
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "dcm ISIS"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0012_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "raw ISIS"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phosphates_mix_002/meas_MID110_bow_isis_15_FID37704.dat"]

p.job_list = [  # p.job["phasing"],
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

p.analyze_job_list = [  p.job["channel_combining"],
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["POI_range_ppm"] = [-5, 5]

# remove 1h channel from data
p.job["channel_combining"]["phasing"] = False
p.job["channel_combining"]["weights"] = [False, True]

p.job["apodizing"]["damping_hz"] = 110

p.job["displaying"]["magnitude_mode"] = False
p.job["displaying"]["range_ppm"] = [-25, 25]

p.analyze_enable = False
p.run()

# %% 13/08/2019 - phosphates_mix_003, optim exc pulse BW
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "1 kHz"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0016_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "2 kHz"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0017_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[2]["legend"] = "3 kHz"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0018_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[3]["legend"] = "4 kHz"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0019_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[4]["legend"] = "5 kHz"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0020_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[5]["legend"] = "6 kHz"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0021_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[6]["legend"] = "7 kHz"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0022_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[7]["legend"] = "8 kHz"
p.dataset[7]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0023_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[8]["legend"] = "9 kHz"
p.dataset[8]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0024_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[9]["legend"] = "10 kHz"
p.dataset[9]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0025_bow-isis-15/original-primary_e09_0001.dcm"]

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

p.analyze_job_list = [  p.job["channel_combining"],
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["POI_range_ppm"] = [-5, 5]
p.job["phasing"]["POI_range_ppm"] = [-2, +2]

# remove 1h channel from data
p.job["channel_combining"]["phasing"] = False
p.job["channel_combining"]["weights"] = [False, True]

p.analyze_enable = True
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [13, 15]
p.job["analyzing-snr"]["n_range_ppm"] = [-20, -10]
p.job["analyzing-snr"]["magnitude_mode"] = True
p.job["analyzing-lw"]["POI_range_ppm"] = [-1, 1]
p.job["analyzing-lw"]["magnitude_mode"] = True

p.job["apodizing"]["damping_hz"] = 110

p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]
p.run()

# %% 13/08/2019 - phosphates_mix_003, optim exc pulse B1
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "50V"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0026_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "100V"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0027_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[2]["legend"] = "150V"
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0028_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[3]["legend"] = "200V"
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0029_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[4]["legend"] = "250V"
p.dataset[4]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0030_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[5]["legend"] = "300V"
p.dataset[5]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0031_bow-isis-15/original-primary_e09_0001.dcm"]

p.dataset[6]["legend"] = "320V"
p.dataset[6]["dcm"]["files"] = ["/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0032_bow-isis-15/original-primary_e09_0001.dcm"]

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

p.analyze_job_list = [  p.job["channel_combining"],
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["POI_range_ppm"] = [-5, 5]
p.job["phasing"]["POI_range_ppm"] = [-2, +2]

# remove 1h channel from data
p.job["channel_combining"]["phasing"] = False
p.job["channel_combining"]["weights"] = [False, True]

p.analyze_enable = True
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [13, 15]
p.job["analyzing-snr"]["n_range_ppm"] = [-20, -10]
p.job["analyzing-snr"]["magnitude_mode"] = True
p.job["analyzing-lw"]["POI_range_ppm"] = [-1, 1]
p.job["analyzing-lw"]["magnitude_mode"] = True

p.job["apodizing"]["damping_hz"] = 110

p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]
p.run()

# %% 21/01/2020 - phantom_31p, optim B1 on levure tube
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "300"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID179_svs_st_vapor_643_optim_FID50085.dat"]

p.dataset[1]["legend"] = "275"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID180_svs_st_vapor_643_optim_FID50086.dat"]

p.dataset[2]["legend"] = "250"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID181_svs_st_vapor_643_optim_FID50087.dat"]

p.dataset[3]["legend"] = "225"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID182_svs_st_vapor_643_optim_FID50088.dat"]

p.dataset[4]["legend"] = "200"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID183_svs_st_vapor_643_optim_FID50089.dat"]

p.dataset[5]["legend"] = "175"
p.dataset[5]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID184_svs_st_vapor_643_optim_FID50090.dat"]

p.dataset[6]["legend"] = "150"
p.dataset[6]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID185_svs_st_vapor_643_optim_FID50091.dat"]

p.dataset[7]["legend"] = "125"
p.dataset[7]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID186_svs_st_vapor_643_optim_FID50092.dat"]

p.dataset[8]["legend"] = "100"
p.dataset[8]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID187_svs_st_vapor_643_optim_FID50093.dat"]

p.dataset[9]["legend"] = "75"
p.dataset[9]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID188_svs_st_vapor_643_optim_FID50094.dat"]

p.dataset[10]["legend"] = "50"
p.dataset[10]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID189_svs_st_vapor_643_optim_FID50095.dat"]

p.dataset[11]["legend"] = "25"
p.dataset[11]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID191_svs_st_vapor_643_optim_FID50097.dat"]

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

p.analyze_job_list = [  p.job["channel_combining"],
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["POI_range_ppm"] = [-1, 1]

# remove 1h channel from data
p.job["channel_combining"]["phasing"] = False
p.job["channel_combining"]["weights"] = [False, True]

p.analyze_enable = True
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [-1, 1]
p.job["analyzing-snr"]["n_range_ppm"] = [-20, -10]
p.job["analyzing-snr"]["magnitude_mode"] = False
p.job["analyzing-lw"]["POI_range_ppm"] = [-1, 1]
p.job["analyzing-lw"]["magnitude_mode"] = True

p.job["apodizing"]["damping_hz"] = 500

p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]
p.run()

# %% 21/01/2020 - phantom_31p, optim B1 on lessive
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "250"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID209_svs_st_vapor_643_optim_FID50115.dat"]

p.dataset[1]["legend"] = "200"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID210_svs_st_vapor_643_optim_FID50116.dat"]

p.dataset[2]["legend"] = "150"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID211_svs_st_vapor_643_optim_FID50117.dat"]

p.dataset[3]["legend"] = "100"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID212_svs_st_vapor_643_optim_FID50118.dat"]

p.dataset[4]["legend"] = "50"
p.dataset[4]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID213_svs_st_vapor_643_optim_FID50119.dat"]

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

p.analyze_job_list = [  p.job["channel_combining"],
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["POI_range_ppm"] = [-1, 1]

# remove 1h channel from data
p.job["channel_combining"]["phasing"] = False
p.job["channel_combining"]["weights"] = [False, True]

p.analyze_enable = True
p.job["analyzing-snr"]["POI_SNR_range_ppm"] = [-1, 1]
p.job["analyzing-snr"]["n_range_ppm"] = [-20, -10]
p.job["analyzing-snr"]["magnitude_mode"] = False
p.job["analyzing-lw"]["POI_range_ppm"] = [-1, 1]
p.job["analyzing-lw"]["magnitude_mode"] = True

p.job["apodizing"]["damping_hz"] = 100

p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]
p.run()

# %% 22/06/2020 - phantom multicomp, ISIS
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.settings["ppm0"] = 0.0

p.dataset[0]["legend"] = "ISIS on Pi+MDAP tube"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/csi-31p-phantom/meas_MID295_bow_isis_15_FID58328.dat"]

p.job_list = [  # p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel_combining"],
                # p.job["concatenate"],
                # p.job["zero_filling"],
                # p.job["physio_analysis"],
                # p.job["data_rejecting"],
                #p.job["realigning"],
                p.job["averaging"],
                p.job["noise_estimation"],
                # p.job["cropping"],
                # p.job["water_removal"],
                p.job["apodizing"],
                p.job["calibrating"],
                p.job["displaying"]]

p.analyze_job_list = [  p.job["channel_combining"],
                        # p.job["zero_filling"],
                        # p.job["realigning"],
                        p.job["averaging"],
                        # p.job["calibrating"]
                        ]

p.job["phasing"]["POI_range_ppm"] = [-1, 1]

# remove 1h channel from data
p.job["channel_combining"]["phasing"] = False
p.job["channel_combining"]["weights"] = [False, True]

p.job["realigning"]["POI_range_ppm"] = [-2, 1.2]
p.job["realigning"]["moving_averages"] = 2
p.job["realigning"]["inter_corr_mode"] = False

p.job["apodizing"]["damping_hz"] = 300

p.job["calibrating"]["POI_shift_range_ppm"] = [-5, 5]
p.job["calibrating"]["POI_shift_true_ppm"] = 0.0

p.job["displaying"]["magnitude_mode"] = True
p.job["displaying"]["range_ppm"] = [-25, 25]

p.analyze_enable = False
p.run()

