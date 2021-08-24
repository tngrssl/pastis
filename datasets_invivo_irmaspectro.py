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

# %% "heart_noWS_missing" reconstruction template
template_name = "heart_noWS_missing"
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
                # p.job["phasing"],
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

p.job["channel_combining"]["phasing"] = True

# measure SNR on water during data rejecting
p.job["data_rejecting"]["POI_SNR_range_ppm"] = [4, 5]

p.save_template(template_name)

# %% "heart" reconstruction template
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
                #p.job["realigning"],
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

# %% 12/07/2021 - noWS data missing
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("heart_noWS_missing")

p.dataset[0]["legend"] = "PRESS 1"
p.dataset[0]["raw"]["files"] = ["/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/meas_MID00618_FID475566_eja_svs_press.dat"]
p.dataset[0]["dcm"]["files"] = ["/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/eja_svs_press_31/IM-0008-0001.dcm",
                                "/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/eja_svs_press_31/IM-0008-0001.dcm"]

p.dataset[1]["legend"] = "PRESS 2"
p.dataset[1]["raw"]["files"] = ["/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/meas_MID00617_FID475565_eja_svs_press.dat",
                                "/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/meas_MID00617_FID475565_eja_svs_press.dat"]
p.dataset[1]["dcm"]["files"] = ["/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/eja_svs_press_32/IM-0009-0001.dcm",
                                "/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/eja_svs_press_32/IM-0009-0001.dcm"]

p.dataset[2]["legend"] = "SVS_SE 2"

p.dataset[2]["dcm"]["files"] = ["/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/svs_se_30_BW2000_256_BH_1RR_35/IM-0012-0001.dcm",
                                "/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/BUCHINGER-WILHELMI FASTING/test MR protocol/Methordirm_Spectrocardiac/20210712/svs_se_30_BW2000_256_BH_1RR_ref_34/IM-0011-0001.dcm"]

p.settings["datasets_indexes"] = [0, 1]
p.run()
p.check_analyze_results(True)

# %% 23/08/2021 - testsurplace and TEST SPECTRO COEUR PRESS vapor
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("heart")

p.dataset[0]["legend"] = "#1 PRESS"
p.dataset[0]["raw"]["files"] = ["/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/testsurplace/meas_MID00093_FID489453_eja_svs_press_VAPOR_100REP.dat",
                                "/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/testsurplace/meas_MID00094_FID489454_eja_svs_press_VAPOROFF_1REP.dat"]

p.dataset[1]["legend"] = "#2 sLASER"
p.dataset[1]["raw"]["files"] = ["/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/testsurplace/TEST SPECTRO COEUR PRESS vapor/meas_MID00123_FID489483_eja_svs_slaser_VAPOR_100REP.dat",
                                "/crmbm/data_seq/users/SSSR/Data/210715_IRMASspectro/testsurplace/TEST SPECTRO COEUR PRESS vapor/meas_MID00124_FID489484_eja_svs_slaser_VAPOROFF_1REP.dat"]

p.run()
p.check_analyze_results(True)
