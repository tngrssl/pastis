#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vivo spinal cord data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.reco as reco
import mrs.db as db
import mrs.log as log
get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = db.data_db("/home/tangir/crmbm/acq_db/sc.pkl")

# display stuff?
display_stuff = False

# standard spinal cord template
p = reco.pipeline()
p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel-combining"],
                # p.job["concatenate"],
                p.job["noise-estimation"],
                p.job["zero-filling"],
                # p.job["physio-analysis"],
                p.job["apodizing"],
                p.job["realigning"],
                p.job["data-rejecting"],
                p.job["averaging"],
                p.job["calibrating"],
                # p.job["water-removal"],
                p.job["cropping"],
                p.job["displaying"]
                ]

p.job["data-rejecting"]["auto_method_list"] = [reco.data_rejection_method.AUTO_AMPLITUDE,
                                                reco.data_rejection_method.AUTO_LINEWIDTH,
                                                reco.data_rejection_method.AUTO_FREQUENCY,
                                                reco.data_rejection_method.AUTO_PHASE]

p.settings["display"] = display_stuff
p.save_template("sc_std_nows")

# standard spinal cord template with concatenate feature (for the first 2 crappy datasets)
p = reco.pipeline("sc_std_nows")
p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                # p.job["FID modulus"],
                p.job["channel-combining"],
                # p.job["concatenate"],
                p.job["noise-estimation"],
                p.job["zero-filling"],
                # p.job["physio-analysis"],
                p.job["apodizing"],
                p.job["realigning"],
                p.job["data-rejecting"],
                p.job["averaging"],
                p.job["calibrating"],
                # p.job["water-removal"],
                p.job["cropping"],
                p.job["displaying"]
                ]

p.job["data-rejecting"]["auto_method_list"] = [   reco.data_rejection_method.AUTO_AMPLITUDE,
                                                        reco.data_rejection_method.AUTO_LINEWIDTH,
                                                        reco.data_rejection_method.AUTO_FREQUENCY,
                                                        reco.data_rejection_method.AUTO_PHASE]

p.settings["display"] = display_stuff
p.save_template("sc_std_nows_concatenate")

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - concatenated STEAM #1 :(
# not included in group, the protocol was not ead, RFC pulses were too short...
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows_concatenate")

p.dataset[0]["legend"] = "STEAM #1"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID101_svs_st_vapor_643_optim_trig_FID29464.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID99_svs_st_vapor_643_optim_trig_FID29462.dat"]

p.dataset[1]["legend"] = "STEAM #1"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID104_svs_st_vapor_643_optim_trig_FID29467.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID103_svs_st_vapor_643_optim_trig_FID29466.dat"]

p.dataset[2]["legend"] = "STEAM #1"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID107_svs_st_vapor_643_optim_trig_FID29470.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID106_svs_st_vapor_643_optim_trig_FID29469.dat"]

p.dataset[3]["legend"] = "STEAM #1"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID110_svs_st_vapor_643_optim_trig_FID29473.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID109_svs_st_vapor_643_optim_trig_FID29472.dat"]

# p.run()
# p.check_analyze_results()
# p.save_datasets(rdb)

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - concatenated sLASER #1 :(
# not included in group, the protocol was not ead, RFC pulses were too short...
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows_concatenate")

p.dataset[0]["legend"] = "sLASER #1"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID113_eja_svs_slaser_optim_trig_FID29476.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID112_eja_svs_slaser_optim_trig_FID29475.dat"]

p.dataset[1]["legend"] = "sLASER #1"
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID116_eja_svs_slaser_optim_trig_FID29479.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID115_eja_svs_slaser_optim_trig_FID29478.dat"]

p.dataset[2]["legend"] = "sLASER #1"
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID119_eja_svs_slaser_optim_trig_FID29482.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID118_eja_svs_slaser_optim_trig_FID29481.dat"]

p.dataset[3]["legend"] = "sLASER #1"
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID122_eja_svs_slaser_optim_trig_FID29485.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID121_eja_svs_slaser_optim_trig_FID29484.dat"]

# p.run()
# p.check_analyze_results()
# p.save_datasets(rdb)

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - sLASER #2 :(
# not included in group, the protocol was not ead, RFC pulses were too short...
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER #2"
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID146_eja_svs_slaser_optim_trig_FID29509.dat",
                                "/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID145_eja_svs_slaser_optim_trig_FID29508.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/291-vs-moelle-spectro-p1/20190315/01_0043_eja-svs-slaser-optim-trig",
                                "/home/tangir/crmbm/acq/291-vs-moelle-spectro-p1/20190315/01_0042_eja-svs-slaser-optim-trig"]

# p.run()
# p.check_analyze_results()
# p.save_datasets(rdb)

# %% 26/06/2019 - 296_ym_p1_brainmoelle - Yasmin :)
# STEAM is shit, sLASER is ok
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "STEAM"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0012_steam-shortte-snr",
                                "/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0011_steam-shortte-snr"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID147_steam_shortTE_SNR+_FID34181.dat",
                                "/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID146_steam_shortTE_SNR+_FID34180.dat"]

p.dataset[1]["legend"] = "sLASER R:N=20:1"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0016_slaser-r-n",
                                "/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0015_slaser-r-n"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID157_slaser_R_N=20+_1_longTE_SNR++++_FID34191.dat",
                                "/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID155_slaser_R_N=20+_1_longTE_SNR++++_FID34189.dat"]

p.job["data-rejecting"]["auto_method_list"] = None
p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 16/07/2019 - 300-pm-p1-moelle - Pelayo :)
# forgot to acquire the ref scans... replace them with non-ws data
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER R:N=25:1"
p.dataset[0]["comment"] = "No REF scan !"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/300-pm-p1-moelle/20190716/01_0010_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/300-pm-p1-moelle/20190716/01_0010_slaser-r-n/original-primary_e09_0001.dcm"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/300-pm-p1-moelle/meas_MID62_slaser_R_N=20+_1_longTE_SNR++++_FID35515.dat",
                                "/home/tangir/crmbm/acq_twix/300-pm-p1-moelle/meas_MID62_slaser_R_N=20+_1_longTE_SNR++++_FID35515.dat"]

p.dataset[1]["legend"] = "sLASER R:N=25:1 trig"
p.dataset[0]["comment"] = "No REF scan !"
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/300-pm-p1-moelle/20190716/01_0011_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/300-pm-p1-moelle/20190716/01_0011_slaser-r-n/original-primary_e09_0001.dcm"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/300-pm-p1-moelle/meas_MID63_slaser_R_N=20+_1_longTE_SNR++++_FID35516.dat",
                                "/home/tangir/crmbm/acq_twix/300-pm-p1-moelle/meas_MID63_slaser_R_N=20+_1_longTE_SNR++++_FID35516.dat"]

p.job["phasing"]["using_ref_data"] = False
p.job["channel-combining"]["using_ref_data"] = False
p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 14/08/2019 - 304-ka-p1-moelle - Karen :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "crappy"
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/304-ka-p1-moelle/20190814/01_0013_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/304-ka-p1-moelle/20190814/01_0013_slaser-r-n/original-primary_e09_0001.dcm"]

p.job["phasing"]["using_ref_data"] = False
p.job["channel-combining"]["using_ref_data"] = False
p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 21/08/2019 - 307-ap-p1-moelle - Ariane :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1 cardiac trig"
p.dataset[0]["resp_bpm"] = 20
p.dataset[0]["heart_bpm"] = 55
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0011_slaser-r-n",
                                "/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0012_slaser-r-n"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID101_slaser_R_N=20+_1_longTE_SNR++++_FID38622.dat",
                                "/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID38623.dat"]

p.dataset[1]["legend"] = "sLASER 20:1 resp trig"
p.dataset[1]["resp_bpm"] = 20
p.dataset[1]["heart_bpm"] = 55
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0013_slaser-r-n",
                                "/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0014_slaser-r-n"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID104_slaser_R_N=20+_1_longTE_SNR++++_FID38625.dat",
                                "/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID105_slaser_R_N=20+_1_longTE_SNR++++_FID38626.dat"]

p.dataset[2]["legend"] = "sLASER 20:1 no trig"
p.dataset[2]["resp_bpm"] = 20
p.dataset[2]["heart_bpm"] = 55
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0015_slaser-r-n",
                                "/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0016_slaser-r-n"]
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID106_slaser_R_N=20+_1_longTE_SNR++++_FID38627.dat",
                                "/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID107_slaser_R_N=20+_1_longTE_SNR++++_FID38628.dat"]

p.dataset[3]["legend"] = "sLASER 10:2 repos. + resp trig"
p.dataset[3]["resp_bpm"] = 20
p.dataset[3]["heart_bpm"] = 55
p.dataset[3]["dcm"]["files"] = ["/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0018_slaser-r-n",
                                "/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0019_slaser-r-n"]
p.dataset[3]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID123_slaser_R_N=10_2_longTE_SNR+++_FID38644.dat",
                                "/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID124_slaser_R_N=10_2_longTE_SNR+++_FID38645.dat"]

p.job["realigning"]["moving_averages"] = 2
p.job["data-rejecting"]["moving_averages"] = 2
p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 27/08/2019 - 308-rs-p1-moelle - Ocha :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1 resp trig"
p.dataset[0]["resp_bpm"] = 15
p.dataset[0]["heart_bpm"] = 65
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0008_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0009_slaser-r-n/original-primary_e09_0001.dcm"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID165_slaser_R_N=20+_1_longTE_SNR++++_FID38910.dat",
                                "/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID166_slaser_R_N=20+_1_longTE_SNR++++_FID38911.dat"]

p.settings["POI_range_ppm"] = [4.5, 4.8]
p.job["realigning"]["moving_averages"] = 2
p.job["data-rejecting"]["auto_method_list"] = [reco.data_rejection_method.AUTO_AMPLITUDE,
                                               reco.data_rejection_method.AUTO_LINEWIDTH,
                                               reco.data_rejection_method.AUTO_FREQUENCY]

p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 29/08/2019 - 310-mg-p1-moelle - Maxime :s
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 13
p.dataset[0]["heart_bpm"] = 70
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/310-mg-p1-moelle/20190829/01_0006_slaser-r-n",
                                "/home/tangir/crmbm/acq/310-mg-p1-moelle/20190829/01_0007_slaser-r-n"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/310-mg-p1-moelle/meas_MID140_slaser_R_N=20+_1_longTE_SNR++++_FID39212.dat",
                                "/home/tangir/crmbm/acq_twix/310-mg-p1-moelle/meas_MID142_slaser_R_N=20+_1_longTE_SNR++++_FID39214.dat"]

# do not realign
if(p.job["realigning"] in p.job_list):
    p.job_list.remove(p.job["realigning"])

p.job["data-rejecting"]["moving_averages"] = 8

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 05/09/2019 - 311-sl-p1-moelle - Simon :))
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 16
p.dataset[0]["heart_bpm"] = 55
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/311-sl-p1-moelle/20190905/01_0020_slaser-r-n",
                                "/home/tangir/crmbm/acq/311-sl-p1-moelle/20190905/01_0018_slaser-r-n"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID128_slaser_R_N=20+_1_longTE_SNR++++_FID39740.dat",
                                "/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID124_slaser_R_N=20+_1_longTE_SNR++++_FID39736.dat"]

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 23/09/2019 - 313-ft-p1-moelle - Fransiska :|
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 12
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/313-ft-p1-moelle/20190923/01_0011_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/313-ft-p1-moelle/20190923/01_0010_slaser-r-n/original-primary_e09_0001.dcm"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/313-ft-p1-moelle/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID41500.dat",
                                "/home/tangir/crmbm/acq_twix/313-ft-p1-moelle/meas_MID65_slaser_R_N=20+_1_longTE_SNR++++_FID41497.dat"]

p.job["realigning"]["moving_averages"] = 2
p.job["data-rejecting"]["moving_averages"] = 2

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 25/09/2019 - 314-yt-p1-moelle - Yolanda :)))
# lw reduced, snr slight loss
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 17
p.dataset[0]["heart_bpm"] = 75
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0010_slaser-r-n",
                                "/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0009_slaser-r-n"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID83_slaser_R_N=20+_1_longTE_SNR++++_FID41681.dat",
                                "/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID81_slaser_R_N=20+_1_longTE_SNR++++_FID41679.dat"]
p.dataset[0]["physio-file"] = "/home/tangir/crmbm/acq_physio/314_YT_P1_MOELLE_2.resp"

p.dataset[1]["legend"] = "sLASER 5:5"
p.dataset[1]["resp_bpm"] = 17
p.dataset[1]["heart_bpm"] = 65
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0012_slaser-r-n",
                                "/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0011_slaser-r-n"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID88_slaser_R_N=5_5+_shortTE_SNR++_FID41686.dat",
                                "/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID86_slaser_R_N=5_5+_shortTE_SNR++_FID41684.dat"]

p.job["realigning"]["moving_averages"] = 2
p.job["data-rejecting"]["moving_averages"] = 2

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 03/10/2019 - 316-ap-p1-moelle - Anissa :)
# ok for me
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 11
p.dataset[0]["heart_bpm"] = 65
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0009_slaser-r-n",
                                "/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0008_slaser-r-n"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID42_slaser_R_N=20+_1_longTE_SNR++++_FID42201.dat",
                                "/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID40_slaser_R_N=20+_1_longTE_SNR++++_FID42199.dat"]
p.dataset[0]["physio-file"] = "/home/tangir/crmbm/acq_physio/316_AP_P1_MOELLE.resp"

p.dataset[1]["legend"] = "sLASER 5:5"
p.dataset[1]["resp_bpm"] = 11
p.dataset[1]["heart_bpm"] = 65
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0012_slaser-r-n",
                                "/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0011_slaser-r-n"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID49_slaser_R_N=5_5+_shortTE_SNR++_FID42208.dat",
                                "/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID47_slaser_R_N=5_5+_shortTE_SNR++_FID42206.dat"]
p.dataset[1]["physio-file"] = "/home/tangir/crmbm/acq_physio/316_AP_P1_MOELLE.resp"

p.settings["POI_range_ppm"] = [4.5, 4.8]
p.settings["POI_shift_range_ppm"] = [4.5, 4.8]
p.settings["POI_LW_range_ppm"] = [4.5, 4.8]

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 17/10/2019 - 319-fc-p1-moelle - Fernando :)
# passing ok, LWdcm = LWraw
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 18
p.dataset[0]["heart_bpm"] = 58
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/319-fc-p1-moelle/20191017/01_0015_slaser-r-n",
                                "/home/tangir/crmbm/acq/319-fc-p1-moelle/20191017/01_0013_slaser-r-n"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID142_slaser_R_N=20+_1_longTE_SNR++++_FID43720.dat",
                                "/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID138_slaser_R_N=20+_1_longTE_SNR++++_FID43716.dat"]

p.dataset[1]["legend"] = "sLASER 10:2"
p.dataset[0]["resp_bpm"] = 16
p.dataset[0]["heart_bpm"] = 59
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/319-fc-p1-moelle/20191017/01_0017_slaser-r-n",
                                "/home/tangir/crmbm/acq/319-fc-p1-moelle/20191017/01_0016_slaser-r-n"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID147_slaser_R_N=10_2_longTE_SNR+++_FID43725.dat",
                                "/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID145_slaser_R_N=10_2_longTE_SNR+++_FID43723.dat"]

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 05/11/2019 - 328-af-p1-moelle - Anne :)
#
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 18
p.dataset[0]["heart_bpm"] = 60
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/328-af-p1-moelle/20191105/01_0011_slaser-r-n",
                                "/home/tangir/crmbm/acq/328-af-p1-moelle/20191105/01_0012_slaser-r-n"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/328-af-p1-moelle/meas_MID64_slaser_R_N=20+_1_longTE_SNR++++_FID45771.dat",
                                "/home/tangir/crmbm/acq_twix/328-af-p1-moelle/meas_MID66_slaser_R_N=20+_1_longTE_SNR++++_FID45773.dat"]

p.job["realigning"]["moving_averages"] = 4
p.job["data-rejecting"]["moving_averages"] = 1

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 08/11/2019 - 329-pi-p1-moelle - Pujalina :)
# dcm better than raw, but ok
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 20
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID169_slaser_R_N=20+_1_longTE_SNR++++_FID46233.dat",
                                "/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID170_slaser_R_N=20+_1_longTE_SNR++++_FID46234.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/329-pi-p1-moelle/20191108/01_0008_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/329-pi-p1-moelle/20191108/01_0009_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "sLASER 10:2"
p.dataset[1]["resp_bpm"] = 23
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID180_slaser_R_N=10_2_longTE_SNR+++_FID46244.dat",
                                "/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID181_slaser_R_N=10_2_longTE_SNR+++_FID46245.dat"]
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/329-pi-p1-moelle/20191108/01_0012_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/329-pi-p1-moelle/20191108/01_0013_slaser-r-n/original-primary_e09_0001.dcm"]

p.job["realigning"]["moving_averages"] = 2
p.job["data-rejecting"]["moving_averages"] = 2

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 26/11/2019 - 333-sc-p1-moelle - Shirley :(
# TWIX data is corrupted ! :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 15
p.dataset[0]["heart_bpm"] = 62
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0008_slaser-r-n",
                                "/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0010_slaser-r-n"]
# p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/333-sc-p1-moelle/meas_MID123_slaser_R_N=20+_1_longTE_SNR++++_FID47359.dat"]

p.dataset[1]["legend"] = "sLASER 20:1 IR"
p.dataset[1]["resp_bpm"] = 15
p.dataset[1]["heart_bpm"] = 62
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0009_slaser-r-n",
                                "/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0010_slaser-r-n"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/333-sc-p1-moelle/meas_MID126_slaser_R_N=20+_1_longTE_SNR++++_FID47362.dat"]

# water peak very small and close to 5ppm artefact
p.settings["POI_range_ppm"] = [4.5, 4.8]
p.settings["POI_shift_range_ppm"] = [4.5, 4.8]
p.settings["POI_LW_range_ppm"] = [4.5, 4.8]

p.job["phasing"]["offset"] = 3.1416

p.settings["datasets_indexes"] = 0
p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 09/12/2019 - 336-nb-p1-moelle - Naouelle :s
# looks ok
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1"
p.dataset[0]["resp_bpm"] = 14
p.dataset[0]["heart_bpm"] = 50
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0010_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0011_slaser-r-n/original-primary_e09_0001.dcm"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID72_slaser_R_N=20+_1_longTE_SNR++++_FID48203.dat",
                                "/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID75_slaser_R_N=20+_1_longTE_SNR++++_FID48206.dat"]

p.dataset[1]["legend"] = "sLASER 20:1 IR"
p.dataset[1]["resp_bpm"] = 14
p.dataset[1]["heart_bpm"] = 50
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0012_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0011_slaser-r-n/original-primary_e09_0001.dcm"]
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID78_slaser_R_N=20+_1_longTE_SNR++++_FID48209.dat",
                                "/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID75_slaser_R_N=20+_1_longTE_SNR++++_FID48206.dat"]

# probably bad ref scan: only NA=1, no phase cycling
# p.job["channel-combining"]["using_ref_data"] = False

p.job["realigning"]["moving_averages"] = 4
p.job["realigning"]["inter_corr_mode"] = True
p.job["data-rejecting"]["moving_averages"] = 2

p.settings["datasets_indexes"] = 0
p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 10/12/2019 - 338-ro-p1-moelle - Rischa :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1 WS"
p.dataset[0]["resp_bpm"] = 20
p.dataset[0]["heart_bpm"] = 85
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/338-ro-p1-moelle/20191210/01_0007_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/338-ro-p1-moelle/20191210/01_0008_slaser-r-n/original-primary_e09_0001.dcm"]
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID112_slaser_R_N=20+_1_longTE_SNR++++_FID48494.dat",
                                "/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID114_slaser_R_N=20+_1_longTE_SNR++++_FID48496.dat"]

p.dataset[1]["legend"] = "sLASER 20:1 noWS"
p.dataset[1]["resp_bpm"] = 20
p.dataset[1]["heart_bpm"] = 85
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID114_slaser_R_N=20+_1_longTE_SNR++++_FID48496.dat",
                                "/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID114_slaser_R_N=20+_1_longTE_SNR++++_FID48496.dat"]

p.dataset[2]["legend"] = "sLASER 20:1 IR"
p.dataset[2]["resp_bpm"] = 20
p.dataset[2]["heart_bpm"] = 85
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID118_slaser_R_N=20+_1_longTE_SNR++++_FID48500.dat",
                                "/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID116_slaser_R_N=20+_1_longTE_SNR++++_FID48498.dat"]

p.settings["datasets_indexes"] = 0
p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 28/01/2019 - 300-pm-p2-moelle - Pelayo P2 :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20:1 (REF with OVS)"
p.dataset[0]["resp_bpm"] = 13
p.dataset[0]["heart_bpm"] = 55
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID50926.dat",
                                "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID69_slaser_R_N=20+_1_longTE_SNR++++_FID50927.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/300-pm-p2-moelle/20200128/01_0011_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/300-pm-p2-moelle/20200128/01_0012_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "STEAM IR (REF with OVS)"
p.dataset[1]["resp_bpm"] = 13
p.dataset[1]["heart_bpm"] = 55
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID62_steam_shortTE_SNR+_FID50920.dat",
                                "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID64_steam_shortTE_SNR+_FID50922.dat"]

p.dataset[2]["legend"] = "sLASER IR 5:5 (REF with OVS)"
p.dataset[2]["resp_bpm"] = 13
p.dataset[2]["heart_bpm"] = 55
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID73_slaser_R_N=5_5+_shortTE_SNR++_FID50931.dat",
                                "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID74_slaser_R_N=5_5+_shortTE_SNR++_FID50932.dat"]

p.settings["datasets_indexes"] = 0
p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 06/02/2019 - 349-ap-p1-moelle - Ahmad Fajar :|
# 3-5ppm region fucked  up
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 10:2 (REF with OVS)"
p.dataset[0]["resp_bpm"] = 18
p.dataset[0]["heart_bpm"] = 51
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID313_slaser_R_N=10_2_longTE_SNR+++_FID51947.dat",
                                "/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID314_slaser_R_N=10_2_longTE_SNR+++_FID51948.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/349-ap-p1-moelle/20200206/01_0006_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/349-ap-p1-moelle/20200206/01_0007_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "sLASER 10:2 (REF without OVS)"
p.dataset[1]["resp_bpm"] = 18
p.dataset[1]["heart_bpm"] = 51
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID313_slaser_R_N=10_2_longTE_SNR+++_FID51947.dat",
                                "/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID315_slaser_R_N=10_2_longTE_SNR+++_FID51949.dat"]
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/349-ap-p1-moelle/20200206/01_0006_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/349-ap-p1-moelle/20200206/01_0008_slaser-r-n/original-primary_e09_0001.dcm"]

p.job["realigning"]["moving_averages"] = 4
p.job["data-rejecting"]["moving_averages"] = 2

p.settings["datasets_indexes"] = 0
p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 24/02/2019 - 355-st-p1-moelle - Steven :)
# raw looks actually better than dcm, but not the estimated snr
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER TE=52ms"
p.dataset[0]["resp_bpm"] = 11
p.dataset[0]["heart_bpm"] = 85
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/355-st-p1-moelle/meas_MID164_slaser_R_N=20+_1_longTE_SNR++++_FID53261.dat",
                                "/home/tangir/crmbm/acq_twix/355-st-p1-moelle/meas_MID166_slaser_R_N=20+_1_longTE_SNR++++_FID53263.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/355-st-p1-moelle/20200224/01_0008_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/355-st-p1-moelle/20200224/01_0009_slaser-r-n/original-primary_e09_0001.dcm"]

p.job["realigning"]["moving_averages"] = 7
p.job["data-rejecting"]["moving_averages"] = 2

p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 04/03/2020 - 304-ka-p2-moelle - Karen P2 :(
# dataset #0 contains big water artefact, data quality is horrible, could not make it happen
# dataset #1, a bit better but still not good enough
# very particular processing: I do not use the ref scan for phasing but the 5ppm artefact which has the best SNR

get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "1st try (30-40Hz water LW)"
p.dataset[0]["resp_bpm"] = 15
p.dataset[0]["heart_bpm"] = 75
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID177_slaser_R_N=20+_1_longTE_SNR++++_FID53952.dat",
                                "/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID178_slaser_R_N=20+_1_longTE_SNR++++_FID53953.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/304-ka-p2-moelle/20200304/01_0016_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/304-ka-p2-moelle/20200304/01_0017_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "2nd try (25Hz water LW)"
p.dataset[1]["resp_bpm"] = 15
p.dataset[1]["heart_bpm"] = 75
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID185_slaser_R_N=20+_1_longTE_SNR++++_FID53960.dat",
                                "/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID187_slaser_R_N=20+_1_longTE_SNR++++_FID53962.dat"]
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/304-ka-p2-moelle/20200304/01_0019_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/304-ka-p2-moelle/20200304/01_0020_slaser-r-n/original-primary_e09_0001.dcm"]


# water peak very small and close to 5ppm artefact
p.job["phasing"]["using_ref_data"] = False

p.settings["POI_range_ppm"] = [5, 5.5]
p.settings["POI_shift_range_ppm"] = [4.5, 4.8]
p.settings["POI_LW_range_ppm"] = [4.5, 4.8]

p.job["realigning"]["moving_averages"] = 2
p.job["realigning"]["inter_corr_mode"] = True

# data rejection is confused with artefact!
p.job["data-rejecting"]["moving_averages"] = 4
p.job["data-rejecting"]["auto_method_list"] = [reco.data_rejection_method.AUTO_AMPLITUDE,
                                               reco.data_rejection_method.AUTO_LINEWIDTH]

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 30/05/2020 - 311-sl-p2-moelle - Simon P2 :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20/1 NA=128 trig"
p.dataset[0]["resp_bpm"] = 22
p.dataset[0]["heart_bpm"] = 55
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID152_slaser_R_N=20+_1_longTE_SNR++++_FID56036.dat",
                                "/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID153_slaser_R_N=20+_1_longTE_SNR++++_FID56037.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/311-sl-p2-moelle/20200529/01_0009_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/311-sl-p2-moelle/20200529/01_0010_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[1]["legend"] = "sLASER 20/1 NA=128 notrig"
p.dataset[1]["resp_bpm"] = 21
p.dataset[1]["heart_bpm"] = 55
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID156_slaser_R_N=20+_1_longTE_SNR++++_FID56040.dat",
                                "/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID157_slaser_R_N=20+_1_longTE_SNR++++_FID56041.dat"]
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/311-sl-p2-moelle/20200529/01_0011_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/311-sl-p2-moelle/20200529/01_0012_slaser-r-n/original-primary_e09_0001.dcm"]

p.job["realigning"]["moving_averages"] = 2
p.job["data-rejecting"]["moving_averages"] = 2
p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 09/06/2020 - 336-nb-p2-moelle - Naouelle P2 :)
# dataset #2 is bad, processed it differently for the 3 others
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 10/2 NA=64 notrig"
p.dataset[0]["resp_bpm"] = 18
p.dataset[0]["heart_bpm"] = 60
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID122_slaser_R_N=10_2_longTE_SNR+++_FID57329.dat",
                                "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID124_slaser_R_N=10_2_longTE_SNR+++_FID57331.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0013_slaser-r-n",
                                "/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0014_slaser-r-n"]

p.job["realigning"]["moving_averages"] = 2
p.job["realigning"]["inter_corr_mode"] = True
p.job["data-rejecting"]["moving_averages"] = 2

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20/1 NA=128 trig"
p.dataset[0]["resp_bpm"] = 17
p.dataset[0]["heart_bpm"] = 58
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID98_slaser_R_N=20+_1_longTE_SNR++++_FID57305.dat",
                                "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID100_slaser_R_N=20+_1_longTE_SNR++++_FID57307.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0008_slaser-r-n",
                                "/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0009_slaser-r-n"]

p.dataset[1]["legend"] = "sLASER 20/1 NA=64 notrig"
p.dataset[1]["resp_bpm"] = 17
p.dataset[1]["heart_bpm"] = 58
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID103_slaser_R_N=20+_1_longTE_SNR++++_FID57310.dat",
                                "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID105_slaser_R_N=20+_1_longTE_SNR++++_FID57312.dat"]
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0010_slaser-r-n",
                                "/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0011_slaser-r-n"]

p.dataset[2]["legend"] = "sLASER 5/5 NA=64 notrig"
p.dataset[2]["resp_bpm"] = 18
p.dataset[2]["heart_bpm"] = 60
p.dataset[2]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID127_slaser_R_N=5_5+_shortTE_SNR++_FID57334.dat",
                                "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID129_slaser_R_N=5_5+_shortTE_SNR++_FID57336.dat"]
p.dataset[2]["dcm"]["files"] = ["/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0015_slaser-r-n",
                                "/home/tangir/crmbm/acq/336-nb-p2-moelle/20200609/01_0016_slaser-r-n"]

p.job["realigning"]["moving_averages"] = 2
p.job["realigning"]["inter_corr_mode"] = False
p.job["data-rejecting"]["moving_averages"] = 2

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 11/06/2020 - 319-fc-p2-moelle - Fernando P2 :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 10/2 NA=64 trig"
p.dataset[0]["resp_bpm"] = 16
p.dataset[0]["heart_bpm"] = 65
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/319-fc-p2-moelle/meas_MID72_slaser_R_N=10_2_longTE_SNR+++_FID57445.dat",
                                "/home/tangir/crmbm/acq_twix/319-fc-p2-moelle/meas_MID74_slaser_R_N=10_2_longTE_SNR+++_FID57447.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/319-fc-p2-moelle/20200611/01_0012_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/319-fc-p2-moelle/20200611/01_0013_slaser-r-n/original-primary_e09_0001.dcm"]

p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 15/06/2020 - 313-ft-p2-moelle - Fransiska P2 :|
# did not pass quality check, raw data looks contaminated with some shit
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20/1"
p.dataset[0]["resp_bpm"] = 17
p.dataset[0]["heart_bpm"] = 73
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/313-ft-p2-moelle/meas_MID239_slaser_R_N=20+_1_longTE_SNR++++_FID57752.dat",
                                "/home/tangir/crmbm/acq_twix/313-ft-p2-moelle/meas_MID240_slaser_R_N=20+_1_longTE_SNR++++_FID57753.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/313-ft-p2-moelle/20200615/01_0011_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/313-ft-p2-moelle/20200615/01_0012_slaser-r-n/original-primary_e09_0001.dcm"]

p.job["realigning"]["moving_averages"] = 2
p.job["data-rejecting"]["moving_averages"] = 4
p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 19/06/2020 - 333-sc-p2-moelle - Shirley P2 :|
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20/1"
p.dataset[0]["resp_bpm"] = 16
p.dataset[0]["heart_bpm"] = 65
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/333-sc-p2-moelle/meas_MID180_slaser_R_N=20+_1_longTE_SNR++++_FID58587.dat",
                                "/home/tangir/crmbm/acq_twix/333-sc-p2-moelle/meas_MID181_slaser_R_N=20+_1_longTE_SNR++++_FID58588.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/333-sc-p2-moelle/20200619/01_0008_slaser-r-n/original-primary_e09_0001.dcm",
                                "/home/tangir/crmbm/acq/333-sc-p2-moelle/20200619/01_0009_slaser-r-n/original-primary_e09_0001.dcm"]

# water peak very small and close to 5ppm artefact
p.job["analyzing-lw"]["POI_range_ppm"] = [4.5, 4.8]
p.job["phasing"]["POI_range_ppm"] = [4.5, 4.8]
p.job["data-rejecting"]["POI_range_ppm"] = [4.5, 4.8]

p.run()
p.check_analyze_results(True)
p.save_datasets(rdb)

# %% 25/06/2020 - 314-yt-p2-moelle - Yolanda P2 :)))
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 20/1 NA=128"
p.dataset[0]["resp_bpm"] = 17
p.dataset[0]["heart_bpm"] = 80
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/314-yt-p2-moelle/meas_MID98_slaser_R_N=20+_1_longTE_SNR++++_FID59064.dat",
                                "/home/tangir/crmbm/acq_twix/314-yt-p2-moelle/meas_MID99_slaser_R_N=20+_1_longTE_SNR++++_FID59065.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/314-yt-p2-moelle/20200625/01_0008_slaser-r-n",
"/home/tangir/crmbm/acq/314-yt-p2-moelle/20200625/01_0009_slaser-r-n"]

p.run()
p.check_analyze_results()
p.save_datasets(rdb)

# %% 25/09/2020 - 349-ap-p2-moelle - Admah Fajar P2 :)))
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline("sc_std_nows")

p.dataset[0]["legend"] = "sLASER 10/2 NA=128"
p.dataset[0]["resp_bpm"] = 17
p.dataset[0]["heart_bpm"] = 70
p.dataset[0]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/349-ap-p2-moelle/meas_MID97_slaser_R_N=10_2_longTE_SNR+++_FID63227.dat",
                                "/home/tangir/crmbm/acq_twix/349-ap-p2-moelle/meas_MID98_slaser_R_N=10_2_longTE_SNR+++_FID63228.dat"]
p.dataset[0]["dcm"]["files"] = ["/home/tangir/crmbm/acq/349-ap-p2-moelle/20200925/01_0007_slaser-r-n",
                                "/home/tangir/crmbm/acq/349-ap-p2-moelle/20200925/01_0008_slaser-r-n"]

# this dataset was interrupted, cannot read it... ;(
p.dataset[1]["legend"] = "sLASER 10/2 NA=96 (repro)"
p.dataset[0]["resp_bpm"] = 17
p.dataset[0]["heart_bpm"] = 70
p.dataset[1]["raw"]["files"] = ["/home/tangir/crmbm/acq_twix/349-ap-p2-moelle/meas_MID102_slaser_R_N=10_2_longTE_SNR+++_FID63232.dat",
                                "/home/tangir/crmbm/acq_twix/349-ap-p2-moelle/meas_MID101_slaser_R_N=10_2_longTE_SNR+++_FID63231.dat"]
p.dataset[1]["dcm"]["files"] = ["/home/tangir/crmbm/acq/349-ap-p2-moelle/20200925/01_0009_slaser-r-n",
                                "/home/tangir/crmbm/acq/349-ap-p2-moelle/20200925/01_0010_slaser-r-n"]

p.settings["datasets_indexes"] = 0
p.job["realigning"]["moving_averages"] = 4
p.job["data-rejecting"]["moving_averages"] = 4

p.run()
p.check_analyze_results()
p.save_datasets(rdb)
