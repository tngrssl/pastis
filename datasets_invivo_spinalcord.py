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
import numpy as np
get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = reco.data_db()

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - concatenated STEAM #1 :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID99_svs_st_vapor_643_optim_trig_FID29462.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID103_svs_st_vapor_643_optim_trig_FID29466.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID106_svs_st_vapor_643_optim_trig_FID29469.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID109_svs_st_vapor_643_optim_trig_FID29472.dat
"""

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID101_svs_st_vapor_643_optim_trig_FID29464.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID104_svs_st_vapor_643_optim_trig_FID29467.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID107_svs_st_vapor_643_optim_trig_FID29470.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID110_svs_st_vapor_643_optim_trig_FID29473.dat
"""

p.display_legends = """
STEAM #1
STEAM #1
STEAM #1
STEAM #1
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 50
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 150
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - concatenated sLASER #1 :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID112_eja_svs_slaser_optim_trig_FID29475.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID115_eja_svs_slaser_optim_trig_FID29478.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID118_eja_svs_slaser_optim_trig_FID29481.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID121_eja_svs_slaser_optim_trig_FID29484.dat
"""

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID113_eja_svs_slaser_optim_trig_FID29476.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID116_eja_svs_slaser_optim_trig_FID29479.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID119_eja_svs_slaser_optim_trig_FID29482.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID122_eja_svs_slaser_optim_trig_FID29485.dat
"""

p.display_legends = """
sLASER #1
sLASER #1
sLASER #1
sLASER #1
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 50
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 150
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - sLASER #2 :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID145_eja_svs_slaser_optim_trig_FID29508.dat
"""

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID146_eja_svs_slaser_optim_trig_FID29509.dat
"""

p.display_legends = """
sLASER #2
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 50
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 150
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 26/06/2019 - 296_ym_p1_brainmoelle - Yasmin
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0012_steam-shortte-snr
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0016_slaser-r-n
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID147_steam_shortTE_SNR+_FID34181.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID157_slaser_R_N=20+_1_longTE_SNR++++_FID34191.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0011_steam-shortte-snr
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190626/02_0015_slaser-r-n
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID146_steam_shortTE_SNR+_FID34180.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID155_slaser_R_N=20+_1_longTE_SNR++++_FID34189.dat
"""

p.display_legends = """
comparing DCM steam
comparing DCM sLASER R:N=20:1
comparing TWIX steam
comparing TWIX sLASER R:N=20:1
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 50
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 150
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 16/07/2019 - 300-pm-p1-moelle - Pelayo :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/300-pm-p1-moelle/20190716/01_0010_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/300-pm-p1-moelle/20190716/01_0011_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/300-pm-p1-moelle/meas_MID62_slaser_R_N=20+_1_longTE_SNR++++_FID35515.dat
/home/tangir/crmbm/acq_twix/300-pm-p1-moelle/meas_MID63_slaser_R_N=20+_1_longTE_SNR++++_FID35516.dat
"""

p.display_legends = """
sLASER R:N=25:1 (DICOM)
sLASER R:N=25:1 trig (DICOM)
sLASER R:N=25:1 (TWIX)
sLASER R:N=25:1 trig (TWIX)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 14/08/2019 - 304-ka-p1-moelle - Karen :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/304-ka-p1-moelle/20190814/01_0013_slaser-r-n/original-primary_e09_0001.dcm
"""

p.display_legends = """
crappy
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["calibrating"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 21/08/2019 - 307-ap-p1-moelle - Ariane :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0011_slaser-r-n
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0013_slaser-r-n
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0015_slaser-r-n
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0018_slaser-r-n
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID101_slaser_R_N=20+_1_longTE_SNR++++_FID38622.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID104_slaser_R_N=20+_1_longTE_SNR++++_FID38625.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID106_slaser_R_N=20+_1_longTE_SNR++++_FID38627.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID123_slaser_R_N=10_2_longTE_SNR+++_FID38644.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0012_slaser-r-n
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0014_slaser-r-n
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0016_slaser-r-n
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0019_slaser-r-n
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID38623.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID105_slaser_R_N=20+_1_longTE_SNR++++_FID38626.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID107_slaser_R_N=20+_1_longTE_SNR++++_FID38628.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID124_slaser_R_N=10_2_longTE_SNR+++_FID38645.dat
"""

p.display_legends = """
sLASER 20:1 cardiac trig
sLASER 20:1 resp trig
sLASER 20:1 no trig
sLASER 10:1 repos. + resp trig
sLASER 20:1 cardiac trig (TWIX)
sLASER 20:1 resp trig (TWIX)
sLASER 20:1 no trig (TWIX)
sLASER 10:1 repos. + resp trig (TWIX)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 25
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 150
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 27/08/2019 - 308-rs-p1-moelle - Ocha
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0008_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID165_slaser_R_N=20+_1_longTE_SNR++++_FID38910.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0009_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID166_slaser_R_N=20+_1_longTE_SNR++++_FID38911.dat
"""

p.display_legends = """
sLASER 20:1 resp trig
sLASER 20:1 resp trig (TWIX)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 75
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 150
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 29/08/2019 - 310-mg-p1-moelle - Maxime :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/310-mg-p1-moelle/meas_MID140_slaser_R_N=20+_1_longTE_SNR++++_FID39212.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/310-mg-p1-moelle/meas_MID142_slaser_R_N=20+_1_longTE_SNR++++_FID39214.dat
"""

p.display_legends = """
sLASER 20:1 SC
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 05/09/2019 - 311-sl-p1-moelle - Simon :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/311-sl-p1-moelle/20190905/01_0020_slaser-r-n
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID128_slaser_R_N=20+_1_longTE_SNR++++_FID39740.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/311-sl-p1-moelle/20190905/01_0018_slaser-r-n
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID124_slaser_R_N=20+_1_longTE_SNR++++_FID39736.dat
"""

p.display_legends = """
sLASER 20:1 (DICOM)
sLASER 20:1 (TWIX)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 23/09/2019 - 313-ft-p1-moelle - Fransiska :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/313-ft-p1-moelle/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID41500.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/313-ft-p1-moelle/meas_MID65_slaser_R_N=20+_1_longTE_SNR++++_FID41497.dat
"""

p.display_legends = """
sLASER 20:1
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["analyzing-snr"]["n_range_ppm"] = [-3, -2]
p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 50
p.jobs["data-rejecting"]["ranges"]["chemical shift (ppm)"] = 0.3
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 100
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 25/09/2019 - 314-yt-p1-moelle - Yolanda :)))
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID83_slaser_R_N=20+_1_longTE_SNR++++_FID41681.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0010_slaser-r-n
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID88_slaser_R_N=5_5+_shortTE_SNR++_FID41686.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0012_slaser-r-n
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID81_slaser_R_N=20+_1_longTE_SNR++++_FID41679.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0009_slaser-r-n
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID86_slaser_R_N=5_5+_shortTE_SNR++_FID41684.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0011_slaser-r-n
"""

p.data_physio_filepaths = """
/home/tangir/crmbm/acq_physio/314_YT_P1_MOELLE_2.resp



"""

p.display_legends = """
sLASER 20:1 (TWIX)
sLASER 20:1 (DICOM)
sLASER 5:5 (TWIX)
sLASER 5:5 (DICOM)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 03/10/2019 - 316-ap-p1-moelle - Anissa :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID42_slaser_R_N=20+_1_longTE_SNR++++_FID42201.dat
/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID49_slaser_R_N=5_5+_shortTE_SNR++_FID42208.dat
/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0009_slaser-r-n
/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0012_slaser-r-n
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID40_slaser_R_N=20+_1_longTE_SNR++++_FID42199.dat
/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID47_slaser_R_N=5_5+_shortTE_SNR++_FID42206.dat
/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0008_slaser-r-n
/home/tangir/crmbm/acq/316-ap-p1-moelle/20191003/01_0011_slaser-r-n
"""

p.data_physio_filepaths = """
/home/tangir/crmbm/acq_physio/316_AP_P1_MOELLE.resp
/home/tangir/crmbm/acq_physio/316_AP_P1_MOELLE.resp


"""

p.display_legends = """
sLASER 20:1 (TWIX)
sLASER 5:5 (TWIX)
sLASER 20:1 (DICOM)
sLASER 5:5 (DICOM)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 17/10/2019 - 319-fc-p1-moelle - Fernando :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID142_slaser_R_N=20+_1_longTE_SNR++++_FID43720.dat",
                    "/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID147_slaser_R_N=10_2_longTE_SNR+++_FID43725.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID138_slaser_R_N=20+_1_longTE_SNR++++_FID43716.dat",
                        "/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID145_slaser_R_N=10_2_longTE_SNR+++_FID43723.dat"]

p.display_legends = """
sLASER 20:1 (TWIX)
sLASER 10:2 (TWIX)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 05/11/2019 - 328-af-p1-moelle - Anne
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/328-af-p1-moelle/meas_MID64_slaser_R_N=20+_1_longTE_SNR++++_FID45771.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/328-af-p1-moelle/meas_MID66_slaser_R_N=20+_1_longTE_SNR++++_FID45773.dat
"""

p.display_legends = """
sLASER 20:1 (TWIX)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]


p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 08/11/2019 - 329-pi-p1-moelle - Pujalina
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID169_slaser_R_N=20+_1_longTE_SNR++++_FID46233.dat
/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID180_slaser_R_N=10_2_longTE_SNR+++_FID46244.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID170_slaser_R_N=20+_1_longTE_SNR++++_FID46234.dat
/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID181_slaser_R_N=10_2_longTE_SNR+++_FID46245.dat
"""

p.display_legends = """
sLASER 20:1 (TWIX)
sLASER 10:2 (TWIX)
"""

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]]

p.jobs["calibrating"]["POI_range_ppm"] = [1.8, 2.4]
p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 26/11/2019 - 333-sc-p1-moelle - Shirley :(
# TWIX data is corrupted ! :(
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = [#"/home/tangir/crmbm/acq_twix/333-sc-p1-moelle/meas_MID123_slaser_R_N=20+_1_longTE_SNR++++_FID47359.dat",
                    #"/home/tangir/crmbm/acq_twix/333-sc-p1-moelle/meas_MID126_slaser_R_N=20+_1_longTE_SNR++++_FID47362.dat",
                    "/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0008_slaser-r-n",
                    "/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0009_slaser-r-n"]

p.data_ref_filepaths = [#"",
                        #"",
                        "/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0010_slaser-r-n",
                        "/home/tangir/crmbm/acq/333-sc-p1-moelle/20191126/01_0010_slaser-r-n"]

p.display_legends = [#"sLASER 20:1 (TWIX)",
                     #"sLASER 20:1 IR (TWIX)",
                     "sLASER 20:1 (DCM)",
                     "sLASER 20:1 IR (DCM)"]

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                # p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.analyze_enable = True
p.run()
p.save(rdb)

# %% 09/12/2019 - 336-nb-p1-moelle - Naouelle :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID72_slaser_R_N=20+_1_longTE_SNR++++_FID48203.dat
/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID78_slaser_R_N=20+_1_longTE_SNR++++_FID48209.dat
/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0010_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0012_slaser-r-n/original-primary_e09_0001.dcm
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID75_slaser_R_N=20+_1_longTE_SNR++++_FID48206.dat
/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID75_slaser_R_N=20+_1_longTE_SNR++++_FID48206.dat
/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0011_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/336-nb-p1-moelle/20191209/01_0011_slaser-r-n/original-primary_e09_0001.dcm
"""

p.display_legends = """
sLASER 20:1 (TWIX)
sLASER 20:1 IR (TWIX)
sLASER 20:1 (DCM)
sLASER 20:1 IR (DCM)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["scaling"]["scaling_factor"] = 1.0
p.jobs["cropping"]["final_npts"] = 4096
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                #p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                #p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]]

p.jobs["data-rejecting"]["auto"] = False
p.jobs["scaling"]["scaling_factor"] = 0.4e12 / 5000
p.jobs["cropping"]["final_npts"] = 4096
p.data_process_only_this_data_index = [0]
p.run()
p.save(rdb)

# %% 10/12/2019 - 338-ro-p1-moelle - Rischa :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID112_slaser_R_N=20+_1_longTE_SNR++++_FID48494.dat
/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID114_slaser_R_N=20+_1_longTE_SNR++++_FID48496.dat
/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID118_slaser_R_N=20+_1_longTE_SNR++++_FID48500.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID114_slaser_R_N=20+_1_longTE_SNR++++_FID48496.dat
/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID114_slaser_R_N=20+_1_longTE_SNR++++_FID48496.dat
/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID116_slaser_R_N=20+_1_longTE_SNR++++_FID48498.dat
"""

p.display_legends = """
sLASER 20:1 WS (TWIX)
sLASER 20:1 noWS (TWIX)
sLASER 20:1 IR (TWIX)
"""

# 1) dealing with WS data
p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = [0]
p.run()
p.save(rdb)

# 2) dealing with no WS data: FID modulus
p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = False
p.jobs["channel-combining"]["phasing"] = False
p.data_process_only_this_data_index = [1]
p.analyze_enable = False
p.run()
p.save(rdb)

# %% 28/01/2019 - 300-pm-p2-moelle - Pelayo P2 :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = [
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID62_steam_shortTE_SNR+_FID50920.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID62_steam_shortTE_SNR+_FID50920.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID50926.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID50926.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID73_slaser_R_N=5_5+_shortTE_SNR++_FID50931.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID73_slaser_R_N=5_5+_shortTE_SNR++_FID50931.dat"]

p.data_ref_filepaths = [
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID64_steam_shortTE_SNR+_FID50922.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID65_steam_shortTE_SNR+_FID50923.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID69_slaser_R_N=20+_1_longTE_SNR++++_FID50927.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID70_slaser_R_N=20+_1_longTE_SNR++++_FID50928.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID74_slaser_R_N=5_5+_shortTE_SNR++_FID50932.dat",
    "/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID75_slaser_R_N=5_5+_shortTE_SNR++_FID50933.dat"]

p.display_legends = """
STEAM IR (REF with OVS)
STEAM IR (REF without OVS)
sLASER 20:1 (REF with OVS)
sLASER 20:1 (REF without OVS)
sLASER IR 5:5 (REF with OVS)
sLASER IR 5:5 (REF without OVS)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 06/02/2019 - 349-ap-p1-moelle - Ahmad Fajar
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = [
    "/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID313_slaser_R_N=10_2_longTE_SNR+++_FID51947.dat",
    "/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID313_slaser_R_N=10_2_longTE_SNR+++_FID51947.dat"]

p.data_ref_filepaths = [
    "/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID314_slaser_R_N=10_2_longTE_SNR+++_FID51948.dat",
    "/home/tangir/crmbm/acq_twix/349-ap-p1-moelle/meas_MID315_slaser_R_N=10_2_longTE_SNR+++_FID51949.dat"]

p.display_legends = """
sLASER 10:2 (REF with OVS)
sLASER 10:2 (REF without OVS)
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.run()
p.save(rdb)

# %% 24/02/2019 - 355-st-p1-moelle - Steven
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/355-st-p1-moelle/meas_MID164_slaser_R_N=20+_1_longTE_SNR++++_FID53261.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/355-st-p1-moelle/meas_MID166_slaser_R_N=20+_1_longTE_SNR++++_FID53263.dat"]

p.display_legends = """
sLASER TE=52ms
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["auto"] = True
p.jobs["data-rejecting"]["ranges"]["chemical shift (ppm)"] = 0.17

p.run()
p.save(rdb)

# %% 04/03/2020 - 304-ka-p2-moelle - Karen P2
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID177_slaser_R_N=20+_1_longTE_SNR++++_FID53952.dat",
                    "/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID185_slaser_R_N=20+_1_longTE_SNR++++_FID53960.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID178_slaser_R_N=20+_1_longTE_SNR++++_FID53953.dat",
                        "/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID187_slaser_R_N=20+_1_longTE_SNR++++_FID53962.dat"]

p.data_physio_filepaths = []

p.display_legends = ["1st try (30-40Hz water LW)",
                     "2nd try (25Hz water LW)"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.run()
p.save(rdb)

# %% 30/05/2020 - 311-sl-p2-moelle - Simon P2
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID152_slaser_R_N=20+_1_longTE_SNR++++_FID56036.dat",
                    "/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID156_slaser_R_N=20+_1_longTE_SNR++++_FID56040.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID153_slaser_R_N=20+_1_longTE_SNR++++_FID56037.dat",
                        "/home/tangir/crmbm/acq_twix/311-sl-p2-moelle/meas_MID157_slaser_R_N=20+_1_longTE_SNR++++_FID56041.dat"]

p.data_physio_filepaths = []

p.display_legends = ["sLASER 20/1 NA=128 trig",
                     "sLASER 20/1 NA=128 notrig"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.run()
p.save(rdb)

# %% 09/06/2020 - 336-nb-p2-moelle - Naouelle P2
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID98_slaser_R_N=20+_1_longTE_SNR++++_FID57305.dat",
                    "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID103_slaser_R_N=20+_1_longTE_SNR++++_FID57310.dat",
                    "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID122_slaser_R_N=10_2_longTE_SNR+++_FID57329.dat",
                    "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID127_slaser_R_N=5_5+_shortTE_SNR++_FID57334.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID100_slaser_R_N=20+_1_longTE_SNR++++_FID57307.dat",
                        "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID105_slaser_R_N=20+_1_longTE_SNR++++_FID57312.dat",
                        "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID124_slaser_R_N=10_2_longTE_SNR+++_FID57331.dat",
                        "/home/tangir/crmbm/acq_twix/336-nb-p2-moelle/meas_MID129_slaser_R_N=5_5+_shortTE_SNR++_FID57336.dat"]

p.data_physio_filepaths = []

p.display_legends = ["sLASER 20/1 NA=128 trig",
                     "sLASER 20/1 NA=64 notrig",
                     "sLASER 10/2 NA=64 notrig",
                     "sLASER 5/1 NA=64 notrig"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.data_process_only_this_data_index = []
p.run()
p.save(rdb)

# %% 11/06/2020 - 319-fc-p2-moelle - Fernando P2
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/319-fc-p2-moelle/meas_MID72_slaser_R_N=10_2_longTE_SNR+++_FID57445.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/319-fc-p2-moelle/meas_MID74_slaser_R_N=10_2_longTE_SNR+++_FID57447.dat"]

p.data_physio_filepaths = []

p.display_legends = ["sLASER 10/2 NA=64 trig"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.run()
p.save(rdb)

# %% 15/06/2020 - 313-ft-p2-moelle - Fransiska P2
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/313-ft-p2-moelle/meas_MID239_slaser_R_N=20+_1_longTE_SNR++++_FID57752.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/313-ft-p2-moelle/meas_MID240_slaser_R_N=20+_1_longTE_SNR++++_FID57753.dat"]

p.data_physio_filepaths = []

p.display_legends = ["sLASER 20/1"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.run()
p.save(rdb)

# %% 19/06/2020 - 333-sc-p2-moelle - Shirley P2
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/333-sc-p2-moelle/meas_MID180_slaser_R_N=20+_1_longTE_SNR++++_FID58587.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/333-sc-p2-moelle/meas_MID181_slaser_R_N=20+_1_longTE_SNR++++_FID58588.dat"]

p.data_physio_filepaths = []

p.display_legends = ["sLASER 20/1"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                #p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.run()
p.save(rdb)

# %% 25/06/2020 - 314-yt-p2-moelle - Yolanda P2
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/314-yt-p2-moelle/meas_MID98_slaser_R_N=20+_1_longTE_SNR++++_FID59064.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/314-yt-p2-moelle/meas_MID99_slaser_R_N=20+_1_longTE_SNR++++_FID59065.dat"]

p.data_physio_filepaths = []

p.display_legends = ["sLASER 20/1 NA=128"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                p.jobs["noise-estimation"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["apodizing"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = 50.0
p.jobs["data-rejecting"]["auto_method"] = reco.data_rejection_method.AUTO_LINEWIDTH
p.run()
p.save(rdb)


