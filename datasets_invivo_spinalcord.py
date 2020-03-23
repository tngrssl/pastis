#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vivo data.

@author: Tangi Roussel
"""
# %% init

from __future__ import division
import matplotlib.pylab as plt
import mrs.metabase as xxx
import mrs.reco as reco
import numpy as np

from IPython import get_ipython
import warnings
warnings.filterwarnings("ignore", ".*GUI is implemented*")
plt.close('all')
get_ipython().magic('clear')

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - STEAM #1 :(

get_ipython().magic('clear')
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

p.phase_enable = False
p.data_concatenate = True
p.phase_display = False
p.remove_water_enable = False
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI,
    xxx.m_Tau]

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - sLASER #1 :(

get_ipython().magic('clear')
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

p.data_concatenate = True
p.phase_display = False
p.remove_water_enable = False
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI, xxx.m_Tau]

# %% 15/03/2019 - 291-vs-moelle-spectro-p1 - sLASER #2 :(

get_ipython().magic('clear')
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

p.phase_display = False
p.remove_water_enable = False
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI, xxx.m_Tau]

# %% 20/06/2019 - 296_ym_p1_brainmoelle - brain - sLASER DICOM and TWIX :)

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0014_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0012_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID79_slaser_R_N=20+_1_longTE_SNR++++_FID33878.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID81_slaser_R_N=5_5+_shortTE_SNR++_FID33880.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0013_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0009_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID80_slaser_R_N=5_5+_shortTE_SNR++_FID33879.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat
"""

p.display_legends = """
brain - sLASER R:N=5:5 (DICOM)
brain - sLASER R:N=25:1 (DICOM)
brain - sLASER R:N=25:1 (TWIX)
brain - sLASER R:N=5:5 (TWIX)
brain - sLASER R:N=25:1 REF (TWIX)
"""

p.display_amp_factor_list = [1, 1, 1e8 * 50800 / 26300 * 59000 / 75000, 1e8 * 50800 / 26300 * 59000 / 75000, 1e8 * 50800 / 26300 * 59000 / 75000]

p.recombine_phasing = False
p.apodize_enable = True
p.apodize_damping_hz = 3
p.remove_water_enable = False

p.analyse_and_reject_enable = False
p.analyse_and_reject_min = [-100, 3, -1, -3.14]
p.analyse_and_reject_max = [100, 50, 1, 3.14]
p.analyse_and_reject_auto = False

p.analyse_snr_evol = True
p.analyse_linewidth_evol = True

p.average_na = 4
p.fid_modulus = False
p.remove_water_enable = False
p.remove_water_hsvd_components = 10
p.remove_water_hsvd_range = [4.6, 4.8]

p.data_process_only_this_data_index = [2]
p.display_range_ppm = [1, 5]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI, xxx.m_Tau]

# %% 20/06/2019 - 296_ym_p1_brainmoelle - brain - STEAM DICOM and TWIX :) but bad phase

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0008_steam-shortte-snr/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID71_steam_shortTE_SNR+_FID33870.dat
"""
p.data_ref_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0006_steam-shortte-snr/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID69_steam_shortTE_SNR+_FID33868.dat
"""

p.display_legends = """
brain - STEAM (DICOM)
brain - STEAM (TWIX)
"""

p.display_amp_factor_list = [1, 1e8 * 50800 / 26300 * 59000 / 75000]

p.phase_enable = True
p.phase_POI_range_ppm = [4, 6]
p.phase_weak_ws_mode = True
p.realign_enable = True
p.realign_moving_averages = 5
p.realign_POI_range_ppm = [1.5, 2.5]
p.apodize_enable = False
p.apodize_damping_hz = 5
p.remove_water_enable = False
p.analyse_linewidth_evol = True
p.analyse_snr_evol = True
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI, xxx.m_Tau]

# %% 25/06/2019 - 296_ym_p1_brainmoelle - brain - in vivo sLASER TWIX FID modulus tests to compare with conventional MRS

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat
"""

p.display_legends = """
brain - sLASER FIDmod VAPOR FA=0deg
"""

p.display_amp_factor_list = []

p.fid_modulus = True
p.phase_enable = False
p.recombine_phasing = False
p.apodize_enable = False
p.apodize_damping_hz = 5
p.remove_water_enable = True
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI, xxx.m_Tau]

# %% 26/06/2019 - 296_ym_p1_brainmoelle - STEAM and sLASER :)

get_ipython().magic('clear')
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
DCM steam
DCM sLASER R:N=20:1
TWIX steam
TWIX sLASER R:N=20:1
"""

p.display_amp_factor_list = [1 / 25000, 6 / 25000, 1 / 7.3e-5, 1 / 7.3e-5]

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 0, -0.04, -3.14]
p.analyse_and_reject_max = [+100, 20, +0.04, +3.14]

p.apodize_enable = True
p.apodize_damping_hz = 10

p.data_process_only_this_data_index = [3]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 16/07/2019 - 300-pm-p1-moelle - sLASER trig or not, DICOM and TWIX :)

get_ipython().magic('clear')
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

p.display_amp_factor_list = [1e-8, 1e-8, 1, 1]

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 1, -0.5, -1.5]
p.analyse_and_reject_max = [+200, 25, +0.5, +1.5]
p.analyse_and_reject_auto = False

p.apodize_enable = True
p.apodize_damping_hz = 5

p.data_process_only_this_data_index = [3]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 14/08/2019 - 304-ka-p1-moelle - sLASER, bad day :(

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/304-ka-p1-moelle/20190814/01_0013_slaser-r-n/original-primary_e09_0001.dcm
"""

p.display_legends = """
crappy
"""

p.analyse_and_reject_enable = True

p.apodize_enable = True
p.apodize_damping_hz = 10

p.remove_water_enable = False
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI, xxx.m_Tau]

# %% 21/08/2019 - 307-ap-p1-moelle - sLASER on Ariane :) ok bof

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/307-ap-p1-moelle/20190821/01_0012_slaser-r-n
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
sLASER 20:1 cardiac trig REF
sLASER 20:1 cardiac trig
sLASER 20:1 resp trig
sLASER 20:1 no trig
sLASER 10:1 repos. + resp trig
sLASER 20:1 cardiac trig (TWIX)
sLASER 20:1 resp trig (TWIX)
sLASER 20:1 no trig (TWIX)
sLASER 10:1 repos. + resp trig (TWIX)
"""

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 0, -0.1, -3.14]
p.analyse_and_reject_max = [100, 40, 0.1, 3.14]

p.apodize_enable = True
p.apodize_damping_hz = 5

p.data_process_only_this_data_index = [6]

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 27/08/2019 - 308-rs-p1-moelle - brain - sLASER on Ocha

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0024_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID210_slaser_R_N=20+_1_longTE_SNR++++_FID38955.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0025_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID211_slaser_R_N=20+_1_longTE_SNR++++_FID38956.dat
"""

p.display_legends = """
brain - sLASER 20:1 resp trig
brain - sLASER 20:1 resp trig (TWIX)
"""

p.phase_enable = True
p.recombine_phasing = True
p.realign_enable = True

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 0, -1, -3.14]
p.analyse_and_reject_max = [100, 50, 1, 3.14]

p.apodize_enable = True
p.apodize_damping_hz = 5
p.calibrate_enable = True

p.analyse_linewidth_enable = True
p.analyse_linewidth_magnitude_mode = False

p.analyse_snr_enable = True
p.analyse_snr_magnitude_mode = False

p.remove_water_enable = False

p.display_amp_factor_list = [1, 1e8]
p.data_process_only_this_data_index = [1]

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI, xxx.m_Tau, xxx.m_Glu]

# %% 27/08/2019 - 308-rs-p1-moelle - sLASER on Ocha :(

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 0, -1, -3.14]
p.analyse_and_reject_max = [100, 50, 1, 3.14]

p.apodize_enable = True
p.apodize_damping_hz = 5

p.data_process_only_this_data_index = [1]

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 29/08/2019 - 310-mg-p1-moelle - Maxime :(

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 0, -1, -3.14]
p.analyse_and_reject_max = [100, 50, 1, 3.14]

p.apodize_enable = True
p.apodize_damping_hz = 15

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 05/09/2019 - 311-sl-p1-moelle - Simon :)

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 0, -1, -3.14]
p.analyse_and_reject_max = [100, 21, 1, 3.14]
p.analyse_and_reject_auto = False

p.apodize_enable = True
p.apodize_damping_hz = 10

p.data_process_only_this_data_index = [1]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 23/09/2019 - 313-ft-p1-moelle - Fransiska :(

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-30, 0, -0.1, -0.5]
p.analyse_and_reject_max = [30, 35, 0.1, 0.5]

p.apodize_enable = True
p.apodize_damping_hz = 20

p.analyse_snr_n_range_ppm = [-4, -3]  # fat !

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 25/09/2019 - 314-yt-p1-moelle - Yolanda :)))

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID81_slaser_R_N=20+_1_longTE_SNR++++_FID41679.dat
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID83_slaser_R_N=20+_1_longTE_SNR++++_FID41681.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0010_slaser-r-n
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID88_slaser_R_N=5_5+_shortTE_SNR++_FID41686.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0012_slaser-r-n
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID81_slaser_R_N=20+_1_longTE_SNR++++_FID41679.dat
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID81_slaser_R_N=20+_1_longTE_SNR++++_FID41679.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0009_slaser-r-n
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID86_slaser_R_N=5_5+_shortTE_SNR++_FID41684.dat
/home/tangir/crmbm/acq/314-yt-p1-moelle/20190925/01_0011_slaser-r-n
"""

p.data_physio_filepaths = """

/home/tangir/crmbm/acq_physio/314_YT_P1_MOELLE_2.resp



"""

p.display_legends = """
sLASER 20:1 REF (TWIX)
sLASER 20:1 (TWIX)
sLASER 20:1 (DICOM)
sLASER 5:5 (TWIX)
sLASER 5:5 (DICOM)
"""

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-50, 0, -0.1, -0.5]
p.analyse_and_reject_max = [50, 20, 0.1, 0.5]

p.apodize_enable = True
p.apodize_damping_hz = 10

p.remove_water_enable = False
p.remove_water_hsvd_components = 5
p.remove_water_hsvd_range = [4.6, 4.8]

p.data_process_only_this_data_index = [1]
s, s_ref = p.run_pipeline_std()

# let's assume we were using super adiabatic pulses (important for simulation/fit)
s.sequence.pulse_rfc_r = 50.0

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]
fit_metabolite_threshold = 0.75  # mmol/kg

# %% 03/10/2019 - 316-ap-p1-moelle - Anissa :)

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_min = [-100, 0, -0.1, -0.9]
p.analyse_and_reject_max = [100, 25, 0.1, 0.9]

p.realign_moving_averages = 1
p.realign_POI_range_ppm = [4.5, 4.8]

p.analyse_linewidth_range_ppm = [4.5, 4.8]

p.apodize_enable = True
p.apodize_damping_hz = 10

p.data_process_only_this_data_index = [0]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 17/10/2019 - 319-fc-p1-moelle - Fernando :)

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_moving_averages = 3
p.analyse_and_reject_min = [-100, 0, -0.1, -1.5]
p.analyse_and_reject_max = [100, 13, 0.1, 1.5]
p.analyse_and_reject_auto = False

p.apodize_enable = True
p.apodize_damping_hz = 10

p.calibrate_POI_range_ppm = [1.8, 2.3]
p.analyse_snr_range_ppm = [1.8, 2.2]
p.analyse_snr_n_range_ppm = [-4, -3]

p.data_process_only_this_data_index = [0]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 05/11/2019 - 328-af-p1-moelle - Anne :(

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-100, 0, -0.1, -1.5]
p.analyse_and_reject_max = [100, 40, 0.1, 1.5]
p.analyse_and_reject_auto = True

p.apodize_enable = True
p.apodize_damping_hz = 10

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 08/11/2019 - 329-pi-p1-moelle - Pujalina

get_ipython().magic('clear')
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

p.analyse_and_reject_enable = True
p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-100, 0, -0.1, -1.5]
p.analyse_and_reject_max = [100, 60, 0.1, 1.5]
p.analyse_and_reject_auto = False

p.apodize_enable = True
p.apodize_damping_hz = 10

p.calibrate_POI_range_ppm = [1.8, 2.3]
p.analyse_snr_range_ppm = [1.8, 2.2]
p.analyse_snr_n_range_ppm = [-4, -3]

p.data_process_only_this_data_index = [0]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 26/11/2019 - 333-sc-p1-moelle - Shirley :(((

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/333-sc-p1-moelle/meas_MID123_slaser_R_N=20+_1_longTE_SNR++++_FID47359.dat
/home/tangir/crmbm/acq_twix/333-sc-p1-moelle/meas_MID126_slaser_R_N=20+_1_longTE_SNR++++_FID47362.dat
"""

p.data_ref_filepaths = """"""

p.display_legends = """
sLASER 20:1 (TWIX)
sLASER 20:1 IR (TWIX)
"""

p.analyse_and_reject_enable = False
p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-100, 0, -0.1, -1.5]
p.analyse_and_reject_max = [100, 60, 0.1, 1.5]
p.analyse_and_reject_auto = False

p.apodize_enable = True
p.apodize_damping_hz = 10

p.analyse_linewidth_range_ppm = [4, 6]
p.analyse_snr_s_range_ppm = [4, 6]

p.data_process_only_this_data_index = [1]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 09/12/2019 - 336-nb-p1-moelle - Naouelle :)

get_ipython().magic('clear')
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

p.display_amp_factor_list = [1e8, 1e8, 1, 1]

p.analyse_and_reject_enable = False
p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-100, 0, -0.1, -1.5]
p.analyse_and_reject_max = [100, 60, 0.1, 1.5]
p.analyse_and_reject_auto = False

p.apodize_enable = True
p.apodize_damping_hz = 10

p.analyse_linewidth_range_ppm = [4, 6]
p.analyse_snr_s_range_ppm = [4, 6]

p.data_process_only_this_data_index = [1]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 10/12/2019 - 338-ro-p1-moelle - Rischa :)

get_ipython().magic('clear')
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

p.analyse_snr_enable = False
p.analyse_linewidth_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 10

p.analyse_linewidth_range_ppm = [4, 6]
p.analyse_snr_s_range_ppm = [4, 6]

p.data_process_only_this_data_index = [0]
s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit + [
    xxx.m_Cho_CH2,
    xxx.m_Cr_CH2,
    xxx.m_NAA_CH2,
    xxx.m_mI]

# %% 23/01/2019 - 347-re-p1-moelle - Renaud - brain - MM

get_ipython().magic('clear')
# plt.close("all")

p = reco.pipeline()

p.data_filepaths = [
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID216_slaser_R_N=10_2_longTE_SNR+++_FID50575.dat",
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID221_slaser_R_N=10_2_longTE_SNR+++_FID50580.dat",
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID224_steam_shortTE_SNR+_FID50583.dat"]

p.data_ref_filepaths = [
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID218_slaser_R_N=10_2_longTE_SNR+++_FID50577.dat",
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID222_slaser_R_N=10_2_longTE_SNR+++_FID50581.dat",
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID226_steam_shortTE_SNR+_FID50585.dat"]

p.display_legends = """
sLASER 10:2 IR TE=40ms (TWIX)
sLASER 10:2 IR TE=30ms (TWIX)
STEAM IR TE=3ms (TWIX)
"""

p.phase_offset = np.pi
p.phase_display = False

p.recombine_phasing = True

p.apodize_enable = True
p.apodize_damping_hz = 10.0
p.calibrate_POI_range_ppm = [4, 5]
p.calibrate_POI_true_ppm = 4.7
p.realign_enable = False

p.analyse_linewidth_range_ppm = [4.0, 5.0]
p.analyse_snr_s_range_ppm = [4.0, 5.0]

p.data_process_only_this_data_index = [1]
s_mm, s_mm_ref = p.run_pipeline_std()

# %% 23/01/2019 - 347-re-p1-moelle - Renaud - brain - Metabolites :)
get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = [
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID228_slaser_R_N=10_2_longTE_SNR+++_FID50587.dat"]

p.data_ref_filepaths = [
    "/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID229_slaser_R_N=10_2_longTE_SNR+++_FID50588.dat"]

p.display_legends = """
sLASER 10:2 TE=30ms (TWIX)
"""

p.phase_offset = np.pi
p.phase_display = False

p.recombine_phasing = True

p.apodize_enable = True
p.apodize_damping_hz = 10.0
p.calibrate_POI_range_ppm = [4, 5]
p.calibrate_POI_true_ppm = 4.7
p.realign_enable = False

p.analyse_linewidth_range_ppm = [4.0, 5.0]
p.analyse_snr_s_range_ppm = [4.0, 5.0]

s_mb, s_mb_ref = p.run_pipeline_std()

s = s_mb - s_mm
s_ref = s_mb_ref

plt.close('all')
s_mm.display_spectrum()
s_mb.display_spectrum()
s.display_spectrum()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit
fit_metabolite_threshold = 0.5
fit_include_MMs = False

# %% 28/01/2019 - 300-pm-p2-moelle - Pelayo :)
get_ipython().magic('clear')
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

p.phase_offset = np.pi
p.phase_display = False

p.recombine_phasing = True

p.apodize_enable = True
p.apodize_damping_hz = 10.0
p.calibrate_POI_range_ppm = [4, 5]
p.calibrate_POI_true_ppm = 4.7
p.realign_enable = False

p.analyse_linewidth_range_ppm = [4.0, 5.0]
p.analyse_snr_s_range_ppm = [4.0, 5.0]

p.data_process_only_this_data_index = [2]

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit
fit_metabolite_threshold = 0.5
fit_include_MMs = False

# %% 06/02/2019 - 349-ap-p1-moelle - Ahmad Fajar :s
get_ipython().magic('clear')
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

p.phase_offset = np.pi
p.phase_display = False

p.recombine_phasing = True

p.analyse_and_reject_enable = True

p.apodize_enable = True
p.apodize_damping_hz = 10.0
p.calibrate_POI_range_ppm = [4, 5]
p.calibrate_POI_true_ppm = 4.7
p.realign_enable = False

p.analyse_linewidth_range_ppm = [4.0, 5.0]
p.analyse_snr_s_range_ppm = [4.0, 5.0]

s, s_ref = p.run_pipeline_std()

fit_metabolites_prefit = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
fit_metabolites = fit_metabolites_prefit
fit_metabolite_threshold = 0.5
fit_include_MMs = False


# %% 24/02/2019 - 355-st-p1-moelle - Steven
get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/355-st-p1-moelle/meas_MID164_slaser_R_N=20+_1_longTE_SNR++++_FID53261.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/355-st-p1-moelle/meas_MID166_slaser_R_N=20+_1_longTE_SNR++++_FID53263.dat"]

p.display_legends = """
sLASER TE=52ms
"""

p.analyse_and_reject_enable = True
p.analyse_and_reject_auto = False

p.apodize_enable = True
p.apodize_damping_hz = 10.0

s, s_ref = p.run_pipeline_std()

# %% 04/03/2020 - 304-ka-p2-moelle - Karen :((( very weird artefact

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID177_slaser_R_N=20+_1_longTE_SNR++++_FID53952.dat",
                    "/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID185_slaser_R_N=20+_1_longTE_SNR++++_FID53960.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID178_slaser_R_N=20+_1_longTE_SNR++++_FID53953.dat",
                        "/home/tangir/crmbm/acq_twix/304-ka-p2-moelle/meas_MID187_slaser_R_N=20+_1_longTE_SNR++++_FID53962.dat"]

p.data_physio_filepaths = []

p.display_legends = ["1st try (30-40Hz water LW)",
                     "2nd try (25Hz water LW)"]

p.phase_enable = False
p.recombine_enable = False
p.zerofill_enable = False
p.realign_enable = False
p.average_enable = False
p.apodize_enable = False
p.calibrate_enable = False
p.display_enable = False
p.analyse_snr_enable = False
p.analyse_linewidth_enable = False
p.analyse_snr_evol = False
p.analyse_linewidth_evol = False

p.data_process_only_this_data_index = [0]
s, s_ref = p.run_pipeline_std()

# average by channel
s_chan = np.mean(s, axis=0)

sw_list = []
sfw_list = []
for c in range(8):
    s_chan[c, :].display_spectrum(1, str(c), [2, 6], 1.0, 0.0, True)
    sw = s_chan[c, :].analyse_snr([4.5, 4.8], [-1, 0], 'real water', False, True, True, [2, 6])
    sfw = s_chan[c, :].analyse_snr([5, 5.5], [-1, 0], 'fake water', False, True, True, [2, 6])
    sw_list.append(sw)
    sfw_list.append(sfw)


sw_list = np.array(sw_list)
sfw_list = np.array(sfw_list)

sw_list_norm = sw_list / np.max(sw_list)
sfw_list_norm = sfw_list / np.max(sw_list)

plt.figure()
plt.plot(np.arange(8), sw_list_norm, label="real water residue from VOI")
plt.plot(np.arange(8), sfw_list_norm, label="some water selected somewhere?")
plt.legend()

plt.figure()
plt.plot(np.arange(8), sw_list/sfw_list, label="real water / fake water")