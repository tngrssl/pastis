#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vivo data.

@author: Tangi Roussel
"""

from __future__ import division
import matplotlib.pylab as plt
import mrs.metabase as xxx
import mrs.sim as sim
import mrs.reco as reco

from IPython import get_ipython
import warnings
warnings.filterwarnings("ignore", ".*GUI is implemented*")
plt.close('all')
get_ipython().magic('clear')

# %% checking what twix could be fitted

get_ipython().magic('clear')

plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/coeur/meas_MID135_slaser_R_N=20+_1_longTE_WSAT_FID47083.dat
T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/coeur/meas_MID144_slaser_R_N=20+_1_longTE_WSAT_FID47092.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/coeur/meas_MID142_steam_shortTE_WSAT_FID47235.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/coeur/meas_MID144_slaser_R_N=20+_1_WSAT_FID47237.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/coeur/meas_MID146_steam_ULTRAshortTE_WSAT_FID47239.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID58_slaser_R_N=20+_1_WSAT_FID47151.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID60_steam_shortTE_WSAT_FID47153.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID61_slaser_R_N=20+_1_WSAT_FID47154.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID65_slaser_R_N=20+_1_WSAT_FID47158.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID69_steam_shortTE_WSAT_FID47162.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID71_slaser_R_N=20+_1_WSAT_REPRO_FID47164.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID72_steam_shortTE_WSAT_REPRO_FID47165.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/coeur/meas_MID94_slaser_R_N=20+_1_longTE_WSAT_FID47730.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/coeur/meas_MID95_slaser_R_N=20+_1_longTE_WSAT_FID47731.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/coeur/meas_MID96_slaser_R_N=20+_1_longTE_WSAT_FID47732.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/coeur/meas_MID97_slaser_R_N=20+_1_longTE_WSAT_FID47733.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID27_slaser_R_N=20+_1_WSAT_FID47663.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID31_steam_shortTE_WSAT_FID47667.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID37_steam_shortTE_WSAT_FID47673.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID38_steam_shortTE_WSAT_96NEX_FID47674.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID40_slaser_R_N=20+_1_WSAT_FID47676.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID234_steam_shortTE_WSAT_PULS_FID48364.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID241_slaser_WSAT_PULS_FID48371.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID246_slaser_WSAT_PULS_FID48376.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID250_slaser_WSAT_PULS_10NEX_FID48380.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID252_slaser_WSAT_PULS_5NEX_FID48382.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID254_slaser_WSAT_PULS_5NEX_512_FID48384.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID256_slaser_WSAT_PULS_NOWSAT_FID48386.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID161_steam_shortTE_WSAT_FID48291.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID165_steam_shortTE_WSAT_FID48295.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID170_slaser_R_N=20+_1_WSAT_FID48300.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID176_steam_shortTE_WSAT_FID48306.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID185_slaser_R_N=20+_1_WSAT2_FID48315.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/coeur/meas_MID132_slaser_WSAT_PULS_5NEX_512_FID48742.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/coeur/meas_MID141_slaser2_WSAT_ECG_5NEX_512_FID48751.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/muscle/meas_MID81_steam_shortTE_WSAT_FID48691.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/muscle/meas_MID82_steam_shortTE_WSAT_FID48692.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/muscle/meas_MID85_slaser_R_N=20+_1_WSAT2_FID48695.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/muscle/meas_MID86_slaser_R_N=20+_1_WSAT2_FID48696.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/muscle/meas_MID87_slaser_R_N=20+_1_WSAT2_FID48697.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/muscle/meas_MID96_steam_shortTE_WSAT2_FID48706.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7t/muscle/meas_MID97_steam_shortTE_WSAT3_2048_FID48707.dat
"""

p.data_ref_filepaths = """"""

p.data_physio_filepaths = """"""

p.display_legends = """
191122_MID135_slaser_R_N=20+_1_longTE_WSAT_FID47083
191122_MID144_slaser_R_N=20+_1_longTE_WSAT_FID47092
191125_MID142_steam_shortTE_WSAT_FID47235
191125_MID144_slaser_R_N=20+_1_WSAT_FID47237
191125_MID146_steam_ULTRAshortTE_WSAT_FID47239
191125_MID58_slaser_R_N=20+_1_WSAT_FID47151
191125_MID60_steam_shortTE_WSAT_FID47153
191125_MID61_slaser_R_N=20+_1_WSAT_FID47154
191125_MID65_slaser_R_N=20+_1_WSAT_FID47158
191125_MID69_steam_shortTE_WSAT_FID47162
191125_MID71_slaser_R_N=20+_1_WSAT_REPRO_FID47164
191125_MID72_steam_shortTE_WSAT_REPRO_FID47165
191129_MID94_slaser_R_N=20+_1_longTE_WSAT_FID47730
191129_MID95_slaser_R_N=20+_1_longTE_WSAT_FID47731
191129_MID96_slaser_R_N=20+_1_longTE_WSAT_FID47732
191129_MID97_slaser_R_N=20+_1_longTE_WSAT_FID47733
191129_MID27_slaser_R_N=20+_1_WSAT_FID47663
191129_MID31_steam_shortTE_WSAT_FID47667
191129_MID37_steam_shortTE_WSAT_FID47673
191129_MID38_steam_shortTE_WSAT_96NEX_FID47674
191129_MID40_slaser_R_N=20+_1_WSAT_FID47676
191209_MID234_steam_shortTE_WSAT_PULS_FID48364
191209_MID241_slaser_WSAT_PULS_FID48371
191209_MID246_slaser_WSAT_PULS_FID48376
191209_MID250_slaser_WSAT_PULS_10NEX_FID48380
191209_MID252_slaser_WSAT_PULS_5NEX_FID48382
191209_MID254_slaser_WSAT_PULS_5NEX_512_FID48384
191209_MID256_slaser_WSAT_PULS_NOWSAT_FID48386
191209_MID161_steam_shortTE_WSAT_FID48291
191209_MID165_steam_shortTE_WSAT_FID48295
191209_MID170_slaser_R_N=20+_1_WSAT_FID48300
191209_MID176_steam_shortTE_WSAT_FID48306
191209_MID185_slaser_R_N=20+_1_WSAT2_FID48315
191212_MID132_slaser_WSAT_PULS_5NEX_512_FID48742
191212_MID141_slaser2_WSAT_ECG_5NEX_512_FID48751
191212_MID81_steam_shortTE_WSAT_FID48691
191212_MID82_steam_shortTE_WSAT_FID48692
191212_MID85_slaser_R_N=20+_1_WSAT2_FID48695
191212_MID86_slaser_R_N=20+_1_WSAT2_FID48696
191212_MID87_slaser_R_N=20+_1_WSAT2_FID48697
191212_MID96_steam_shortTE_WSAT2_FID48706
191212_MID97_steam_shortTE_WSAT3_2048_FID48707
"""

p.data_coil_nChannels = 1
p.phase_enable = True
p.phase_POI_range_ppm = [1, 2]

p.analyse_and_reject_enable = False
p.analyse_linewidth_enable = False
p.analyse_snr_enable = False

p.realign_enable = False
p.apodize_enable = True
p.apodize_damping_hz = 5
p.calibrate_enable = True
p.calibrate_POI_range_ppm = [1, 2]
p.calibrate_POI_true_ppm = 1.4

p.data_process_only_this_data_index = [7, 8, 9, 10, 11, 18, 19, 20, 30]
s, s_ref = p.run_pipeline_std()


# %% top10 exploitable twix data

get_ipython().magic('clear')
plt.close("all")

# file paths to water-suppressed raw data
data_filepaths = """
T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID112_slaser_R_N=20+_1_longTE_SNR++++_FID47060.dat
T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID114_steam_shortTE_SNR+_FID47062.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID71_slaser_R_N=20+_1_WSAT_REPRO_FID47164.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID72_steam_shortTE_WSAT_REPRO_FID47165.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID38_steam_shortTE_WSAT_96NEX_FID47674.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID40_slaser_R_N=20+_1_WSAT_FID47676.dat
Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0010_slaser-r-n/original-primary_e09_0001.dcm
Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0011_steam-shortte-wsat3/original-primary_e09_0001.dcm
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID170_slaser_R_N=20+_1_WSAT_FID48300.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID176_steam_shortTE_WSAT_FID48306.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID87_slaser_R_N=20+_1_WSAT2_FID48697.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID82_steam_shortTE_WSAT_FID48692.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID112_slaser_R_N=20+_1_WSAT2_FID50231.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID106_steam_shortTE_WSAT_FID50225.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID198_slaser_R_N=20+_1_WSAT2_FID50317.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID194_steam_shortTE_WSAT_FID50313.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID213_slaser_20_1_LONGTE_FID50332.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID247_steam_shortTE_1024_NOOVS_FID48377.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/coeur/meas_MID141_slaser2_WSAT_ECG_5NEX_512_FID48751.dat
T:/users/SSSR/Data/2020_1H_MRS/200206_1HMRS_zb001//meas_MID131_eja_svs_slaser_VAPOR_FID23607.dat
T:/users/SSSR/Data/2020_1H_MRS/200206_1HMRS_zb001//meas_MID150_eja_svs_slaser_VAPOR_FID23626.dat
T:/users/SSSR/Data/2020_1H_MRS/200206_1HMRS_zb001//meas_MID118_SVSSE_16EX_OVS_WSAT_FID23594.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID112_steam_shortTE_WSAT_FID52195.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID135_steam_shortTE_WSAT_FID52218.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID117_slaser_R_N=20+_1_WSAT2_FID52200.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID137_slaser_R_N=20+_1_WSAT2_FID52220.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID161_slaser_R_N=20+_1_WSAT_LONGTE2_FID52244.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID83_steam_shortTE_WSAT_FID52327.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID93_steam_shortTE_WSAT_VOX2_FID52337.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID88_slaser_R_N=20+_1_WSAT2_FID52332.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID96_slaser_R_N=20+_1_WSAT2_FID52340.dat
T:/users/JS/200228_phantoms/meas_MID35_steam_[10]_WSAT_FID53587.dat
T:/users/JS/200228_phantoms/meas_MID30_slaser_[10]_WSAT_FID53582.dat
T:/users/JS/200228_phantoms/meas_MID54_steam_[25]_WSAT_FID53606.dat
T:/users/JS/200228_phantoms/meas_MID59_slaser_[25]_WSAT_FID53611.dat
T:/users/JS/200228_phantoms/meas_MID75_steam_[50]_WSAT_FID53627.dat
T:/users/JS/200228_phantoms/meas_MID80_slaser_[50]_WSAT_FID53632.dat
"""

# file paths to non water-suppressed raw data
data_ref_filepaths = """
T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID116_slaser_R_N=20+_1_longTE_SNR_NOWSAT_FID47064.dat
T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID115_steam_shortTE_SNR+_NOWSAT_FID47063.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID70_slaser_R_N=20+_1_V1_NOWSAT_FID47163.dat
T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID67_steam_shortTE_NOWSAT_FID47160.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID39_steam_shortTE_NOSAT_FID47675.dat
T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID41_slaser_R_N=20+_1_NOSAT_FID47677.dat
Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0010_slaser-r-n/original-primary_e09_0001.dcm
Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0012_steam-shortte-nowsat3/original-primary_e09_0001.dcm
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID172_slaser_R_N=20+_1_NOWSAT_FID48302.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID180_steam_shortTE_NOWSAT_FID48310.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID89_slaser_R_N=20+_1_NOWSAT2OVS_FID48699.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID83_steam_shortTE_NOWSAT_FID48693.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID113_slaser_R_N=20+_1_NOWSAT2_FID50232.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID107_steam_shortTE_NOWSAT_FID50226.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID199_slaser_R_N=20+_1_NOWSAT2_FID50318.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID195_steam_shortTE_NOWSAT_FID50314.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID213_slaser_20_1_LONGTE_FID50332.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID255_steam_1024_NOOVS_NOWSAT_FID48385.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/coeur/meas_MID142_slaser2_NOWSAT_ECG_5NEX_512_FID48752.dat
T:/users/SSSR/Data/2020_1H_MRS/200206_1HMRS_zb001/meas_MID134_eja_svs_slaser_NOVAPOR_FID23610.dat
T:/users/SSSR/Data/2020_1H_MRS/200206_1HMRS_zb001//meas_MID134_eja_svs_slaser_NOVAPOR_FID23610.dat
T:/users/SSSR/Data/2020_1H_MRS/200206_1HMRS_zb001//meas_MID121_SVSSE_16EX_OVS_NOWSAT_FID23597.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID113_steam_shortTE_NOWSAT_FID52196.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID136_steam_shortTE_NOWSAT_FID52219.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID118_slaser_R_N=20+_1_NOWSAT2_FID52201.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID138_slaser_R_N=20+_1_NOWSAT2_FID52221.dat
T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID161_slaser_R_N=20+_1_WSAT_LONGTE2_FID52244.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID84_steam_shortTE_NOWSAT_FID52328.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID94_steam_shortTE_NOWSAT_VOX2_FID52338.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID89_slaser_R_N=20+_1_NOWSAT2_FID52333.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID97_slaser_R_N=20+_1_NOWSAT2_FID52341.dat
T:/users/JS/200228_phantoms/meas_MID36_steam_[10]_NOWSAT_FID53588.dat
T:/users/JS/200228_phantoms/meas_MID31_slaser_[10]_NOWSAT_FID53583.dat
T:/users/JS/200228_phantoms/meas_MID55_steam_[25]_NOWSAT_FID53607.dat
T:/users/JS/200228_phantoms/meas_MID60_slaser_[25]_NOWSAT_FID53612.dat
T:/users/JS/200228_phantoms/meas_MID82_steam_[50]_NOWSAT_FID53634.dat
T:/users/JS/200228_phantoms/meas_MID81_slaser_[50]_NOWSAT_FID53633.dat
"""

# legend captions
display_legends = """
0_191122_331
1_191122_331
2_191125_332
3_191125_332
4_191129_334_steam
5_191129_334_slaser
6_191203_DICOM
7_191203_DICOM
8_191209_337_slaser
9_191209_337_steam
10_191212_339_slaser
11_191212_339_steam
12_200122_345
13_200122_345
14_200122_346
15_200122_346
16_200122_346_longTE
17_191209_coeur
18_191212_coeur
19_200206_1H_MRS_001_slaser07
20_200206_1H_MRS_001_slaser26
21_200206_1H_MRS_001 PRESS
22_351_steam
23_351_steam_vox2
24_351_slaser
25_351_slaser_vox2
26_351_slaser_longTE
27_352_steam
28_352_steam_vox2
29_352_slaser
30_352_slaser_vox2
31_phantom10_steam
32_phantom10_sLASER
33_phantom25_steam
34_phantom25_sLASER
35_phantom50_steam
36_phantom50_sLASER
"""

# the data lists above, which data index should we process (starting from 0)
# [] means process all
index_to_process = 22


# first process the non water suppressed data
p = reco.pipeline()
p.data_filepaths = data_ref_filepaths
p.data_ref_filepaths = """"""
p.data_physio_filepaths = """"""
p.display_legends = display_legends

# we are using a single channel RF coil
p.data_coil_nChannels = 1
# self-phasing the water spectrum using the first point of FID
p.phase_enable = True
p.phase_weak_ws_mode = True
# let's estimate water linewidth to have an idea of the shim
p.analyse_linewidth_enable = True
p.analyse_linewidth_range_ppm = [4, 5.2]
# no SNR estimation
p.analyse_snr_enable = False
# no spectrum realignment
p.realign_enable = False
# let's denoize a little (5Hz exponential apodization)
p.apodize_damping_hz = 5
# let's calibrate ppm scale so that water peak is at 4.7ppm
p.calibrate_POI_range_ppm = [4, 5.2]
p.calibrate_POI_true_ppm = 4.7
# do not display final processed water spectrum, it is boring
p.display_enable = False
# run the process pipeline
p.data_process_only_this_data_index = [index_to_process]
s_ref, tmp = p.run_pipeline_std()

# now process the water-suppressed (WS) data
p = reco.pipeline()
p.data_filepaths = data_filepaths
p.data_ref_filepaths = """"""
p.data_physio_filepaths = """"""
p.display_legends = display_legends
# we are using a single channel RF coil
p.data_coil_nChannels = 1
# reject bad data
p.analyse_and_reject_enable = True
p.analyse_and_reject_POI_range_ppm = [4, 5.5]
p.analyse_and_reject_min = [-100, 0, -0.5, -3.14]
p.analyse_and_reject_max = [100, 500, 0.5, 3.14]
# self-phasing the WS spectrum using the the big lipid peak at ~1.5ppm
p.phase_enable = True
# p.phase_POI_range_ppm = [1, 2]  # find this peak in this ppm range
# linewidth estimation here (Cr)
p.analyse_linewidth_enable = False
p.analyse_linewidth_range_ppm = [2.9, 3.1]
# but let's estimate the SNR of the Cr peak at 3ppm just to have a rough idea of the quality of the data
p.analyse_snr_enable = True
p.analyse_snr_s_range_ppm = [2.95, 3.10]  # signal ppm range
p.analyse_snr_n_range_ppm = [-3, -1]  # noise ppm range
# but let's estimate the SNR of the lipid peak at 1.2ppm just to have a rough idea of the quality of the data
# p.analyse_snr_enable = False
# p.analyse_snr_s_range_ppm = [1, 1.5]  # signal ppm range
# p.analyse_snr_n_range_ppm = [-3, -1]  # noise ppm range
# no spectrum realignment for now
p.realign_enable = False
# let's denoize a little (5Hz exponential apodization)
p.apodize_damping_hz = 15
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.calibrate_POI_range_ppm = [4, 5.5]
p.calibrate_POI_true_ppm = 4.7
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.calibrate_POI_range_ppm = [1, 2]
# p.calibrate_POI_true_ppm = 1.4
# display ppm range
p.display_range_ppm = [0, 5]
# run the process pipeline
p.data_process_only_this_data_index = [index_to_process]
s, tmp = p.run_pipeline_std()

# some fitting parameters
metabolites_fit = [
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Alcar,
    xxx.m_Carni,
    xxx.m_Carno]

lipids_fit = [
    xxx.m_Lip1,
    xxx.m_Lip2,
    xxx.m_Lip3,
    xxx.m_Lip4]

fit_result_csv_filename = "fit_results.csv"
water_concentration = 87700.0

# to run the fit, switch to the pipeline_fit_muscle.py script cell by cell

# %% man phase
s2 = s * np.exp(1j * 0.5)
s2.display_spectrum()



