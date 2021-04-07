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

# %% process data
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.dataset[0]["legend"] = "0_191122_331"
p.dataset[0]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID112_slaser_R_N=20+_1_longTE_SNR++++_FID47060.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID116_slaser_R_N=20+_1_longTE_SNR_NOWSAT_FID47064.dat"]

p.dataset[1]["legend"] = "1_191122_331"
p.dataset[1]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID114_steam_shortTE_SNR+_FID47062.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191122_NUTRI7T/muscle/meas_MID115_steam_shortTE_SNR+_NOWSAT_FID47063.dat"]

p.dataset[2]["legend"] = "2_191125_332_slaser"
p.dataset[2]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID71_slaser_R_N=20+_1_WSAT_REPRO_FID47164.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID70_slaser_R_N=20+_1_V1_NOWSAT_FID47163.dat"]

p.dataset[3]["legend"] = "3_191125_332"
p.dataset[3]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID72_steam_shortTE_WSAT_REPRO_FID47165.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191125_NUTRI7T/muscle/meas_MID67_steam_shortTE_NOWSAT_FID47160.dat"]

p.dataset[4]["legend"] = "4_191129_334_steam"
p.dataset[4]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID38_steam_shortTE_WSAT_96NEX_FID47674.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID39_steam_shortTE_NOSAT_FID47675.dat"]

p.dataset[5]["legend"] = "5_191129_334_slaser"
p.dataset[5]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID40_slaser_R_N=20+_1_WSAT_FID47676.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191129_NUTRI7T/muscle/meas_MID41_slaser_R_N=20+_1_NOSAT_FID47677.dat"]

p.dataset[6]["legend"] = "6_191203_DICOM"
p.dataset[6]["dcm"]["files"] = ["Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0010_slaser-r-n/original-primary_e09_0001.dcm",
                                "Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0010_slaser-r-n/original-primary_e09_0001.dcm"]

p.dataset[7]["legend"] = "7_191203_DICOM"
p.dataset[7]["dcm"]["files"] = ["Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0011_steam-shortte-wsat3/original-primary_e09_0001.dcm",
                                "Y:/data/users/js/magnetom/nutri7t/335-am-p1-coeur/20191203/01_0012_steam-shortte-nowsat3/original-primary_e09_0001.dcm"]

p.dataset[8]["legend"] = "8_191209_337_slaser"
p.dataset[8]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID170_slaser_R_N=20+_1_WSAT_FID48300.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID172_slaser_R_N=20+_1_NOWSAT_FID48302.dat"]

p.dataset[9]["legend"] = "9_191209_337_steam"
p.dataset[9]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID176_steam_shortTE_WSAT_FID48306.dat",
                                "T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/muscle/meas_MID180_steam_shortTE_NOWSAT_FID48310.dat"]

p.dataset[10]["legend"] = "10_191212_339_slaser"
p.dataset[10]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID87_slaser_R_N=20+_1_WSAT2_FID48697.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID89_slaser_R_N=20+_1_NOWSAT2OVS_FID48699.dat"]

p.dataset[11]["legend"] = "11_191212_339_steam"
p.dataset[11]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID82_steam_shortTE_WSAT_FID48692.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/muscle/meas_MID83_steam_shortTE_NOWSAT_FID48693.dat"]

p.dataset[12]["legend"] = "12_200122_345_slaser"
p.dataset[12]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID112_slaser_R_N=20+_1_WSAT2_FID50231.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID113_slaser_R_N=20+_1_NOWSAT2_FID50232.dat"]

p.dataset[13]["legend"] = "13_200122_345_steam"
p.dataset[13]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID106_steam_shortTE_WSAT_FID50225.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID107_steam_shortTE_NOWSAT_FID50226.dat"]

p.dataset[14]["legend"] = "14_200122_346_slaser"
p.dataset[14]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID198_slaser_R_N=20+_1_WSAT2_FID50317.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID199_slaser_R_N=20+_1_NOWSAT2_FID50318.dat"]

p.dataset[15]["legend"] = "15_200122_346_steam"
p.dataset[15]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID194_steam_shortTE_WSAT_FID50313.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID195_steam_shortTE_NOWSAT_FID50314.dat"]

p.dataset[16]["legend"] = "16_200122_346_longTE"
p.dataset[16]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID203_slaser_20_1_LONGTE_FID50322.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID203_slaser_20_1_LONGTE_FID50322.dat"]

p.dataset[17]["legend"] = "17_191209_coeur"
p.dataset[17]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID247_steam_shortTE_1024_NOOVS_FID48377.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID255_steam_1024_NOOVS_NOWSAT_FID48385.dat"]

p.dataset[18]["legend"] = "18_191212_coeur"
p.dataset[18]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/coeur/meas_MID141_slaser2_WSAT_ECG_5NEX_512_FID48751.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/coeur/meas_MID142_slaser2_NOWSAT_ECG_5NEX_512_FID48752.dat"]

p.dataset[19]["legend"] = "19_200206_1H_MRS_001_slaser07"
p.dataset[19]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID131_eja_svs_slaser_VAPOR_FID23607.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID134_eja_svs_slaser_NOVAPOR_FID23610.dat"]

p.dataset[20]["legend"] = "20_200206_1H_MRS_001_slaser26"
p.dataset[20]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID150_eja_svs_slaser_VAPOR_FID23626.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001//meas_MID134_eja_svs_slaser_NOVAPOR_FID23610.dat"]

p.dataset[21]["legend"] = "21_200206_1H_MRS_001 PRESS"
p.dataset[21]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID118_SVSSE_16EX_OVS_WSAT_FID23594.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001//meas_MID121_SVSSE_16EX_OVS_NOWSAT_FID23597.dat"]

p.dataset[22]["legend"] = "22_351_steam"
p.dataset[22]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID112_steam_shortTE_WSAT_FID52195.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID113_steam_shortTE_NOWSAT_FID52196.dat"]

p.dataset[23]["legend"] = "23_351_steam_vox2"
p.dataset[23]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID135_steam_shortTE_WSAT_FID52218.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID136_steam_shortTE_NOWSAT_FID52219.dat"]

p.dataset[24]["legend"] = "24_351_slaser"
p.dataset[24]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID117_slaser_R_N=20+_1_WSAT2_FID52200.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID118_slaser_R_N=20+_1_NOWSAT2_FID52201.dat"]

p.dataset[25]["legend"] = "25_351_slaser_vox2"
p.dataset[25]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID137_slaser_R_N=20+_1_WSAT2_FID52220.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID138_slaser_R_N=20+_1_NOWSAT2_FID52221.dat"]

p.dataset[26]["legend"] = "26_351_slaser_longTE"
p.dataset[26]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID161_slaser_R_N=20+_1_WSAT_LONGTE2_FID52244.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200211_nutri7t351/meas_MID161_slaser_R_N=20+_1_WSAT_LONGTE2_FID52244.dat"]

p.dataset[27]["legend"] = "27_352_steam"
p.dataset[27]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID83_steam_shortTE_WSAT_FID52327.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID84_steam_shortTE_NOWSAT_FID52328.dat"]

p.dataset[28]["legend"] = "28_352_steam_vox2"
p.dataset[28]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID93_steam_shortTE_WSAT_VOX2_FID52337.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID94_steam_shortTE_NOWSAT_VOX2_FID52338.dat"]

p.dataset[29]["legend"] = "29_352_slaser"
p.dataset[29]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID88_slaser_R_N=20+_1_WSAT2_FID52332.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID89_slaser_R_N=20+_1_NOWSAT2_FID52333.dat"]

p.dataset[30]["legend"] = "30_352_slaser_vox2"
p.dataset[30]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID96_slaser_R_N=20+_1_WSAT2_FID52340.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID97_slaser_R_N=20+_1_NOWSAT2_FID52341.dat"]

p.dataset[31]["legend"] = "31_phantom10_steam"
p.dataset[31]["raw"]["files"] = ["T:/users/JS/200228_phantoms/meas_MID35_steam_[10]_WSAT_FID53587.dat",
                                 "T:/users/JS/200228_phantoms/meas_MID36_steam_[10]_NOWSAT_FID53588.dat"]

p.dataset[32]["legend"] = "32_phantom10_sLASER"
p.dataset[32]["raw"]["files"] = ["T:/users/JS/200228_phantoms/meas_MID30_slaser_[10]_WSAT_FID53582.dat",
                                 "T:/users/JS/200228_phantoms/meas_MID31_slaser_[10]_NOWSAT_FID53583.dat"]

p.dataset[33]["legend"] = "33_phantom25_steam"
p.dataset[33]["raw"]["files"] = ["T:/users/JS/200228_phantoms/meas_MID54_steam_[25]_WSAT_FID53606.dat",
                                 "T:/users/JS/200228_phantoms/meas_MID55_steam_[25]_NOWSAT_FID53607.dat"]

p.dataset[34]["legend"] = "34_phantom25_sLASER"
p.dataset[34]["raw"]["files"] = ["T:/users/JS/200228_phantoms/meas_MID59_slaser_[25]_WSAT_FID53611.dat",
                                 "T:/users/JS/200228_phantoms/meas_MID60_slaser_[25]_NOWSAT_FID53612.dat"]

p.dataset[35]["legend"] = "35_phantom50_steam"
p.dataset[35]["raw"]["files"] = ["T:/users/JS/200228_phantoms/meas_MID75_steam_[50]_WSAT_FID53627.dat",
                                 "T:/users/JS/200228_phantoms/meas_MID82_steam_[50]_NOWSAT_FID53634.dat"]

p.dataset[36]["legend"] = "36_phantom50_sLASER"
p.dataset[36]["raw"]["files"] = ["T:/users/JS/200228_phantoms/meas_MID80_slaser_[50]_WSAT_FID53632.dat",
                                 "T:/users/JS/200228_phantoms/meas_MID81_slaser_[50]_NOWSAT_FID53633.dat"]

p.dataset[37]["legend"] = "37_1HMRS_002_slaser1"
p.dataset[37]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID55_eja_svs_slaser_VAPOR_BH3_FID25884.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID50_eja_svs_slaser_NOVAPOR_BH2_FID25879.dat"]

p.dataset[38]["legend"] = "38_1HMRS_002_slaser2"
p.dataset[38]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID49_eja_svs_slaser_VAPOR_BH2_FID25878.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID50_eja_svs_slaser_NOVAPOR_BH2_FID25879.dat"]

p.dataset[39]["legend"] = "39_1HMRS_002_steam1"
p.dataset[39]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID38_eja_svs_steam_FID25867.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID52_eja_svs_steam_NOVAPOR_BH2_FID25881.dat"]

p.dataset[40]["legend"] = "40_1HMRS_002_steam2"
p.dataset[40]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID51_eja_svs_steam_VAPOR_BH2_FID25880.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID52_eja_svs_steam_NOVAPOR_BH2_FID25881.dat"]

p.dataset[41]["legend"] = "41_1HMRS_002_press1"
p.dataset[41]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID44_eja_svs_press_FID25873.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID54_eja_svs_press_NOVAPOR_BH2_FID25883.dat"]

p.dataset[42]["legend"] = "42_1HMRS_002_press2"
p.dataset[42]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID53_eja_svs_press_VAPOR_BH2_FID25882.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID54_eja_svs_press_NOVAPOR_BH2_FID25883.dat"]

p.dataset[43]["legend"] = "43_1HMRS_002_slaserFR"
p.dataset[43]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID46_eja_svs_slaser_VAPOR_FB_FID25875.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID50_eja_svs_slaser_NOVAPOR_BH2_FID25879.dat"]

p.dataset[44]["legend"] = "44_1HMRS_002_steamFR"
p.dataset[44]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID47_eja_svs_steam_VAPO_FB_FID25876.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID52_eja_svs_steam_NOVAPOR_BH2_FID25881.dat"]

p.dataset[45]["legend"] = "45_1HMRS_002_pressFR"
p.dataset[45]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID48_eja_svs_press_VAPOR_FB_FID25877.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/meas_MID54_eja_svs_press_NOVAPOR_BH2_FID25883.dat"]

p.dataset[46]["legend"] = "46_29_352_slaser_longTE"
p.dataset[46]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID108_slaser_20_1_TE200_FID52352.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID108_slaser_20_1_TE200_FID52352.dat"]

p.dataset[47]["legend"] = "47_345_slaser_TE350"
p.dataset[47]["raw"]["files"] = ["T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID128_slaser_20_1_LONGTE_FID50247.dat",
                                 "T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID128_slaser_20_1_LONGTE_FID50247.dat"]

p.dataset[48]["legend"] = "48_358_slaser_coeur_FR"
p.dataset[48]["raw"]["files"] = ["T:/users/JS/2020_coeur_mrs7T/200311_JB_358/meas_MID62_slaser_WSAT_PULS_100NEX_FID54390.dat",
                                 "T:/users/JS/2020_coeur_mrs7T/200311_JB_358/meas_MID63_slaser_WATER_PULS_10NEX_FID54391.dat"]

p.dataset[49]["legend"] = "49_1HMRS_003_slaser1"
p.dataset[49]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID146_eja_svs_slaser_VAPOR_BH_FID27306.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID147_eja_svs_slaser_NOVAPOR_BH_FID27307.dat"]

p.dataset[50]["legend"] = "50_1HMRS_003_slaser2"
p.dataset[50]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID206_eja_svs_slaser_VAPOR_BH_FID27366.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID147_eja_svs_slaser_NOVAPOR_BH_FID27307.dat"]

p.dataset[51]["legend"] = "51_1HMRS_003_steam1"
p.dataset[51]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID199_eja_svs_steam_VAPOR_BH_FID27359.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID149_eja_svs_steam_NOVAPOR_BH_FID27309.dat"]

p.dataset[52]["legend"] = "52_1HMRS_003_steam2"
p.dataset[52]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID208_eja_svs_steam_VAPOR_BH_FID27368.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID149_eja_svs_steam_NOVAPOR_BH_FID27309.dat"]

p.dataset[53]["legend"] = "53_1HMRS_003_press1"
p.dataset[53]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID200_eja_svs_press_VAPOR_BH_FID27360.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID153_eja_svs_press_NOVAPOR_BH_FID27313.dat"]

p.dataset[54]["legend"] = "54_1HMRS_003_press2"
p.dataset[54]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID209_eja_svs_press_VAPOR_BH_FID27369.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID153_eja_svs_press_NOVAPOR_BH_FID27313.dat"]

p.dataset[55]["legend"] = "55_1HMRS_003_slaserFR"
p.dataset[55]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID201_eja_svs_slaser_VAPOR_FB_FID27361.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID147_eja_svs_slaser_NOVAPOR_BH_FID27307.dat"]

p.dataset[56]["legend"] = "56_1HMRS_003_steamFR"
p.dataset[56]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID202_eja_svs_steam_VAPOR_FB_FID27362.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID149_eja_svs_steam_NOVAPOR_BH_FID27309.dat"]

p.dataset[57]["legend"] = "57_1HMRS_003_pressFR"
p.dataset[57]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID203_eja_svs_press_VAPOR_FB_FID27363.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID153_eja_svs_press_NOVAPOR_BH_FID27313.dat"]

p.dataset[58]["legend"] = "58_1HMRS_004_slaser1"
p.dataset[58]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID85_eja_svs_slaser_VAPOR_BH2_FID27450.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID86_eja_svs_slaser_NOVAPOR_BH2_FID27451.dat"]

p.dataset[59]["legend"] = "59_1HMRS_004_slaser2"
p.dataset[59]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID77_eja_svs_slaser_VAPOR_BH2_FID27442.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID78_eja_svs_slaser_NOVAPOR_BH2_FID27443.dat"]

p.dataset[60]["legend"] = "60_1HMRS_004_steam1"
p.dataset[60]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID71_eja_svs_steam_VAPOR_BH_FID27436.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID72_eja_svs_steam_NOVAPOR_BH_FID27437.dat"]

p.dataset[61]["legend"] = "61_1HMRS_004_steam2"
p.dataset[61]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID91_eja_svs_steam_VAPOR_BH_FID27456.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID92_eja_svs_steam_NOVAPOR_BH_FID27457.dat"]

p.dataset[62]["legend"] = "62_1HMRS_004_press1"
p.dataset[62]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID88_eja_svs_press_VAPOR_BH_FID27453.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID89_eja_svs_press_NOVAPOR_BH_FID27454.dat"]

p.dataset[63]["legend"] = "63_1HMRS_004_press2"
p.dataset[63]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID74_eja_svs_press_VAPOR_BH_FID27439.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID75_eja_svs_press_NOVAPOR_BH_FID27440.dat"]

p.dataset[64]["legend"] = "64_1HMRS_004_slaserFR"
p.dataset[64]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID80_eja_svs_slaser_VAPOR_FB_FID27445.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID86_eja_svs_slaser_NOVAPOR_BH2_FID27451.dat"]

p.dataset[65]["legend"] = "65_1HMRS_004_steamFR"
p.dataset[65]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID84_eja_svs_steam_VAPOR_FB_FID27449.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID92_eja_svs_steam_NOVAPOR_BH_FID27457.dat"]

p.dataset[66]["legend"] = "66_1HMRS_004_pressFR"
p.dataset[66]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID82_eja_svs_press_VAPOR_FB_FID27447.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID75_eja_svs_press_NOVAPOR_BH_FID27440.dat"]

p.dataset[67]["legend"] = "67_1HMRS_005_slaser1"
p.dataset[67]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID257_eja_svs_slaser_VAPOR_BH_FID27622.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID238_eja_svs_slaser_NOVAPOR_BH_FID27603.dat"]

p.dataset[68]["legend"] = "68_1HMRS_005_slaser2"
p.dataset[68]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID246_eja_svs_slaser_VAPOR_BH_FID27611.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID220_eja_svs_slaser_NOVAPOR_BH_FID27585.dat"]

p.dataset[69]["legend"] = "69_1HMRS_005_steam1"
p.dataset[69]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID250_eja_svs_steam_VAPOR_BH_FID27615.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID251_eja_svs_steam_NOVAPOR_BH_FID27616.dat"]

p.dataset[70]["legend"] = "70_1HMRS_005_steam2"
p.dataset[70]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID223_eja_svs_steam_VAPOR_BH_FID27588.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID224_eja_svs_steam_NOVAPOR_BH_FID27589.dat"]

p.dataset[71]["legend"] = "71_1HMRS_005_press1"
p.dataset[71]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID226_eja_svs_press_VAPOR_BH_FID27591.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID227_eja_svs_press_NOVAPOR_BH_FID27592.dat"]

p.dataset[72]["legend"] = "72_1HMRS_005_press2"
p.dataset[72]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID253_eja_svs_press_VAPOR_BH_FID27618.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID254_eja_svs_press_NOVAPOR_BH_FID27619.dat"]

p.dataset[73]["legend"] = "73_1HMRS_005_slaserFR"
p.dataset[73]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID259_eja_svs_slaser_VAPOR_FB_FID27624.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID248_eja_svs_slaser_NOVAPOR_BH_FID27613.dat"]

p.dataset[74]["legend"] = "74_1HMRS_005_steamFR"
p.dataset[74]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID233_eja_svs_steam_VAPOR_FB_FID27598.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID251_eja_svs_steam_NOVAPOR_BH_FID27616.dat"]

p.dataset[75]["legend"] = "75_1HMRS_005_pressFR"
p.dataset[75]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID235_eja_svs_press_VAPOR_FB_FID27600.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID254_eja_svs_press_NOVAPOR_BH_FID27619.dat"]

p.dataset[76]["legend"] = "76_1HMRS_006_slaser1"
p.dataset[76]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID287_eja_svs_slaser_VAPOR_BH_FID27652.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID288_eja_svs_slaser_NOVAPOR_BH_FID27653.dat"]

p.dataset[77]["legend"] = "77_1HMRS_006_slaser2"
p.dataset[77]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID311_eja_svs_slaser_VAPOR_BH_FID27676.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID288_eja_svs_slaser_NOVAPOR_BH_FID27653.dat"]

p.dataset[78]["legend"] = "78_1HMRS_006_steam1"
p.dataset[78]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID291_eja_svs_steam_VAPOR_BH_FID27656.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID292_eja_svs_steam_NOVAPOR_BH_FID27657.dat"]

p.dataset[79]["legend"] = "79_1HMRS_006_steam2"
p.dataset[79]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID290_eja_svs_steam_VAPOR_BH_FID27655.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID292_eja_svs_steam_NOVAPOR_BH_FID27657.dat"]

p.dataset[80]["legend"] = "80_1HMRS_006_press1"
p.dataset[80]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID312_eja_svs_press_VAPOR_BH_FID27677.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID295_eja_svs_press_NOVAPOR_BH_FID27660.dat"]

p.dataset[81]["legend"] = "81_1HMRS_006_press2"
p.dataset[81]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID294_eja_svs_press_VAPOR_BH_FID27659.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID295_eja_svs_press_NOVAPOR_BH_FID27660.dat"]

p.dataset[82]["legend"] = "82_1HMRS_006_slaserFR"
p.dataset[82]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID297_eja_svs_slaser_VAPOR_FB_FID27662.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID288_eja_svs_slaser_NOVAPOR_BH_FID27653.dat"]

p.dataset[83]["legend"] = "83_1HMRS_006_steamFR"
p.dataset[83]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID298_eja_svs_steam_VAPOR_FB_FID27663.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID292_eja_svs_steam_NOVAPOR_BH_FID27657.dat"]

p.dataset[84]["legend"] = "84_1HMRS_006_pressFR"
p.dataset[84]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID301_eja_svs_press_VAPOR_FB_FID27666.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID295_eja_svs_press_NOVAPOR_BH_FID27660.dat"]

p.dataset[85]["legend"] = "85_1HMRS_004_steam3"
p.dataset[85]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID57_eja_svs_steam_VAPOR_BH_FID27422.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID92_eja_svs_steam_NOVAPOR_BH_FID27457.dat"]

p.dataset[86]["legend"] = "86_1HMRS_007_slaser1"
p.dataset[86]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID62_eja_svs_slaser_VAPOR_BH_FID29816.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID65_eja_svs_slaser_NOVAPOR_BH_FID29819.dat"]

p.dataset[87]["legend"] = "87_1HMRS_007_slaser2"
p.dataset[87]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID77_eja_svs_slaser_VAPOR_BH_REP_FID29831.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID65_eja_svs_slaser_NOVAPOR_BH_FID29819.dat"]

p.dataset[88]["legend"] = "88_1HMRS_007_steam1"
p.dataset[88]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID66_eja_svs_steam_NOVAPOR_BH_FID29820.dat"]

p.dataset[89]["legend"] = "89_1HMRS_007_steam2"
p.dataset[89]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID66_eja_svs_steam_NOVAPOR_BH_FID29820.dat"]

p.dataset[90]["legend"] = "90_1HMRS_007_press1"
p.dataset[90]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID57_eja_svs_press_VAPOR_BH_FID29811.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID64_eja_svs_press_NOVAPOR_BH_FID29818.dat"]

p.dataset[91]["legend"] = "91_1HMRS_007_press2"
p.dataset[91]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID75_eja_svs_press_VAPOR_BH_REP_FID29829.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID64_eja_svs_press_NOVAPOR_BH_FID29818.dat"]

p.dataset[92]["legend"] = "92_1HMRS_007_slaserFB"
p.dataset[92]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID68_eja_svs_slaser_VAPOR_FB_FID29822.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID65_eja_svs_slaser_NOVAPOR_BH_FID29819.dat"]

p.dataset[93]["legend"] = "93_1HMRS_007_steamFB"
p.dataset[93]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID73_eja_svs_steam_VAPOR_FB_FID29827.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID66_eja_svs_steam_NOVAPOR_BH_FID29820.dat"]

p.dataset[94]["legend"] = "94_1HMRS_007_pressFB"
p.dataset[94]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID72_eja_svs_press_VAPOR_FB_FID29826.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_007js2/meas_MID64_eja_svs_press_NOVAPOR_BH_FID29818.dat"]

p.dataset[95]["legend"] = "95_1HMRS_008_slaser1"
p.dataset[95]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID101_eja_svs_press_VAPOR_BH_FID29855.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID106_eja_svs_slaser_NOVAPOR_BH_FID29860.dat"]

p.dataset[96]["legend"] = "96_1HMRS_008_slaser2"
p.dataset[96]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID129_eja_svs_slaser_VAPOR_BH_FID29883.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID122_eja_svs_slaser_NOVAPOR_BH_FID29876.dat"]

p.dataset[97]["legend"] = "97_1HMRS_008_steam1"
p.dataset[97]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID123_eja_svs_steam_NOVAPOR_BH_FID29877.dat"]

p.dataset[98]["legend"] = "98_1HMRS_008_steam2"
p.dataset[98]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/"]

p.dataset[99]["legend"] = "99_1HMRS_008_press1"
p.dataset[99]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID118_eja_svs_press_VAPOR_BH_FID29872.dat",
                                 "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID121_eja_svs_press_NOVAPOR_BH_FID29875.dat"]

p.dataset[100]["legend"] = "100_1HMRS_008_press2"
p.dataset[100]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID131_eja_svs_press_VAPOR_BH_FID29885.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID121_eja_svs_press_NOVAPOR_BH_FID29875.dat"]

p.dataset[101]["legend"] = "101_1HMRS_008_slaserFB"
p.dataset[101]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID124_eja_svs_slaser_VAPOR_FB_FID29878.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID122_eja_svs_slaser_NOVAPOR_BH_FID29876.dat"]

p.dataset[102]["legend"] = "102_1HMRS_008_steamFB"
p.dataset[102]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID126_eja_svs_steam_VAPOR_FB_FID29880.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/"]

p.dataset[103]["legend"] = "103_1HMRS_008_pressFB"
p.dataset[103]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID128_eja_svs_press_VAPOR_FB_FID29882.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_008sr2/meas_MID121_eja_svs_press_NOVAPOR_BH_FID29875.dat"]

p.dataset[104]["legend"] = "104_1HMRS_009_slaser1"
p.dataset[104]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID134_eja_svs_slaser_VAPOR_BH_FID29305.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID120_eja_svs_slaser_NOVAPOR_FID29291.dat"]

p.dataset[105]["legend"] = "105_1HMRS_009_slaser2"
p.dataset[105]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID139_eja_svs_slaser_VAPOR_BH_FID29310.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID120_eja_svs_slaser_NOVAPOR_FID29291.dat"]

p.dataset[106]["legend"] = "106_1HMRS_009_steam1"
p.dataset[106]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID122_eja_svs_steam_NOVAPOR_FID29293.dat"]

p.dataset[107]["legend"] = "107_1HMRS_009_steam2"
p.dataset[107]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID122_eja_svs_steam_NOVAPOR_FID29293.dat"]

p.dataset[108]["legend"] = "108_1HMRS_009_press1"
p.dataset[108]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID111_eja_svs_press_VAPOR_BH_FID29282.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID123_eja_svs_press_novapor_FID29294.dat"]

p.dataset[109]["legend"] = "109_1HMRS_009_press2"
p.dataset[109]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID140_eja_svs_press_VAPOR_BH_FID29311.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID123_eja_svs_press_novapor_FID29294.dat"]

p.dataset[110]["legend"] = "110_1HMRS_009_slaserFB"
p.dataset[110]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID124_eja_svs_slaser_VAPOR_FB_FID29295.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID120_eja_svs_slaser_NOVAPOR_FID29291.dat"]

p.dataset[111]["legend"] = "111_1HMRS_009_steamFB"
p.dataset[111]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID125_eja_svs_steam_VAPOR_FB_FID29296.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID122_eja_svs_steam_NOVAPOR_FID29293.dat"]

p.dataset[112]["legend"] = "112_1HMRS_009_pressFB"
p.dataset[112]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID127_eja_svs_press_VAPOR_FB_FID29298.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_pv009/meas_MID123_eja_svs_press_novapor_FID29294.dat"]

p.dataset[113]["legend"] = "113_1HMRS_010_slaser1"
p.dataset[113]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID197_eja_svs_slaser_VAPOR_BHVOX2_FID29368.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID191_eja_svs_slaser_NOVAPOR_FID29362.dat"]

p.dataset[114]["legend"] = "114_1HMRS_010_slaser2"
p.dataset[114]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID202_eja_svs_slaser_VAPOR_BHVOX2REPRO_FID29373.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID191_eja_svs_slaser_NOVAPOR_FID29362.dat"]

p.dataset[115]["legend"] = "115_1HMRS_010_steam1"
p.dataset[115]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID193_eja_svs_steam_NOVAPOR_FID29364.dat"]

p.dataset[116]["legend"] = "116_1HMRS_010_steam2"
p.dataset[116]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID193_eja_svs_steam_NOVAPOR_FID29364.dat"]

p.dataset[117]["legend"] = "117_1HMRS_010_press1"
p.dataset[117]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID189_eja_svs_press_VAPOR_BHVOX2_FID29360.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID192_eja_svs_press_NOVAPOR_FID29363.dat"]

p.dataset[118]["legend"] = "118_1HMRS_010_press2"
p.dataset[118]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID203_eja_svs_press_VAPOR_BHVOX2REPRO_FID29374.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID192_eja_svs_press_NOVAPOR_FID29363.dat"]

p.dataset[119]["legend"] = "119_1HMRS_010_slaserFB"
p.dataset[119]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID200_eja_svs_slaser_VAPOR_FB_FID29371.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID191_eja_svs_slaser_NOVAPOR_FID29362.dat"]

p.dataset[120]["legend"] = "120_1HMRS_010_steamFB"
p.dataset[120]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID196_eja_svs_steam_VAPOR_FB_FID29367.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID193_eja_svs_steam_NOVAPOR_FID29364.dat"]

p.dataset[121]["legend"] = "121_1HMRS_010_pressFB"
p.dataset[121]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID195_eja_svs_press_VAPOR_FB_FID29366.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202605_1hmrs_lp010/meas_MID192_eja_svs_press_NOVAPOR_FID29363.dat"]

p.dataset[122]["legend"] = "122_1HMRS_phantom001_slaser"
p.dataset[122]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID61_eja_svs_slaser_NOVAPOR_VOX4_FID28498.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID61_eja_svs_slaser_NOVAPOR_VOX4_FID28498.dat"]

p.dataset[123]["legend"] = "123_1HMRS_phantom001_steam"
p.dataset[123]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID63_eja_svs_steam_NOVAPOR_VOX4_FID28500.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID63_eja_svs_steam_NOVAPOR_VOX4_FID28500.dat"]

p.dataset[124]["legend"] = "124_1HMRS_phantom001_press"
p.dataset[124]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID64_eja_svs_press_NOVAPOR_VOX4_FID28501.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID64_eja_svs_press_NOVAPOR_VOX4_FID28501.dat"]

p.dataset[125]["legend"] = "125_1hmrs_phantom001_steamTE64"
p.dataset[125]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID66_eja_svs_steam_NOVAPOR_VOX4_TE64_FID28503.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID66_eja_svs_steam_NOVAPOR_VOX4_TE64_FID28503.dat"]

p.dataset[126]["legend"] = "126_1hmrs_phantom001_pressTE64"
p.dataset[126]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID65_eja_svs_press_NOVAPOR_VOX4_TE64_FID28502.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/201405_1hmrs_phantom001/meas_MID65_eja_svs_press_NOVAPOR_VOX4_TE64_FID28502.dat"]

p.dataset[127]["legend"] = "127_1hmrs_011_slaser1"
p.dataset[127]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID220_eja_svs_slaser_VAPOR_BHVOX2_FID29974.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID228_eja_svs_slaser_NOVAPOR_BH_FID29982.dat"]

p.dataset[128]["legend"] = "128_1hmrs_011_slaser2"
p.dataset[128]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID223_eja_svs_slaser_VAPOR_BHVOX2_FID29977.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID228_eja_svs_slaser_NOVAPOR_BH_FID29982.dat"]

p.dataset[129]["legend"] = "129_1hmrs_011_steam1"
p.dataset[129]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/"]

p.dataset[130]["legend"] = "130_1hmrs_011_steam2"
p.dataset[130]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/"]

p.dataset[131]["legend"] = "131_1hmrs_011_press1"
p.dataset[131]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID221_eja_svs_press_VAPOR_BH8VOX2_FID29975.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID212_eja_svs_press_NOVAPOR_BH_FID29966.dat"]

p.dataset[132]["legend"] = "132_1hmrs_011_press2"
p.dataset[132]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID224_eja_svs_press_VAPOR_BH8VOX2_FID29978.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID227_eja_svs_press_NOVAPOR_BH_FID29981.dat"]

p.dataset[133]["legend"] = "133_1hmrs_011_slaserFB"
p.dataset[133]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID215_eja_svs_slaser_VAPOR_FB_FID29969.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID213_eja_svs_slaser_NOVAPOR_BH_FID29967.dat"]

p.dataset[134]["legend"] = "134_1hmrs_011_steamFB"
p.dataset[134]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID219_eja_svs_steam_VAPOR_FB_FID29973.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID214_eja_svs_steam_NOVAPOR_BH_FID29968.dat"]

p.dataset[135]["legend"] = "135_1hmrs_011_pressFB"
p.dataset[135]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID217_eja_svs_press_VAPOR_FB_FID29971.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/202905_1hmrs_011tr/meas_MID212_eja_svs_press_NOVAPOR_BH_FID29966.dat"]

p.dataset[136]["legend"] = "136_1hmrs_012_slaser1"
p.dataset[136]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID41_eja_svs_slaser_VAPOR_BH_FID30852.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID64_eja_svs_slaser_NOVAPOR_BH_FID30875.dat"]

p.dataset[137]["legend"] = "137_1hmrs_012_slaser2"
p.dataset[137]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID74_eja_svs_slaser_VAPOR_BH_REPRO_FID30885.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID80_eja_svs_slaser_NOVAPOR_BH_REPRO_FID30891.dat"]
p.dataset[137]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/012-vc/20200605/01_0015_eja-svs-press-vapor-bh/original-primary_e09_0001.dcm",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID80_eja_svs_slaser_NOVAPOR_BH_REPRO_FID30891.dat"]

p.dataset[138]["legend"] = "138_1hmrs_012_steam1"
p.dataset[138]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID63_eja_svs_steam_NOVAPOR_BH_FID30874.dat"]

p.dataset[139]["legend"] = "139_1hmrs_012_steam2"
p.dataset[139]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID63_eja_svs_steam_NOVAPOR_BH_FID30874.dat"]

p.dataset[140]["legend"] = "140_1hmrs_012_press1"
p.dataset[140]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID45_eja_svs_press_VAPOR_BH_FID30856.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID61_eja_svs_press_NOVAPOR_BH_FID30872.dat"]

p.dataset[141]["legend"] = "141_1hmrs_012_press2"
p.dataset[141]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID51_eja_svs_press_VAPOR_BH_FID30862.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID61_eja_svs_press_NOVAPOR_BH_FID30872.dat"]

p.dataset[142]["legend"] = "142_1hmrs_012_slaserFB"
p.dataset[142]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID65_eja_svs_slaser_VAPOR_FB_FID30876.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID64_eja_svs_slaser_NOVAPOR_BH_FID30875.dat"]

p.dataset[143]["legend"] = "143_1hmrs_012_steamFB"
p.dataset[143]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID63_eja_svs_steam_NOVAPOR_BH_FID30874.dat"]

p.dataset[144]["legend"] = "144_1hmrs_012_pressFB"
p.dataset[144]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID67_eja_svs_press_VAPOR_FB_FID30878.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_012vc/meas_MID61_eja_svs_press_NOVAPOR_BH_FID30872.dat"]

p.dataset[145]["legend"] = "145_1hmrs_013_slaser1"
p.dataset[145]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID106_eja_svs_slaser_VAPOR_BH_FID30917.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID116_eja_svs_slaser_NOVAPOR_BH_FID30927.dat"]

p.dataset[146]["legend"] = "146_1hmrs_013_slaser2"
p.dataset[146]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID130_eja_svs_slaser_VAPOR_BH_REPRO2_FID30941.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID132_eja_svs_slaser_NOVAPOR_BH_FID30943.dat"]

p.dataset[147]["legend"] = "147_1hmrs_013_steam1"
p.dataset[147]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/"]

p.dataset[148]["legend"] = "148_1hmrs_013_steam2"
p.dataset[148]["dcm"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/"]

p.dataset[149]["legend"] = "149_1hmrs_013_press1"
p.dataset[149]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID114_eja_svs_press_VAPOR_BH_FID30925.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID115_eja_svs_press_NOVAPOR_BH_FID30926.dat"]

p.dataset[150]["legend"] = "150_1hmrs_013_press2"
p.dataset[150]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID128_eja_svs_press_VAPOR_BH_REPRO_FID30939.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID131_eja_svs_press_NOVAPOR_BH_FID30942.dat"]

p.dataset[151]["legend"] = "151_1hmrs_013_slaserFB"
p.dataset[151]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID120_eja_svs_slaser_VAPOR_FB_FID30931.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID116_eja_svs_slaser_NOVAPOR_BH_FID30927.dat"]

p.dataset[152]["legend"] = "152_1hmrs_013_steamFB"
p.dataset[152]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID124_eja_svs_steam_VAPOR_FB_FID30935.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID129_eja_svs_steam_VAPOR_BH_REPRO_FID30940.dat"]

p.dataset[153]["legend"] = "153_1hmrs_013_pressFB"
p.dataset[153]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID122_eja_svs_press_VAPOR_FB_FID30933.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200506_1hmrs_013zb/meas_MID115_eja_svs_press_NOVAPOR_BH_FID30926.dat"]

p.dataset[154]["legend"] = None
p.dataset[154]["dcm"]["files"] = [None, None]

p.dataset[155]["legend"] = "155_dicom"
p.dataset[155]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/007-js2/20200528/01_0059_eja-svs-press-vapor-bh-rep/original-primary_e09_0001.dcm",
                                  "Y:/data/users/js/verio/1h_mrs/009-pv/20200525/01_0015_eja-svs-slaser-vapor-bh/original-primary_e09_0001.dcm"]

p.dataset[156]["legend"] = "156_1hmrs_014_slaser1"
p.dataset[156]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID99_eja_svs_slaser_VAPOR_BH_FID33268.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID102_eja_svs_slaser_NOVAPOR_BH_FID33271.dat"]

p.dataset[157]["legend"] = "157_1hmrs_014_slaser2"
p.dataset[157]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID108_eja_svs_slaser_VAPOR_BH_REPRO_FID33277.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID114_eja_svs_slaser_NOVAPOR_BH_TO_USE_FID33283.dat"]

p.dataset[158]["legend"] = "158_1hmrs_014_steam1"
p.dataset[158]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID98_eja_svs_steam_VAPOR_BH_FID33267.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID104_eja_svs_steam_NOVAPOR_BH_FID33273.dat"]

p.dataset[159]["legend"] = "159_1hmrs_014_steam2"
p.dataset[159]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID112_eja_svs_steam_VAPOR_BH_FID33281.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID113_eja_svs_steam_NOVAPOR_BH_FID33282.dat"]

p.dataset[160]["legend"] = "160_1hmrs_014_press1"
p.dataset[160]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID97_eja_svs_press_VAPOR_BH_FID33266.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID103_eja_svs_press_NOVAPOR_BH_FID33272.dat"]

p.dataset[161]["legend"] = "161_1hmrs_014_press2"
p.dataset[161]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID109_eja_svs_press_VAPOR_BH_REPRO_FID33278.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID111_eja_svs_press_NOVAPOR_BH_FID33280.dat"]

p.dataset[162]["legend"] = "162_1hmrs_014_slaserFB"
p.dataset[162]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID105_eja_svs_slaser_VAPOR_FB_FID33274.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID102_eja_svs_slaser_NOVAPOR_BH_FID33271.dat"]

p.dataset[163]["legend"] = "163_1hmrs_014_steamFB"
p.dataset[163]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID107_eja_svs_steam_VAPOR_FB_FID33276.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID104_eja_svs_steam_NOVAPOR_BH_FID33273.dat"]

p.dataset[164]["legend"] = "164_1hmrs_014_pressFB"
p.dataset[164]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID106_eja_svs_press_VAPOR_FB_FID33275.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200701_1hmrs_014al/meas_MID103_eja_svs_press_NOVAPOR_BH_FID33272.dat"]

p.dataset[165]["legend"] = "165_1hmrs_015_slaser1"
p.dataset[165]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID48_eja_svs_slaser_VAPOR_BH_FID33500.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID56_eja_svs_slaser_NOVAPOR_BH_FID33508.dat"]

p.dataset[166]["legend"] = "166_1hmrs_015_slaser2"
p.dataset[166]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID53_eja_svs_slaser_VAPOR_BH_FID33505.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID56_eja_svs_slaser_NOVAPOR_BH_FID33508.dat"]

p.dataset[167]["legend"] = "167_1hmrs_015_steam1"
p.dataset[167]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID52_eja_svs_steam_VAPOR_BH_FID33504.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID57_eja_svs_steam_NOVAPOR_BH_FID33509.dat"]

p.dataset[168]["legend"] = "168_1hmrs_015_steam2"
p.dataset[168]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID67_eja_svs_steam_VAPOR_BH_REPRO_FID33519.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID70_eja_svs_steam_NOVAPOR_BH_FID33522.dat"]

p.dataset[169]["legend"] = "169_1hmrs_015_press1"
p.dataset[169]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID50_eja_svs_press_VAPOR_BH_FID33502.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID55_eja_svs_press_NOVAPOR_BH_FID33507.dat"]

p.dataset[170]["legend"] = "170_1hmrs_015_press2"
p.dataset[170]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID66_eja_svs_press_VAPOR_BH_REPRO_FID33518.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID69_eja_svs_press_NOVAPOR_BH_FID33521.dat"]

p.dataset[171]["legend"] = "171_1hmrs_015_slaserFB"
p.dataset[171]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID58_eja_svs_slaser_VAPOR_FB_FID33510.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID56_eja_svs_slaser_NOVAPOR_BH_FID33508.dat"]

p.dataset[172]["legend"] = "172_1hmrs_015_steamFB"
p.dataset[172]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID62_eja_svs_steam_VAPOR_FB_FID33514.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID57_eja_svs_steam_NOVAPOR_BH_FID33509.dat"]

p.dataset[173]["legend"] = "173_1hmrs_015_pressFB"
p.dataset[173]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID60_eja_svs_press_VAPOR_FB_FID33512.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_015fk/meas_MID55_eja_svs_press_NOVAPOR_BH_FID33507.dat"]

p.dataset[174]["legend"] = "174_1hmrs_016_slaserTUBE1"
p.dataset[174]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0008_eja-svs-slaser-vapor-bh/original-primary_e09_0001.dcm",
                                  "Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0010_eja-svs-slaser-novapor-tube1/original-primary_e09_0001.dcm"]
p.dataset[174]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID115_eja_svs_slaser_VAPOR_BH_FID33567.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID116_eja_svs_slaser_NOVAPOR_TUBE1_FID33568.dat"]

p.dataset[175]["legend"] = "175_1hmrs_016_steamTUBE1"
p.dataset[175]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID120_eja_svs_steam_VAPOR_TUBE1_FID33572.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID121_eja_svs_steam_NOVAPOR_TUBE1_FID33573.dat"]
p.dataset[175]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0016_eja-svs-steam-vapor-tube1/original-primary_e09_0001.dcm",
                                  "Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0018_eja-svs-steam-novapor-tube1/original-primary_e09_0001.dcm"]

p.dataset[176]["legend"] = "176_1hmrs_016_pressTUBE1"
p.dataset[176]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID118_eja_svs_press_VAPOR_TUBE1_FID33570.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID119_eja_svs_press_NOVAPOR_TUBE1_FID33571.dat"]
p.dataset[176]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0012_eja-svs-press-vapor-tube1/original-primary_e09_0001.dcm",
                                  "Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0014_eja-svs-press-novapor-tube1/original-primary_e09_0001.dcm"]

p.dataset[177]["legend"] = "177_1hmrs_016_slaserTUBE2"
p.dataset[177]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID129_eja_svs_slaser_VAPOR_TUBE2to_use_FID33581.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID131_eja_svs_slaser_NOVAPOR_TUBE2to_use++_FID33583.dat"]
p.dataset[177]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0024_eja-svs-slaser-vapor-tube2to-use/original-primary_e09_0001.dcm",
                                  "Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0026_eja-svs-slaser-novapor-tube2to-use/original-primary_e09_0001.dcm"]

p.dataset[178]["legend"] = "178_1hmrs_016_steamTUBE2"
p.dataset[178]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID134_eja_svs_steam_VAPOR_TUBE2_to_use_FID33586.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID135_eja_svs_steam_NOVAPOR_TUBE2_to_use_FID33587.dat"]
p.dataset[178]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0034_eja-svs-steam-vapor-tube2-to-use/original-primary_e09_0001.dcm",
                                  "Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0032_eja-svs-press-novapor-tube2-to-use/original-primary_e09_0001.dcm"]

p.dataset[179]["legend"] = "179_1hmrs_016_pressTUBE2"
p.dataset[179]["raw"]["files"] = ["C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID132_eja_svs_press_VAPOR_TUBE2_to_use_FID33584.dat",
                                  "C:/Users/jsourdon/Desktop/2020_1H_MRS/200702_1hmrs_016_ConcCr/meas_MID133_eja_svs_press_NOVAPOR_TUBE2_to_use_FID33585.dat"]
p.dataset[179]["dcm"]["files"] = ["Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0030_eja-svs-press-vapor-tube2-to-use/original-primary_e09_0001.dcm",
                                  "Y:/data/users/js/verio/1h_mrs/016-conc-cr/20200702/01_0036_eja-svs-steam-novapor-tube2-to-use/original-primary_e09_0001.dcm"]

# the data lists above, which data index should we process (starting from 0)
# [] means process all
index_to_process = 178

# --- process the water-suppressed (WS) data ---

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                p.job["channel_combining"],
                p.job["noise_estimation"],
                p.job["zero_filling"],
                p.job["apodizing"],
                #p.job["data_rejecting"],
                #p.job["data_rejecting"],
                #p.job["realigning"],
                p.job["averaging"],
                p.job["cropping"],
                #p.job["water_removal"],
                p.job["calibrating"],
                #p.job["phasing (suspect)"],
                p.job["displaying"]
                ]

p.analyze_enable = True
p.analyze_job_list = [  p.job["channel_combining"],
                        p.job["zero_filling"],
                        #p.job["realigning"],
                        p.job["averaging"],
                        p.job["calibrating"]
                        ]

p.settings["pkl_filepath"] = "/home/tangir/crmbm/acq_db/muscle.pkl"

# self-phasing the WS spectrum using the the big lipid peak at ~1.5ppm (Here you can choose btw Lipid and Cr for phasing)
p.job["phasing"]["POI_range_ppm"] = [4, 5]  # find this peak in this ppm range
p.job["phasing"]["offset"] = 0.0  # manual phase offset
p.job["phasing"]["average_per_channel_mode"] = False
p.job["phasing"]["first_point_fid_mode"] = False
p.job["phasing"]["using_ref_data"] = False
p.job["phasing"]["display"] = False

p.job["channel_combining"]["phasing"] = True
p.job["channel_combining"]["using_ref_data"] = True

# reject bad data
p.job["data_rejecting"]["moving_averages"] = 1
p.job["data_rejecting"]["POI_range_ppm"] = [4, 5.2]
# rejection ranges
#p.job["data_rejecting"]["ranges"]["amplitude (%)"] = 50  # do not reject on amplitude changes
p.job["data_rejecting"]["ranges"]["linewidth (Hz)"] = [-1, 80]  # max linewidth acceptable
p.job["data_rejecting"]["ranges"]["chemical shift (ppm)"] = 0.5  # +/- ppm
#p.job["data_rejecting"]["ranges"]["phase std. factor (%)"] = 60.0  # stan's phase +/- 60% of std criteria (see doi:10.1002/jmri.26802)
# auto rejection based on linewidth?
p.job["data_rejecting"]["auto_method_list"] = [reco.data_rejection_method.AUTO_AMPLITUDE, reco.data_rejection_method.AUTO_LINEWIDTH, reco.data_rejection_method.AUTO_FREQUENCY, reco.data_rejection_method.AUTO_PHASE]
# minimum allowed SNR change (%) when adjusting the linewidth criteria
p.job["data_rejecting"]["auto_allowed_snr_change"] = 0.0

# frequency realignement
p.job["realigning"]["POI_range_ppm"] = [4, 5.2]
p.job["realigning"]["inter_corr_mode"] = True

# select number of excitations to average
# p.job["averaging"]["na"] = 8
# let's denoize a little (5Hz exponential apodization)
p.job["apodizing"]["damping_hz"] = 10
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.job["calibrating"]["POI_true_ppm"] = 4.7
p.job["calibrating"]["POI_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.job["calibrating"]["POI_true_ppm"] = 1.4
# p.job["calibrating"]["POI_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.job["analyzing-snr"]["s_range_ppm"] = [2.9, 3.1]  # signal ppm range
p.job["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.job["analyzing-lw"]["range_ppm"] = [1, 1.5]

# display ppm range
p.job["displaying"]["range_ppm"] = [0.5, 5]
# run the process pipeline
p.settings["datasets_indexes"] = [index_to_process]
p.run()
# save this to db file
p.save_datasets()

reco.remove_grids_from_all_figs()

# %% need to fix a phase problem? manual phasing (to change: 0.0)
p._data_list[0] = p._data_list[0] * np.exp(1j * 0.0)
p._data_list[0].display_spectrum_1d()
# save this to db file
p.save_datasets()

# %% print data rejection stuff + NOVAPOR water lw
if(p._data_list[0].data_rejection is None):
    print(">>> No data rejection jowjow!")
else:
    n_initial = p._data_list[0].data_rejection[0]['Pre-rejection']['na']
    n_final = p._data_list[0].data_rejection[-1]['Post-rejection']['na']
    n_rejected = n_initial - n_final
    rejection_rate = n_rejected / n_initial * 100.0
    print(">>> Data rejected = %d/%d (%.2f%%)" % (n_rejected, n_initial, rejection_rate))

water_lw = p._data_list[0].data_ref.analyze_linewidth_1d([4.5, 5])
print(">>> Water LW = %.2fHz" % water_lw)
