#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vivo data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.aliases as xxx
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

# %% process data
get_ipython().magic("clear")
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
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID203_slaser_20_1_LONGTE_FID50322.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID247_steam_shortTE_1024_NOOVS_FID48377.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/coeur/meas_MID141_slaser2_WSAT_ECG_5NEX_512_FID48751.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID131_eja_svs_slaser_VAPOR_FID23607.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID150_eja_svs_slaser_VAPOR_FID23626.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID118_SVSSE_16EX_OVS_WSAT_FID23594.dat
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
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID55_eja_svs_slaser_VAPOR_BH3_FID25884.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID49_eja_svs_slaser_VAPOR_BH2_FID25878.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID38_eja_svs_steam_FID25867.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID51_eja_svs_steam_VAPOR_BH2_FID25880.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID44_eja_svs_press_FID25873.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID53_eja_svs_press_VAPOR_BH2_FID25882.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID46_eja_svs_slaser_VAPOR_FB_FID25875.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID47_eja_svs_steam_VAPO_FB_FID25876.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID48_eja_svs_press_VAPOR_FB_FID25877.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID108_slaser_20_1_TE200_FID52352.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID128_slaser_20_1_LONGTE_FID50247.dat
T:/users/JS/2020_coeur_mrs7T/200311_JB_358/meas_MID62_slaser_WSAT_PULS_100NEX_FID54390.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID146_eja_svs_slaser_VAPOR_BH_FID27306.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID206_eja_svs_slaser_VAPOR_BH_FID27366.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID199_eja_svs_steam_VAPOR_BH_FID27359.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID208_eja_svs_steam_VAPOR_BH_FID27368.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID200_eja_svs_press_VAPOR_BH_FID27360.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID209_eja_svs_press_VAPOR_BH_FID27369.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID201_eja_svs_slaser_VAPOR_FB_FID27361.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID202_eja_svs_steam_VAPOR_FB_FID27362.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID203_eja_svs_press_VAPOR_FB_FID27363.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID85_eja_svs_slaser_VAPOR_BH2_FID27450.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID77_eja_svs_slaser_VAPOR_BH2_FID27442.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID71_eja_svs_steam_VAPOR_BH_FID27436.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID91_eja_svs_steam_VAPOR_BH_FID27456.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID88_eja_svs_press_VAPOR_BH_FID27453.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID74_eja_svs_press_VAPOR_BH_FID27439.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID80_eja_svs_slaser_VAPOR_FB_FID27445.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID84_eja_svs_steam_VAPOR_FB_FID27449.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID82_eja_svs_press_VAPOR_FB_FID27447.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID257_eja_svs_slaser_VAPOR_BH_FID27622.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID246_eja_svs_slaser_VAPOR_BH_FID27611.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID250_eja_svs_steam_VAPOR_BH_FID27615.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID223_eja_svs_steam_VAPOR_BH_FID27588.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID226_eja_svs_press_VAPOR_BH_FID27591.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID253_eja_svs_press_VAPOR_BH_FID27618.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID259_eja_svs_slaser_VAPOR_FB_FID27624.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID233_eja_svs_steam_VAPOR_FB_FID27598.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID235_eja_svs_press_VAPOR_FB_FID27600.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID287_eja_svs_slaser_VAPOR_BH_FID27652.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID311_eja_svs_slaser_VAPOR_BH_FID27676.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID291_eja_svs_steam_VAPOR_BH_FID27656.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID290_eja_svs_steam_VAPOR_BH_FID27655.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID312_eja_svs_press_VAPOR_BH_FID27677.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID294_eja_svs_press_VAPOR_BH_FID27659.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID297_eja_svs_slaser_VAPOR_FB_FID27662.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID298_eja_svs_steam_VAPOR_FB_FID27663.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID301_eja_svs_press_VAPOR_FB_FID27666.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID57_eja_svs_steam_VAPOR_BH_FID27422.dat
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
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t346/muscle/meas_MID203_slaser_20_1_LONGTE_FID50322.dat
T:/users/SSSR/Data/2020_NUTRI7T/191209_nutri7T/coeur/meas_MID255_steam_1024_NOOVS_NOWSAT_FID48385.dat
T:/users/SSSR/Data/2020_NUTRI7T/191212_nutri7T/coeur/meas_MID142_slaser2_NOWSAT_ECG_5NEX_512_FID48752.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001/meas_MID134_eja_svs_slaser_NOVAPOR_FID23610.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001//meas_MID134_eja_svs_slaser_NOVAPOR_FID23610.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200206_1HMRS_zb001//meas_MID121_SVSSE_16EX_OVS_NOWSAT_FID23597.dat
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
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID50_eja_svs_slaser_NOVAPOR_BH2_FID25879.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID50_eja_svs_slaser_NOVAPOR_BH2_FID25879.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID52_eja_svs_steam_NOVAPOR_BH2_FID25881.dat
/home/tangir/desktop/200303_1HMRS_pd002:meas_MID52_eja_svs_steam_NOVAPOR_BH2_FID25881.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID54_eja_svs_press_NOVAPOR_BH2_FID25883.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID54_eja_svs_press_NOVAPOR_BH2_FID25883.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID50_eja_svs_slaser_NOVAPOR_BH2_FID25879.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID52_eja_svs_steam_NOVAPOR_BH2_FID25881.dat
/home/tangir/desktop/200303_1HMRS_pd002/meas_MID54_eja_svs_press_NOVAPOR_BH2_FID25883.dat
T:/users/SSSR/Data/2020_NUTRI7T/200213_nutri7t352/meas_MID108_slaser_20_1_TE200_FID52352.dat
T:/users/SSSR/Data/2020_NUTRI7T/200122_nutri7t345/muscle/meas_MID128_slaser_20_1_LONGTE_FID50247.dat
T:/users/JS/2020_coeur_mrs7T/200311_JB_358/meas_MID63_slaser_WATER_PULS_10NEX_FID54391.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID147_eja_svs_slaser_NOVAPOR_BH_FID27307.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID147_eja_svs_slaser_NOVAPOR_BH_FID27307.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID149_eja_svs_steam_NOVAPOR_BH_FID27309.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID149_eja_svs_steam_NOVAPOR_BH_FID27309.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID153_eja_svs_press_NOVAPOR_BH_FID27313.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID153_eja_svs_press_NOVAPOR_BH_FID27313.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID147_eja_svs_slaser_NOVAPOR_BH_FID27307.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID149_eja_svs_steam_NOVAPOR_BH_FID27309.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200312_1HMRS_ma003/meas_MID153_eja_svs_press_NOVAPOR_BH_FID27313.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID86_eja_svs_slaser_NOVAPOR_BH2_FID27451.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID78_eja_svs_slaser_NOVAPOR_BH2_FID27443.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID72_eja_svs_steam_NOVAPOR_BH_FID27437.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID92_eja_svs_steam_NOVAPOR_BH_FID27457.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID89_eja_svs_press_NOVAPOR_BH_FID27454.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID75_eja_svs_press_NOVAPOR_BH_FID27440.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID86_eja_svs_slaser_NOVAPOR_BH2_FID27451.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID92_eja_svs_steam_NOVAPOR_BH_FID27457.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID75_eja_svs_press_NOVAPOR_BH_FID27440.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID238_eja_svs_slaser_NOVAPOR_BH_FID27603.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID220_eja_svs_slaser_NOVAPOR_BH_FID27585.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID251_eja_svs_steam_NOVAPOR_BH_FID27616.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID224_eja_svs_steam_NOVAPOR_BH_FID27589.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID227_eja_svs_press_NOVAPOR_BH_FID27592.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID254_eja_svs_press_NOVAPOR_BH_FID27619.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID248_eja_svs_slaser_NOVAPOR_BH_FID27613.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID251_eja_svs_steam_NOVAPOR_BH_FID27616.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/200313_1HMRS_tr005/meas_MID254_eja_svs_press_NOVAPOR_BH_FID27619.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID288_eja_svs_slaser_NOVAPOR_BH_FID27653.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID288_eja_svs_slaser_NOVAPOR_BH_FID27653.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID292_eja_svs_steam_NOVAPOR_BH_FID27657.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID292_eja_svs_steam_NOVAPOR_BH_FID27657.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID295_eja_svs_press_NOVAPOR_BH_FID27660.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID295_eja_svs_press_NOVAPOR_BH_FID27660.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID288_eja_svs_slaser_NOVAPOR_BH_FID27653.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID292_eja_svs_steam_NOVAPOR_BH_FID27657.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_sr006/meas_MID295_eja_svs_press_NOVAPOR_BH_FID27660.dat
C:/Users/jsourdon/Desktop/2020_1H_MRS/201303_1hmrs_af004/meas_MID92_eja_svs_steam_NOVAPOR_BH_FID27457.dat
"""

# legend captions
display_legends = """
0_191122_331
1_191122_331
2_191125_332_slaser
3_191125_332
4_191129_334_steam
5_191129_334_slaser
6_191203_DICOM
7_191203_DICOM
8_191209_337_slaser
9_191209_337_steam
10_191212_339_slaser
11_191212_339_steam
12_200122_345_slaser
13_200122_345_steam
14_200122_346_slaser
15_200122_346_steam
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
37_1HMRS_002_slaser1
38_1HMRS_002_slaser2
39_1HMRS_002_steam1
40_1HMRS_002_steam2
41_1HMRS_002_press1
42_1HMRS_002_press2
43_1HMRS_002_slaserFR
44_1HMRS_002_steamFR
45_1HMRS_002_pressFR
46_29_352_slaser_longTE
47_345_slaser_TE350
48_358_slaser_coeur_FR
49_1HMRS_003_slaser1
50_1HMRS_003_slaser2
51_1HMRS_003_steam1
52_1HMRS_003_steam2
53_1HMRS_003_press1
54_1HMRS_003_press2
55_1HMRS_003_slaserFR
56_1HMRS_003_steamFR
57_1HMRS_003_pressFR
58_1HMRS_004_slaser1
59_1HMRS_004_slaser2
60_1HMRS_004_steam1
61_1HMRS_004_steam2
62_1HMRS_004_press1
63_1HMRS_004_press2
64_1HMRS_004_slaserFR
65_1HMRS_004_steamFR
66_1HMRS_004_pressFR
67_1HMRS_005_slaser1
68_1HMRS_005_slaser2
69_1HMRS_005_steam1
70_1HMRS_005_steam2
71_1HMRS_005_press1
72_1HMRS_005_press2
73_1HMRS_005_slaserFR
74_1HMRS_005_steamFR
75_1HMRS_005_pressFR
76_1HMRS_006_slaser1
77_1HMRS_006_slaser2
78_1HMRS_006_steam1
79_1HMRS_006_steam2
80_1HMRS_006_press1
81_1HMRS_006_press2
82_1HMRS_006_slaserFR
83_1HMRS_006_steamFR
84_1HMRS_006_pressFR
85_1HMRS_004_steam3
"""

# the data lists above, which data index should we process (starting from 0)
# [] means process all
index_to_process = 45

# --- process the water-suppressed (WS) data ---
p = reco.pipeline()
p.data_filepaths = data_filepaths
p.data_ref_filepaths = data_ref_filepaths
p.data_physio_filepaths = []
p.display_legends = display_legends
# we are using a single channel RF coil
p.data_coil_nChannels = 1

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
                p.jobs["displaying"]
                ]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]]

# self-phasing the WS spectrum using the the big lipid peak at ~1.5ppm (Here you can choose btw Lipid and Cr for phasing)
p.jobs["phasing"]["POI_range_ppm"] = [4, 5.2]  # find this peak in this ppm range
p.jobs["phasing"]["offset"] = 0.0  # manual phase offset
p.jobs["phasing"]["using_ref_data"] = False

# reject bad data
p.jobs["data-rejecting"]["moving_averages"] = 1
p.jobs["data-rejecting"]["POI_range_ppm"] = [4, 5.2]
# rejection ranges
p.jobs["data-rejecting"]["ranges"]["amplitude (%)"] = None  # do not reject on amplitude changes
p.jobs["data-rejecting"]["ranges"]["linewidth (Hz)"] = 20  # max linewidth acceptable
p.jobs["data-rejecting"]["ranges"]["chemical shift (ppm)"] = 0.5  # +/- ppm
p.jobs["data-rejecting"]["ranges"]["phase std. factor (%)"] = 60.0  # stan's phase +/- 60% of std criteria (see doi:10.1002/jmri.26802)
# auto rejection based on linewidth?
p.jobs["data-rejecting"]["auto"] = True
# minimum allowed SNR change (%) when adjusting the linewidth criteria
p.jobs["data-rejecting"]["auto_allowed_snr_change"] = -10.0

# frequency realignement
p.jobs["realigning"]["POI_range_ppm"] = [4, 5.2]
# select number of excitations to average
# p.jobs["averaging"]["na"] = 8
# let's denoize a little (5Hz exponential apodization)
p.jobs["apodizing"]["damping_hz"] = 10
# let's calibrate ppm scale so that water residue is at 4.7ppm
p.jobs["calibrating"]["POI_true_ppm"] = 4.7
p.jobs["calibrating"]["POI_range_ppm"] = [4, 5.2]
# or if the residue is too small, let's use the lipid at 1.5ppm
# p.jobs["calibrating"]["POI_true_ppm"] = 1.4
# p.jobs["calibrating"]["POI_range_ppm"] = [1, 2]

# snr and linewidth estimation on Cr peak
p.jobs["analyzing-snr"]["s_range_ppm"] = [2.5, 3.5]  # signal ppm range
p.jobs["analyzing-snr"]["n_range_ppm"] = [-3, -1]  # noise ppm range
p.jobs["analyzing-lw"]["range_ppm"] = [2.5, 3.5]

# display ppm range
p.jobs["displaying"]["range_ppm"] = [0, 5]
# run the process pipeline
p.data_process_only_this_data_index = [index_to_process]
p.run()
# save this to db file
p.save(rdb)

# %% need to fix a phase problem? manual phasing
p._data_list[0] = p._data_list[0] * np.exp(1j * 0.0)
p._data_list[0].display_spectrum_1d()
# save this to db file
p.save(rdb)

# %% print data rejection stuff
n_initial = p._data_list[0].na_pre_data_rejection
n_final = p._data_list[0].na_post_data_rejection
n_rejected = n_initial - n_final
rejection_rate = n_rejected / n_initial * 100.0

print(">>> Data rejected = %d/%d (%.2f%%)" % (n_rejected, n_initial, rejection_rate))

# %% replace
"C:/Users/jsourdon/Desktop/2020_1H_MRS/200303_1HMRS_pd002/"
# by
"/home/tangir/desktop/200303_1HMRS_pd002"
