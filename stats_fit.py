#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A script used for group study analysis.

@author: Tangi Roussel
"""

from __future__ import division

import os
import pickle
import datetime
import numpy as np
import mrs.metabase as xxx
import mrs.sim as sim

from IPython import get_ipython
import warnings
warnings.filterwarnings("ignore", ".*GUI is implemented*")
get_ipython().magic('clear')

result_folder = './20191102_quantification.results'

# init stat vectors
big_list_params_fit_all_linked_T2s_abs_water = []
big_list_params_fit_all_linked_T2s_abs_cre = []
big_list_params_fit_all_linked_T2s_ratios = []
big_list_params_prefit_abs = []
big_list_metabolites_prefit = []
big_list_metabolites_fit = []
big_list_fname = []
big_list_birthyear = []
big_list_sex = []

# browse pickle files and extract data
for fname in os.listdir(result_folder):
    if fname.endswith(".pkl"):
        with open(result_folder + '/' + fname, 'rb') as f:
            [params_fit_all, params_fit_all_linked_T2s, params_fit_all_linked_T2s_abs_cre, params_fit_all_linked_T2s_abs_water, params_fit_all_linked_T2s_ratios, params_fit_singlets, params_prefit, params_prefit_abs, params_prefit_ref, params_ref_fit, params_ref_fit_T2s, metabolites_fit, metabolites_prefit, link_singlets_s_m_parts, fname, birthyear, sex] = pickle.load(f)

            # store
            big_list_params_fit_all_linked_T2s_abs_water.append(params_fit_all_linked_T2s_abs_water)
            big_list_params_fit_all_linked_T2s_abs_cre.append(params_fit_all_linked_T2s_abs_cre)
            big_list_params_fit_all_linked_T2s_ratios.append(params_fit_all_linked_T2s_ratios)
            big_list_params_prefit_abs.append(params_prefit_abs)
            big_list_metabolites_prefit.append(metabolites_prefit)
            big_list_metabolites_fit.append(metabolites_fit)
            big_list_fname.append(fname)
            big_list_birthyear.append(birthyear)
            big_list_sex.append(sex)


# convert lists to numpy
big_np_params_fit_all_linked_T2s_abs_water = np.array(big_list_params_fit_all_linked_T2s_abs_water)
big_np_params_fit_all_linked_T2s_abs_cre = np.array(big_list_params_fit_all_linked_T2s_abs_cre)
big_np_params_fit_all_linked_T2s_ratios = np.array(big_list_params_fit_all_linked_T2s_ratios)
big_np_params_prefit_abs = np.array(big_list_params_prefit_abs)
big_np_metabolites_prefit = np.array(big_list_metabolites_prefit)
big_np_metabolites_fit = np.array(big_list_metabolites_fit)
big_np_fname = np.array(big_list_fname)

index_data2keep = [1, 2, 3, 4, 6]
# remove 311_sl (wrong TE) and 307_ap (low SNR) to reduce CV...

# stats about sex/age
age_mean = np.mean(datetime.datetime.now().year - big_list_birthyear)
age_std = np.std(datetime.datetime.now().year - big_list_birthyear)
print("> subjects' age = %.2f +- %.2f years old" % (age_mean, age_std))

sex_male = np.sum(big_list_sex == 2)
sex_female = np.sum(big_list_sex == 1)
print("> subjects' sex: male = %d | female = %d" % (sex_male, sex_female))

# stats abs. rel. to water
metabolites_fit = np.unique(big_list_metabolites_fit)
big_np_params_fit_all_linked_T2s_abs_water_mean = big_np_params_fit_all_linked_T2s_abs_water[index_data2keep, :].mean(axis=0)
big_np_params_fit_all_linked_T2s_abs_water_std = big_np_params_fit_all_linked_T2s_abs_water[index_data2keep, :].std(axis=0)

# stats abs. rel. to cre
metabolites_fit = np.unique(big_list_metabolites_fit)
big_np_params_fit_all_linked_T2s_abs_cre_mean = big_np_params_fit_all_linked_T2s_abs_cre[index_data2keep, :].mean(axis=0)
big_np_params_fit_all_linked_T2s_abs_cre_std = big_np_params_fit_all_linked_T2s_abs_cre[index_data2keep, :].std(axis=0)

# stats ratios X/Cre
metabolites_fit = np.unique(big_list_metabolites_fit)
big_np_params_fit_all_linked_T2s_ratios_mean = big_np_params_fit_all_linked_T2s_ratios[index_data2keep, :].mean(axis=0)
big_np_params_fit_all_linked_T2s_ratios_std = big_np_params_fit_all_linked_T2s_ratios[index_data2keep, :].std(axis=0)

# need a simulation machine
simtool = sim.metabolite_database()
simtool.initialize()

# convert numpy to params object
big_params_params_fit_all_linked_T2s_abs_water_mean = sim.params(simtool)
big_params_params_fit_all_linked_T2s_abs_water_std = sim.params(simtool)
#
big_params_params_fit_all_linked_T2s_abs_cre_mean = sim.params(simtool)
big_params_params_fit_all_linked_T2s_abs_cre_std = sim.params(simtool)
#
big_params_params_fit_all_linked_T2s_ratios_mean = sim.params(simtool)
big_params_params_fit_all_linked_T2s_ratios_std = sim.params(simtool)

big_params_params_fit_all_linked_T2s_abs_water_mean[:] = big_np_params_fit_all_linked_T2s_abs_water_mean
big_params_params_fit_all_linked_T2s_abs_water_std[:] = big_np_params_fit_all_linked_T2s_abs_water_std
#
big_params_params_fit_all_linked_T2s_abs_cre_mean[:] = big_np_params_fit_all_linked_T2s_abs_cre_mean
big_params_params_fit_all_linked_T2s_abs_cre_std[:] = big_np_params_fit_all_linked_T2s_abs_cre_std
#
big_params_params_fit_all_linked_T2s_ratios_mean[:] = big_np_params_fit_all_linked_T2s_ratios_mean
big_params_params_fit_all_linked_T2s_ratios_std[:] = big_np_params_fit_all_linked_T2s_ratios_std

# human range
p_human_min = sim.params(simtool).set_default_human_min(False)
p_human_max = sim.params(simtool).set_default_human_max(False)
p_human_ratios_min = p_human_min._correct_normalize_concentrations(1.0 / p_human_min[xxx.m_Cre_303ppm, xxx.p_cm])
p_human_ratios_max = p_human_max._correct_normalize_concentrations(1.0 / p_human_max[xxx.m_Cre_303ppm, xxx.p_cm])

# disp
sim.disp_bargraph(1000, [p_human_min,
                         big_params_params_fit_all_linked_T2s_abs_water_mean,
                         big_params_params_fit_all_linked_T2s_abs_cre_mean,
                         p_human_max],
                  [p_human_min * 0.0,
                   big_params_params_fit_all_linked_T2s_abs_water_std,
                   big_params_params_fit_all_linked_T2s_abs_cre_std,
                   p_human_max * 0.0],
                  ['Human MIN',
                   'Fit rel. to Water',
                   'Fit rel. to Cre',
                   'Human MAX'], False, False, False, False, xxx.p_cm, metabolites_fit, 0.2)

big_params_params_fit_all_linked_T2s_ratios_mean.print()
big_params_params_fit_all_linked_T2s_ratios_std.print()

im = xxx.m_Cho_NoTriMethyl
mm = big_params_params_fit_all_linked_T2s_ratios_mean[im, 0]
stdd = big_params_params_fit_all_linked_T2s_ratios_std[im, 0]
cvv = stdd / mm * 100.0
print("> tCho/tCr = %.2f (%.2f) CV=%.2f" % (mm, stdd, cvv))

im = xxx.m_NAA_s
mm = big_params_params_fit_all_linked_T2s_ratios_mean[im, 0]
stdd = big_params_params_fit_all_linked_T2s_ratios_std[im, 0]
cvv = stdd / mm * 100.0
print("> tNAA/Cr = %.2f (%.2f) CV=%.2f" % (mm, stdd, cvv))

im = xxx.m_Myo
mm = big_params_params_fit_all_linked_T2s_ratios_mean[im, 0]
stdd = big_params_params_fit_all_linked_T2s_ratios_std[im, 0]
cvv = stdd / mm * 100.0
print("> mI/Cr = %.2f (%.2f) CV=%.2f" % (mm, stdd, cvv))
