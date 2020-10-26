#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script to be called once some data was processed using the MRSData2 class and is is currently in memory.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.fit as fit
import mrs.sim as sim
import mrs.aliases as xxx
import mrs.db as db
import mrs.log as log
import numpy as np

import pdb

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = db.data_db("/home/tangir/crmbm/acq_db/20201023_brain.pkl")
# rdb = db.data_db("/home/tangir/crmbm/acq_db/20201021_brain.pkl")

# %% select datasets


def sel_func(dataset_entry):
    r = (True)
    return(r)


# add a _fit suffix to the db file
rdb_sel = rdb.select_datasets(sel_func)
rdb_sel.link_to_file(rdb.db_file[:-4] + "_fit" + ".pkl")

# %% fit strategies

fit_stategies = {}

# easy
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_Cr_CH3,
    xxx.m_Cho_CH3,
    xxx.m_Water])

fit_stategies["easy"] = metabolites_fit

# medium
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_mI,
    xxx.m_Water])

fit_stategies["medium"] = metabolites_fit

# difficult
metabolites_fit = np.sort([
    xxx.m_NAA_CH3,
    xxx.m_NAA_CH2,
    xxx.m_Cr_CH3,
    xxx.m_Cr_CH2,
    xxx.m_Cho_CH3,
    xxx.m_Cho_CH2,
    xxx.m_Gln,
    xxx.m_Glu,
    xxx.m_mI,
    xxx.m_Tau,
    xxx.m_Water])

fit_stategies["hard"] = metabolites_fit

water_concentration = 55000.0  # mmol/kg

# %% run the fits!

# browse though datasets
for this_patient in rdb_sel:
    for this_study in rdb_sel[this_patient]:
        for this_scan in rdb_sel[this_patient][this_study]:
            # get the data
            this_dataset = rdb_sel[this_patient][this_study][this_scan]["dataset"]
            this_data = this_dataset["raw"]["data"]
            if(this_data is None):
                # no raw data? ok get the dicom
                this_data = this_scan["dataset"]["dcm"]["data"]

            # display the data
            this_data.display_spectrum_1d()

            # perform area integration and ref fit

            # metabolite db
            metabolites_fit = list(fit_stategies.values())[0]
            meta_bs = sim.metabolite_basis_set()
            meta_bs.initialize()

            # sequence
            # take sequence from dataset, usually sLASER
            seq = this_data.sequence
            # seq.pulse_rfc_real_shape_enable = False
            # seq.pulse_exc_duration = 0.001
            # seq.pulse_rfc_duration = 0.001
            # seq.spoiler_duration = 0.001

            # make it simple: SE sequence
            seq = sim.mrs_seq_press(seq.te, seq.tr, seq.na, seq.ds, seq.nuclei, seq.npts, seq.fs, seq.f0)
            seq.initialize(meta_bs)

            # %% area integration

            # for non water-suppressed data
            prefittool = fit.prefit_tool(this_data.data_ref, seq)
            prefittool.area_integration_peak_ranges = [0.5]
            prefittool.area_integration_peaks = [xxx.m_Water]
            prefittool.initialize()
            params_ref_area, params_ref_area_pnorm = prefittool.run()

            # for water-suppressed data
            prefittool = fit.prefit_tool(this_data, seq)
            prefittool.area_integration_peak_ranges = [0.1, 0.1, 0.1]
            prefittool.area_integration_peaks = [xxx.m_NAA_CH3, xxx.m_Cr_CH3, xxx.m_Cho_CH3]
            prefittool.initialize()
            params_area, params_area_pnorm = prefittool.run()

            # %% fit non water-suppressed data

            fittool = fit.fit_tool(this_data.data_ref, seq)
            # min/max fitting bounds
            fittool.params_min = fittool.params_min.set_default_water_min()
            fittool.params_max = fittool.params_max.set_default_water_max()
            # water concentration parameter
            fittool.params_min[xxx.m_Water, xxx.p_cm] = 0
            fittool.params_max[xxx.m_Water, xxx.p_cm] = 1000000.0
            # water linewidth
            fittool.params_min[xxx.m_Water, xxx.p_dd] = 10.0
            fittool.params_max[xxx.m_Water, xxx.p_dd] = 500.0
            # water frequency shift
            fittool.params_min[xxx.m_Water, xxx.p_df] = -50.0
            fittool.params_max[xxx.m_Water, xxx.p_df] = +50.0
            # water phase
            fittool.params_min[xxx.m_Water, xxx.p_dp] = -0.1
            fittool.params_max[xxx.m_Water, xxx.p_dp] = +0.1

            # initial fit values
            fittool.params_init = fittool.params_min.copy()
            fittool.params_init[xxx.m_Water, xxx.p_cm] = 1.0
            fittool.params_init[xxx.m_Water, xxx.p_dd] = 100.0
            fittool.params_init[xxx.m_Water, xxx.p_df] = 0.0
            fittool.params_init[xxx.m_Water, xxx.p_dp] = 0.0

            # numerical optimization parameters
            fittool.optim_xtol = 1e-9
            fittool.optim_ftol = 1e-9
            fittool.optim_gtol = 1e-9
            fittool.display_range_ppm = [3, 7]

            # run the fit
            fittool.initialize()
            [params_ref_fit, optim_result_ref] = fittool.run()

            # use various fit strategies
            for this_fit_strategy in fit_stategies:

                # %% prepare simulation machine

                # metabolite db
                metabolites_fit = fit_stategies[this_fit_strategy]
                meta_bs = sim.metabolite_basis_set()
                meta_bs.initialize()

                # sequence
                # take sequence from dataset, usually sLASER
                seq = this_data.sequence
                # seq.pulse_rfc_real_shape_enable = False
                # seq.pulse_exc_duration = 0.001
                # seq.pulse_rfc_duration = 0.001
                # seq.spoiler_duration = 0.001

                # make it simple: SE sequence
                seq = sim.mrs_seq_press(seq.te, seq.tr, seq.na, seq.ds, seq.nuclei, seq.npts, seq.fs, seq.f0)
                seq.initialize(meta_bs)

                # %% fit water-suppressed data

                seq.initialize(meta_bs)

                fittool = fit.fit_tool(this_data, seq)

                # default fitting bounds from muscle template
                fittool.params_min = fittool.params_min.set_default_min()
                fittool.params_max = fittool.params_max.set_default_max()

                # ranges for concentration
                fittool.params_min[:, xxx.p_cm] = 0.09
                fittool.params_max[:, xxx.p_cm] = 300.0
                fittool.params_max[xxx.m_Water, xxx.p_cm] = 100000.0

                # linewidth bounds for metabolites
                # use water linewidth
                fittool.params_min[xxx.m_All_MBs, xxx.p_dd] = 5
                fittool.params_max[xxx.m_All_MBs, xxx.p_dd] = params_ref_fit[xxx.m_Water, xxx.p_dd] + 5

                # initial concentrations
                fittool.params_init[:, xxx.p_cm] = 0.0
                fittool.params_init[metabolites_fit, xxx.p_cm] = 0.1

                # frequency shifts for metabolites
                fittool.params_min[xxx.m_All_MBs, xxx.p_df] = -15.0
                fittool.params_max[xxx.m_All_MBs, xxx.p_df] = 15.0

                # phase bounds for metabolites and lipids
                fittool.params_min[:, xxx.p_dp] = -0.1
                fittool.params_max[:, xxx.p_dp] = +0.1

                # linklock
                fittool.params_init.linklock[:] = 1
                fittool.params_init.linklock[metabolites_fit[0], :] = [0, -200, 0, -100]
                for im in metabolites_fit[1:]:
                    fittool.params_init.linklock[im, :] = [0, 200, 0, 100]

                # leave water phase free
                fittool.params_init.linklock[xxx.m_Water, :] = [0, 0, 0, 0]

                # numerical optimization parameters
                fittool.display_range_ppm = [0, 6]
                fittool.display_frequency = 2

                # run the fit
                fittool.initialize()
                [params_fit_final, optim_result] = fittool.run()

                # %% save the fit pipeline to the db
                fit_results = {"params_ref_area": params_ref_area,
                                "params_ref_area_pnorm": params_ref_area_pnorm,
                                "params_area": params_area,
                                "params_area_pnorm": params_area_pnorm,
                                "params_ref_fit": params_ref_fit,
                                "params_fit_final": params_fit_final,
                                "optim_result": optim_result,
                                "fittool": fittool}

                rdb_sel.save_fit(this_scan, this_fit_strategy, fit_results)
