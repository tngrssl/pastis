#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A fit_tool class that deals with the quantification of MRS data based on the scipy's least_squares optimisation function.

@author: Tangi Roussel
"""

import numpy as np
from mrs import sim
from mrs import log
from mrs import aliases as xxx
import matplotlib.pylab as plt
import scipy.optimize as optimize
import matplotlib._color_data as mcd
from enum import Enum
import time

import pdb

PARS_LONG_NAMES = ["Concentration [mmol/kg]", "Linewidth factor [Hz]", "Frequency shift [Hz]", "Phase [rd]"]
PARS_SHORT_NAMES = ["cm", "dd", "df", "dp"]


class fit_plot_type(Enum):
    """The enum fit_plot_type describes the type of plots that can be displayed real-time during the fit."""

    BARGRAPH_CM = 1
    BARGRAPH_DD = 2
    BARGRAPH_DF = 3
    BARGRAPH_DP = 4
    COST_FUNCTION = 5
    CORR_MATRIX = 6
    TIME_DOMAIN = 7


import os
import pickle
from datetime import datetime, timedelta
from shutil import copyfile


def _dummy_check_func(d, rp, fp, dp):
    """
    Return true if you want to select this dataset/pipeline when getting data from database using the get_datasets method. This is only a dummy example. Please feel free to recode a function with the same prototype to select for example only 3T scans, or short-TE scans, etc.

    Parameters
    ----------
    d : MRSData2 object
        Function with two arguments (MRSData2 object, pipeline object) that you should code and which returns True if you want to get it this dataset/pipeline out of the database.
    rp : reco pipeline object
    fp : fit pipeline object
    dp : dict of parameter results

    Returns
    -------
    r : boolean
        True if we keep this dataset/pipeline
    """
    r = (d.te > 0.0 and
         d.na > 0)  # that is stupid, just an example
    return(r)


class prefit_tool:
    """The prefit_tool class is used to perfom area integration of peaks in the spectrum."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, data, seq=None):
        """
        Construct a prefit_tool object.

        Parameters
        ----------
        data : MRSData2 numpy array [timepoints]
            MRS data to fit stored in a MRSData2 object
        seq : mrs_sequence object
            Sequence used for the acquisition of the data. By default, the sequence used to acquire the data.
        """
        # data signal to process
        self.data = data

        # sequence used to acquire and therefore to modelize the signal
        if(seq is None):
            self._seq = self.data.sequence
        else:
            self._seq = seq

        # initialize the sequence now if needed
        if(not self.seq.ready):
            self.seq.initialize()

        # metabolite db
        self._meta_bs = self.seq.meta_bs

        # metabolites to integrate
        self.area_integration_peaks = [xxx.m_Water]
        # ppm range to look for peak maximum
        self.area_integration_peak_ranges = [0.4]  # ppm
        # internal parameters
        self._peak_names = []
        self._peak_ppms = []
        self._peak_nprots = []
        self._peak_base_width = []
        self._peak_areas = []
        self._peak_areas_norm = []

        # display options
        self.display_enable = True
        self.display_fig_index = 1000
        self.display_range_ppm = [0, 6]  # ppm

        # to know if the object is initialized
        self._ready = False

        # freeze
        self.__isfrozen = True

    @property
    def ready(self):
        """
        Property get function for _ready.

        Returns
        -------
        self._ready : bool
            to tell if the object if initialized or not
        """
        return(self._ready)

    @property
    def seq(self):
        """
        Property get function for _seq.

        Returns
        -------
        self._seq : mrs_sequence object
            Sequence used for the acquisition of the data
        """
        return(self._seq)

    @property
    def meta_bs(self):
        """
        Property get function for meta_bs.

        Returns
        -------
        self._meta_bs : metabolite_basis_set() object
            Metabolite database to use for simulation
        """
        return(self._meta_bs)

    def initialize(self):
        """Initialize the prefit procedure."""
        log.info("initializing...")

        # clear all
        self._peak_names = []
        self._peak_ppms = []
        self._peak_nprots = []

        # first find chemical shifts for those peaks
        meta_keys = list(self.meta_bs.keys())
        integration_metagroup_keys = [meta_keys[ind] for ind in self.area_integration_peaks]
        for this_metagroup_key in integration_metagroup_keys:
            meta_found_singulet = None
            for this_meta in list(self.meta_bs[this_metagroup_key]["metabolites"].items()):
                # find first meta in this metagroup
                this_meta_name = this_meta[0]
                this_meta_ppm = this_meta[1]["ppm"][0]
                this_meta_nprots = float(len(this_meta[1]["ppm"]))
                this_meta_j = this_meta[1]["J"]

                # check that it is a singlet: only one meta, all Js are zeroes and only one unique ppm
                if(np.all(this_meta_j == 0.0) and len(set(this_meta[1]["ppm"])) == 1):
                    meta_found_singulet = this_meta

            if(meta_found_singulet is None):
                log.error("> one of the metabolites you chose for peak integration does not have a singlet, aborting area integration...")
            else:
                # store
                self._peak_names.append(this_meta_name)
                self._peak_ppms.append(this_meta_ppm)
                self._peak_nprots.append(this_meta_nprots)

        # I am ready now
        self._ready = True

    def run(self):
        """
        Run the area integration pipeline.

        Returns
        -------
        pars : params object
            Estimated relative metabolic concentration using area integration
        pars_norm : params object
            Estimated relative metabolic concentration using area integration, normalized by proton number
        """
        log.info_line________________________()
        log.info("running peak integration...")
        log.info_line________________________()

        # ready or not, here I come
        if(not self.ready):
            log.error("this prefit_tool object was not initialized!")

        # init
        s = self.data.copy()
        ppm = s.frequency_axis_ppm()

        # find maximum peak in range and its chemical shift
        sf = s.spectrum()
        sf_real = np.real(sf)

        if(self.display_enable):
            fig = plt.figure(self.display_fig_index)
            fig.clf()
            axs = fig.subplots()
            fig.canvas.set_window_title("mrs.fit.prefit_tool.run")
            axs.plot(self.data.frequency_axis_ppm(), self.data.spectrum().real, "k-", label="data")

        # now for each peak, integrate the area
        self._peak_areas = []
        self._peak_areas_norm = []
        pars = sim.params(self.meta_bs)
        pars[:] = 0.0
        pars_pnorm = sim.params(self.meta_bs)
        pars_pnorm[:] = 0.0

        log.info("integrating peak for...")
        log.info_line________________________()
        for (this_peak_index, this_peak_range, this_peak_name, this_peak_ppm, this_peak_np) in zip(self.area_integration_peaks, self.area_integration_peak_ranges, self._peak_names, self._peak_ppms, self._peak_nprots):
            log.info("[%s] theoretically at %.2fppm in theory" % (this_peak_name, this_peak_ppm))

            # first find closest peak in range
            this_peak_search_range = [this_peak_ppm - this_peak_range / 2.0, this_peak_ppm + this_peak_range / 2.0]
            _, this_peak_ppm, _, this_peak_lw_hz, _, _ = self.data._analyze_peak_1d(this_peak_search_range)
            log.debug("found it in the spectrum at %.2fppm" % this_peak_ppm)

            # integrate area

            # according to https://doi.org/10.1002/nbm.4257
            # peak integration area should be 2 * peak linewidth
            this_peak_lw_ppm = this_peak_lw_hz / self.data.f0
            if(this_peak_lw_ppm > this_peak_range):
                this_peak_lw_ppm = this_peak_range

            ippm_peak_range = (ppm > (this_peak_ppm - this_peak_lw_ppm)) & (ppm < (this_peak_ppm + this_peak_lw_ppm))
            sf_to_integrate = sf_real[ippm_peak_range]
            ppm_to_integrate = ppm[ippm_peak_range]
            sf_to_integrate = sf_to_integrate[::-1]
            ppm_to_integrate = ppm_to_integrate[::-1]
            this_peak_area = np.trapz(sf_to_integrate, ppm_to_integrate)
            log.info("got a peak area of %f" % this_peak_area)
            log.info("and that is %f when normalized to number of 1H (%.0f)" % (this_peak_area / this_peak_np, this_peak_np))
            log.info_line________________________()

            # store
            self._peak_areas.append(this_peak_area)
            self._peak_areas_norm.append(this_peak_area / this_peak_np)
            pars[this_peak_index, 0] = this_peak_area
            pars_pnorm[this_peak_index, 0] = this_peak_area / this_peak_np

            # display the area
            if(self.display_enable):
                axs.fill_between(ppm_to_integrate, 0, sf_to_integrate, label=this_peak_name)

        # finalize figure
        if(self.display_enable):
            axs.set_xlim(self.display_range_ppm[1], self.display_range_ppm[0])
            axs.set_xlabel('chemical shift (ppm)')
            axs.set_ylabel('spectrum')
            axs.grid('on')
            axs.legend()
            plt.pause(0.5)

        # return a param vectors
        return(pars, pars_pnorm)


class fit_tool:
    """The fit_tool class is used to adjust the MRS model implemented in the mrs.sim.metabolite_database class on some experimentally acquired data."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, data, seq=None):
        """
        Construct a prefit_tool object.

        Parameters
        ----------
        data : MRSData2 numpy array [timepoints]
            MRS data to fit stored in a MRSData2 object
        seq : mrs_sequence object
            Sequence used for the acquisition of the data. By default, the sequence used to acquire the data.
        """
        # data signal to process
        self.data = data

        # sequence used to acquire and therefore to modelize the signal
        if(seq is None):
            self._seq = self.data.sequence
        else:
            self._seq = seq

        # initialize the sequence now if needed
        if(not self.seq.ready):
            self.seq.initialize()

        # metabolite db
        self._meta_bs = self.seq.meta_bs

        # --- defaut lower bound parameters ---
        self.params_min = sim.params(self.meta_bs)
        self.params_min.set_default_min()
        self.params_min[xxx.m_All_MBs, xxx.p_dd] = 5.0
        self.params_min[xxx.m_All_MBs, xxx.p_df] = -10.0
        self.params_min[xxx.m_All_MBs, xxx.p_dp] = -np.pi / 4.0

        # --- default upper boud parameters ---
        self.params_max = sim.params(self.meta_bs)
        self.params_max.set_default_max()
        self.params_max[xxx.m_All_MBs, xxx.p_dd] = 50.0
        self.params_max[xxx.m_All_MBs, xxx.p_df] = +10.0
        self.params_max[xxx.m_All_MBs, xxx.p_dp] = +np.pi / 4.0

        # --- default initial parameters ---
        self.params_init = (self.params_min + self.params_max) / 2.0
        # with early minimal concentrations
        self.params_init[xxx.m_All_MBs, xxx.p_cm] = self.params_min[xxx.m_All_MBs, xxx.p_cm] * 1.1
        # with early minimal dmaping
        self.params_init[xxx.m_All_MBs, xxx.p_dd] = self.params_min[xxx.m_All_MBs, xxx.p_dd] * 1.1

        # link-lock vectors should be the same for initial, minimum and maximum parameter sets: link them here
        self.params_min._linklock = self.params_init.linklock
        self.params_max._linklock = self.params_init.linklock

        # private option to know if we are dealing with a water fit
        self._water_only = False

        # --- optimization options ---
        # should we use jacobin information during the fit or not?
        self.optim_jacobian = True
        # ppm range (By default, use metabolite_basis_set ppm range)
        self.optim_ppm_range = self.seq.meta_bs.ppm_range
        # stop fit if parameter changes go below this tolerance
        self.optim_xtol = 1e-9
        # stop fit if (error?) function changes go below this tolerance
        self.optim_ftol = 1e-9
        # stop fit if gradient (?) changes go below this tolerance
        self.optim_gtol = 1e-9
        # count number of calls to error function
        self._model_call_count = 0
        # record cost function
        self._cost_function = []
        # time
        self._fit_time = None
        # FQN noise region
        self.fqn_noise_range = [-2, -1]

        # --- display options ---
        self.display_enable = True
        # figure index
        self.display_fig_index = 2000
        # ppm range
        self.display_range_ppm = [1, 6]  # ppm
        # display every n calls to error function
        self.display_frequency = 20
        # describes the type of plots showed in the 4 subplots during the fit
        self.display_subplots_types = [fit_plot_type.BARGRAPH_CM, fit_plot_type.BARGRAPH_DD, fit_plot_type.BARGRAPH_DF, fit_plot_type.CORR_MATRIX]

        # --- display options ---
        # should we show the CRBs error bars during the fit?
        self.display_CRBs = True

        # to know if the object is initialized
        self._ready = False

        # freeze
        self.__isfrozen = True

    @property
    def ready(self):
        """
        Property get function for _ready.

        Returns
        -------
        self._ready : bool
            to tell if the object if initialized or not
        """
        return(self._ready)

    @property
    def seq(self):
        """
        Property get function for _seq.

        Returns
        -------
        self._seq : mrs_sequence object
            Sequence used for the acquisition of the data
        """
        return(self._seq)

    @property
    def meta_bs(self):
        """
        Property get function for meta_bs.

        Returns
        -------
        self._meta_bs : metabolite_basis_set() object
            Metabolite database to use for simulation
        """
        return(self._meta_bs)

    def _jac(self, params_free):
        """
        Return the jacobian vector needed for optimization.

        Parameters
        ----------
        params_free : 1D numpy array of free parameters
            An array of fit parameter values which are free according to the link-lock (LL) array.

        Returns
        -------
        j_free : numpy array [number of free parameters,number of time points]
            Jacobian model used to numerical optimization
        """
        # link-lock stuff
        params = self.params_init.copy()
        params = params.toFullParams(params_free)

        # calling jacobian model
        j_full = self.seq._jac(params)

        # reducing it to free params
        j_free = j_full[params.linklock <= 0, :]
        j_free = np.real(j_free)

        # fix for complex numbers
        j_free_cmplx_odd_even = np.zeros([j_free.shape[0], j_free.shape[1] * 2], dtype=np.float64)
        j_free_cmplx_odd_even[:, 0::2] = j_free.real
        j_free_cmplx_odd_even[:, 1::2] = j_free.imag

        # transpose
        j_free_cmplx_odd_even_tp = np.transpose(j_free_cmplx_odd_even)
        return(j_free_cmplx_odd_even_tp)

    def _minimizeThis(self, params_free, fit_terminated=False):
        """
        Return the difference between the model and the data, also called residue. The user should not call this function. it is the job of the optimization function that will try to mininize the residue using least-squares and gradient descent stuff, etc.

        Parameters
        ----------
        params_free : 1D numpy array of free parameters
            An array of fit parameter values which are free according to the link-lock (LL) array
        fit_terminated : bool
            If fit terminated, change figure title :)

        Returns
        -------
        diff_cmplx_odd_even : numpy array
            Residue with real and imaginary values interleaved (only way here to deal with numerical optimization with complex numbers!)
        """
        # link-lock stuff
        params = self.params_init.copy()
        params = params.toFullParams(params_free)

        # and boum the model
        mod = self.seq._model(params)
        # and the difference with the data
        diff = mod - self.data

        # fit criteria
        rsq_t, rsq_f = self.get_Rsq(params)
        fqn = self.get_FQN(params)

        # cost
        current_err = np.sum(np.abs(np.power(diff, 2)))
        self._cost_function.append(current_err)

        # iterate counter
        self._model_call_count = self._model_call_count + 1

        # display
        if(self.display_enable and (fit_terminated or np.mod(self._model_call_count, self.display_frequency) == 0)):
            fig = plt.figure(self.display_fig_index)
            fig.clf()
            axs = fig.subplots(2, 4)
            fig.canvas.set_window_title("mrs.fit.lsqfit._minimizeThis")

            # display real-time bargraphs
            # real-time CRBs
            if(self.display_CRBs):
                params_CRBs_abs, _ = self.get_CRBs(params)
                params._errors = params_CRBs_abs

            params_val_list = [self.params_min, self.params_max, params]
            params_leg_list = ['lower bounds', 'upper bounds', 'current']

            # let's display the subplots as required by self.display_subplots_types
            axs_list = [axs[0, 2], axs[0, 3], axs[1, 2], axs[1, 3]]
            for ax, plt_type in zip(axs_list, self.display_subplots_types):
                if(plt_type == fit_plot_type.BARGRAPH_CM):
                    disp_bargraph(ax, params_val_list, params_leg_list, True, True, self._water_only, True, xxx.p_cm)
                elif(plt_type == fit_plot_type.BARGRAPH_DD):
                    disp_bargraph(ax, params_val_list, params_leg_list, True, True, self._water_only, True, xxx.p_dd)
                elif(plt_type == fit_plot_type.BARGRAPH_DF):
                    disp_bargraph(ax, params_val_list, params_leg_list, True, True, self._water_only, True, xxx.p_df)
                elif(plt_type == fit_plot_type.BARGRAPH_DP):
                    disp_bargraph(ax, params_val_list, params_leg_list, True, True, self._water_only, True, xxx.p_dp)
                elif(plt_type == fit_plot_type.COST_FUNCTION):
                    ax.plot(self._cost_function, 'k-')
                    ax.set_xlabel('iterations')
                    ax.set_ylabel('LSQ fit residue')
                    ax.grid('on')
                elif(plt_type == fit_plot_type.CORR_MATRIX):
                    self.disp_corr_mat(ax, params)
                elif(plt_type == fit_plot_type.TIME_DOMAIN):
                    ax.plot(self.data.time_axis(), np.real(self.data), 'k-', linewidth=0.5)
                    ax.plot(mod.time_axis(), np.real(mod), 'k-', linewidth=2)
                    ax.plot(diff.time_axis(), np.real(diff), 'g-', linewidth=0.5)
                    ax.set_xlabel('time (s)')
                    ax.set_ylabel('FID (real)')
                    ax.grid('on')

            ax = plt.subplot(1, 2, 1)
            disp_fit(ax, self.data, params, self.seq, True, True, None, self._water_only, self.display_range_ppm)
            current_fit_time = time.time() - self._fit_time
            if(fit_terminated):
                ax.set_title("Fit terminated! %ds elapsed \n iteration #%d | residue = %.2E | R2 = %.2f/%.2f | FQN = %.2f" % (current_fit_time, self._model_call_count, current_err, rsq_t, rsq_f, fqn))
            else:
                ax.set_title("Running fit... %ds elapsed \n iteration #%d | residue = %.2E | R2 = %.2f/%.2f | FQN = %.2f" % (current_fit_time, self._model_call_count, current_err, rsq_t, rsq_f, fqn))
            fig.subplots_adjust(left=0.02, bottom=0.1, right=0.98, top=0.95, wspace=0.2, hspace=None)
            plt.pause(0.01)

            # save PNG
            # plt.savefig("%0.4d.png" % self._model_call_count)

            log.debug(str(self._model_call_count) + "th model call!")
            params.print(True, True)

        # fix for complex numbers
        diff_cmplx_odd_even = np.zeros(diff.size * 2, dtype=np.float64)
        diff_cmplx_odd_even[0::2] = diff.real
        diff_cmplx_odd_even[1::2] = diff.imag

        return(diff_cmplx_odd_even)

    def initialize(self):
        """Initialize the fit procedure."""
        log.info("initializing...")

        # do we need to adapt our ppm range?
        if(self.optim_ppm_range != self.seq.meta_bs.ppm_range):
            # we need to reinitialize the metabolite basis set
            self.seq.meta_bs.ppm_range = self.optim_ppm_range
            self.seq.meta_bs.initialize()
            # and rerun sequence !
            self.seq.initialize()

        # ckeck consistency lower<>init<>upper bounds
        log.info("checking consistency with parameter bounds...")
        for l in range(self.params_init.shape[0]):
            for c in range(self.params_init.shape[1]):
                # min vs max
                if(self.params_min[l, c] >= self.params_max[l, c] and self.params_init.linklock[l, c] != 1):
                    log.error(" One lower bound (%f) is equal/greater than a upper bound (%f) for [%s] at index (%d,%d)!" % (self.params_min[l, c], self.params_max[l, c], self.params_init.get_meta_names()[l], l, c))
                # min vs init
                if(self.params_min[l, c] >= self.params_init[l, c] and self.params_init.linklock[l, c] != 1):
                    log.error(" One initial value (%f) is equal/lower than the lower bound value (%f) for [%s] at index (%d,%d)!" % (self.params_init[l, c], self.params_min[l, c], self.params_init.get_meta_names()[l], l, c))
                # max vs init
                if(self.params_max[l, c] <= self.params_init[l, c] and self.params_init.linklock[l, c] != 1):
                    log.error(" One initial value (%f) is equal/greater than the upper bound value (%f) for [%s] at index (%d,%d)!" % (self.params_init[l, c], self.params_max[l, c], self.params_init.get_meta_names()[l], l, c))

        # checking if LL is not broken
        if(not self.params_init.check()):
            log.error("> the link-lock vector looks broken!")

        # checking if we are dealing with a reference scan fit (water only)
        where_LL = np.where(self.params_init.linklock[:, 0] == 0)
        self._water_only = (len(where_LL[0]) == 1 and where_LL[0][0] == xxx.m_Water)
        if(self._water_only):
            log.info("mmmh, looks like you are fitting a water reference scan!")

        # I am ready now
        self._ready = True

    def run(self):
        """
        Run the fit.

        Returns
        -------
        params_fit : params object
            Estimated metabolic parameters using the fitting algorithm
        optim_result : 'OptimizeResult' object returned by scipy.optimize.least_squares
            Numerical optimization parameters, includes the R^2 coefficients [rsq_t and rsq_f] and the FQN coefficient [fqn]
        """
        log.info_line________________________()
        log.info("running fit...")
        log.info_line________________________()

        # ready or not, here I come
        if(not self.ready):
            log.error("> This fit_tool object was not initialized!")

        # dummy full>free>full>free conversion to apply master/slave rules
        params_free = self.params_init.toFreeParams()
        params_full = self.params_init.toFullParams(params_free)
        params_free = params_full.toFreeParams()

        params_free = self.params_init.toFreeParams()
        params_min_free = self.params_min.toFreeParams()
        params_max_free = self.params_max.toFreeParams()
        self._model_call_count = 0

        # 3, 2, 1
        self._fit_time = time.time()
        self._model_call_count = 0
        self._cost_function = []

        # run it
        log.info("calling least-square optimizer...")
        if(self.optim_jacobian):
            optim_result = optimize.least_squares(self._minimizeThis, params_free, bounds=(params_min_free, params_max_free), jac=self._jac, xtol=self.optim_xtol, gtol=self.optim_gtol, ftol=self.optim_ftol)
        else:
            optim_result = optimize.least_squares(self._minimizeThis, params_free, bounds=(params_min_free, params_max_free), xtol=self.optim_xtol, gtol=self.optim_gtol, ftol=self.optim_ftol)

        # nth final display
        self._minimizeThis(optim_result.x, True)
        self._fit_time = time.time() - self._fit_time

        # results
        log.info("optimization terminated!")
        # final pars
        params_fit = self.params_init.copy()
        params_fit = params_fit.toFullParams(optim_result.x)
        # final CRBs
        params_CRBs_abs, _ = self.get_CRBs(params_fit)
        params_fit._errors = params_CRBs_abs

        # add fit criteria to optim_result
        rsq_t, rsq_f = self.get_Rsq(params_fit)
        optim_result.rsq_t = rsq_t
        optim_result.rsq_f = rsq_f
        optim_result.fqn = self.get_FQN(params_fit)

        # switch back to not ready
        self._ready = False

        return(params_fit, optim_result)

    def get_CRBs(self, params):
        """
        Return the Cramér-Rao Lower Bounds.

        Parameters
        ----------
        params : params object
            Array of simulation parameters

        Returns
        -------
        params_CRBs_abs : params object
            Array of simulation parameters, absolute CRBs
        params_CRBs_rel : params object
            Array of simulation parameters, relative CRBs (%)
        """
        # calling jacobian model
        j_full = self.seq._jac(params)

        # reducing it to free params
        j_free = j_full[params.linklock <= 0, :]
        j_free = np.real(j_free)

        # Fisher (number of free metabolites*number of params, time points)
        F = np.dot(j_free, np.transpose(j_free))

        # covariance
        try:
            covar = np.linalg.inv(np.real(F))
        except:
            # it happens that F is a singular matrix
            # I believe it is just bad luck but I am not really sure how to handle this error...
            # for now, I just had some little salt to F
            log.warning("weird singular matrix bug while inverting Fisher matrix during CRB calculation!")
            F[0, 0] = F[0, 0] + 1e-6
            covar = np.linalg.inv(np.real(F))

        # CRBs
        CRBs = np.zeros(covar.shape[0])
        for i in range(covar.shape[0]):
            CRBs[i] = np.sqrt(covar[i, i])

        # original noise level before apodization
        noise_std = self.data.noise_level

        # normalize to noise: absolute CRBs
        CRBs_abs = CRBs * noise_std

        # normalize to parameter values: relative CRBs (%)
        params_free = params.toFreeParams()
        CRBs_rel = (CRBs_abs / params_free) * 100.0

        # put back to full parameter form
        params_CRBs_abs = self.params_init.copy()
        params_CRBs_abs = params_CRBs_abs.toFullParams(CRBs_abs)
        params_CRBs_rel = self.params_init.copy()
        params_CRBs_rel = params_CRBs_rel.toFullParams(CRBs_rel)

        return(params_CRBs_abs, params_CRBs_rel)

    def get_Rsq(self, params):
        """
        Return the R^2 coefficient in time and frequency domain.

        Parameters
        ----------
        params : params object
            Array of simulation parameters

        Returns
        -------
        r2_t : float
            R^2 correlation coefficient between data and fit in TIME DOMAIN
        r2_f: float
            R^2 correlation coefficient between data and fit in FREQUENCY DOMAIN
        """
        # the model
        mod = self.seq._model(params)

        # and the R coeff in TIME DOMAIN
        c = np.corrcoef(self.data, mod)
        r_t = c[0, 1]
        r2_t = np.abs(r_t**2)

        # and the R coeff in FREQUENCY DOMAIN
        c = np.corrcoef(self.data.spectrum(), mod.spectrum())
        r_f = c[0, 1]
        r2_f = np.abs(r_f**2)

        return(r2_t, r2_f)

    def get_FQN(self, params):
        """
        Return the FQN coefficient. See https://doi.org/10.1002/nbm.4257 : Quantitatively, this can be expressed using the fit quality number (FQN), which is the ratio of the variance in the fit residual divided by the variance in the pure spectral noise. For an ideal fit, the FQN should be close to 1.0, and the FQN/SNR ratio should be much less than 1.

        Parameters
        ----------
        params : params object
            Array of simulation parameters

        Returns
        -------
        fqn : float
            FQN coefficient of the fit
        """
        # residue
        mod = self.seq._model(params)
        diff = mod - self.data

        # our model is in TIME DOMAIN
        # but FQN only makes sense in FREQUENCY DOMAIN
        # because the time domain data here is very probalby apodized
        # the pure noise variance is therefore impossible to estimate
        # that's fine with me, let's do it FREQUENCY DOMAIN

        # estimate noise variance in user specified spectral region
        ppm = self.data.frequency_axis_ppm()
        sf = np.squeeze(self.data.spectrum())
        sf_analyze = np.real(sf)
        ippm_noise_range = (self.fqn_noise_range[0] < ppm) & (ppm < self.fqn_noise_range[1])
        data_noise_var = np.std(sf_analyze[ippm_noise_range])

        # variance of fit residual
        residue_var = np.std(np.real(diff.spectrum()))

        fqn = residue_var / data_noise_var
        return(fqn)

    def disp_corr_mat(self, ax, params, mIndex_list=None, pIndex_list=[xxx.p_cm], color_bar=False):
        """
        Return and display the parameter correlation matrix.

        Parameters
        ----------
        ax : matplotlib axis
            Axis to use for plotting
        params : params object
            Array of simulation parameters
        mIndex_list : list of int
            Metabolite indexes to display in bargraph
        pIndex_list : list of int
            Parameter indexes to display in bargraph
        color_bar : boolean
            Add a color (True) or not (False)

        Returns
        -------
        params_cor : numpy array
            Correlation matrix
        """
        # calling jacobian model
        j_full = self.seq._jac(params)

        # probably do not want to look at the full matrix but only some metabolites/parameters
        mat_linklock_mIndex = np.full(params.linklock.shape, False, dtype=bool)
        mat_linklock_mIndex[mIndex_list, :] = True
        mat_linklock_pIndex = np.full(params.linklock.shape, False, dtype=bool)
        mat_linklock_pIndex[:, pIndex_list] = True
        mat_linklock_mask = np.logical_and(np.logical_and(mat_linklock_mIndex, mat_linklock_pIndex), (params.linklock <= 0))

        # build labels vector the same way
        cor_labels = []
        meta_names = params.get_meta_names()
        for i, m in enumerate(meta_names):
            for j, p in enumerate(PARS_SHORT_NAMES):
                if(mat_linklock_mask[i, j]):
                    cor_labels.append(m + "|" + p)

        if(np.all(mat_linklock_mask == False)):
            log.error(">> mrs.fit_tool.disp_CorrMat: nothing to display, please adjust the mIndex_list/pIndex_list parameters!")

        # reducing it to free params
        j_free = j_full[mat_linklock_mask, :]
        j_free = np.real(j_free)

        # Fisher (number of free metabolites*number of params, time points)
        F = np.dot(j_free, np.transpose(j_free))

        # covariance
        covar = np.linalg.inv(np.real(F))

        # parameter correlation
        params_cor = np.zeros(covar.shape)
        for i in range(params_cor.shape[0]):
            for j in range(params_cor.shape[1]):
                params_cor[i, j] = covar[i, j] / (np.sqrt(covar[i, i]) * np.sqrt(covar[j, j]))

        ax.cla()
        im = ax.imshow(params_cor, clim=(-1, +1))
        ax.set_xticks(np.arange(len(cor_labels)))
        ax.set_yticks(np.arange(len(cor_labels)))
        ax.set_xticklabels(cor_labels, rotation=90, fontsize=6)
        ax.set_yticklabels(cor_labels, fontsize=6)
        if(color_bar):
            plt.colorbar(im)

        return(params_cor)


def disp_bargraph(ax, params_val_list, params_leg_list, colored_LLs=True, bMM=False, bWater=False, bFitMode=False, pIndex=xxx.p_cm, mIndex_list=None, width=0.3):
    """
    Plot a bargraph of concentrations.

    Parameters
    ----------
    ax : matplotlib axis
        Axis to use for plotting
    params_val_list : list of params objects
        List of parameter arrays to display as bars
    params_leg_list : list of params objects
        List of legend caption for each bar
    colored_LLs : boolean
        Shows link-lock connections using colors (True) or not (False)
    bMM : boolean
        Includes macromolecular parameters (True) or not (False)
    bWater : boolean
        Includes water parameters (True) or not (False)
    bFitMode : boolean
        If (True), plot the first two bars in black (bounds), the last one with colors (fit), at the same x, no legends
    pIndex : int
        Parameter index to display in bargraph
    mIndex_list : list of int
        Metabolite indexes to display in bargraph
    width : float (optional)
        Bar width
    """
    # LLcolor_names=[name for name in mcd.XKCD_COLORS]
    LLcolor_names = ['xkcd:grey', 'xkcd:blue', 'xkcd:red', 'xkcd:violet', 'xkcd:green',
                     'xkcd:goldenrod', 'xkcd:crimson', 'xkcd:salmon', 'xkcd:wheat', 'xkcd:lightblue']

    # init
    lbl = PARS_LONG_NAMES[pIndex]

    # build logical mask
    meta_mask = np.full(params_val_list[0].shape, True, dtype=bool)

    # filter out the fixed LLs
    meta_mask[params_val_list[0].linklock[:, xxx.p_cm] == 1, :] = False

    # filter out MMs?
    if(not bMM):
        meta_mask[xxx.m_All_MMs, :] = False

    # filter out water?
    if(not bWater):
        meta_mask[xxx.m_Water, :] = False

    if(mIndex_list is not None):
        meta_mask_m = meta_mask.copy()
        meta_mask_m[:] = False
        meta_mask_m[mIndex_list, :] = True
        meta_mask = np.logical_and(meta_mask, meta_mask_m)

    # parameter filter
    meta_mask_p = meta_mask.copy()
    meta_mask_p[:] = False
    meta_mask_p[:, pIndex] = True
    meta_mask = np.logical_and(meta_mask, meta_mask_p)

    # check that we still have something to display
    if(np.all(meta_mask == False)):
        log.error("nothing to display! Please check disp_bargraph parameters like pIndex, mIndex_list, bWater and bMM!")

    # filter data now
    params_val_list = [p[meta_mask] for p in params_val_list]
    params_LL_list = [p.linklock[meta_mask] for p in params_val_list]
    params_std_list = [p._errors[meta_mask] for p in params_val_list]
    meta_names = params_val_list[0].get_meta_names()
    meta_names = [meta_names[i] for i in range(len(meta_names)) if meta_mask[i, pIndex]]

    # how many metabolites to display?
    n = np.sum(meta_mask[:])

    # prepare bars
    nBars = len(params_val_list)
    pos_bars = np.arange(n)
    pos_shift = np.linspace(-width * (nBars - 1) / 2.0, +width * (nBars - 1) / 2.0, nBars)
    pv_min = params_val_list[0].min()
    pv_max = params_val_list[0].max()
    if(bFitMode):
        pos_shift = pos_shift * 0.0

    # iterate in lists and plot bars yeah
    for (k, pv, pll, ps, pl, pos) in zip(range(len(params_val_list)), params_val_list, params_LL_list, params_std_list, params_leg_list, pos_shift):
        bar_list = ax.bar(pos_bars + pos, pv, width, yerr=ps, label=pl)

        # record min/max for later use
        if(np.min(pv) < pv_min):
            pv_min = np.min(pv)
        if(np.max(pv) > pv_max):
            pv_max = np.max(pv)

        for im_bar, this_bar, this_LL in zip(range(len(bar_list)), bar_list, pll):
            if(colored_LLs):
                this_LLcolor_name = LLcolor_names[np.mod(k + int(np.abs(this_LL)), len(LLcolor_names))]
            else:
                this_LLcolor_name = LLcolor_names[np.mod(k, len(LLcolor_names))]

            this_LLcolor = mcd.XKCD_COLORS[this_LLcolor_name].upper()
            if(bFitMode and k < 2):
                this_bar.set_color('grey')
            else:
                this_bar.set_color(this_LLcolor)

    ax.set_ylabel(lbl)
    ax.set_ylim([pv_min, pv_max])
    ax.set_xticks(pos_bars)
    ax.set_xticklabels(meta_names, rotation=90)
    ax.grid('on')
    if(not bFitMode):
        ax.legend()


def disp_fit(ax, data, params, seq, LL_exluding=True, LL_merging=False, mIndex_list=None, water_only=False, display_range=[1, 5]):
    """
    Plot a bargraph of concentrations.

    Parameters
    ----------
    ax : matplotlib axis
        Axis to use for plotting
    params : params object
        List of parameter arrays to display as bars
    seq : mrs_sequence
        Virtual MRS sequence object
    LL_excluding : boolean
        Shows only free or master spectra (True) or not (False)
    LL_merging : boolean
        Shows concentration link-lock connections by merging spectra (True) or not (False)
    mIndex_list : list
        Indexes of metabolite to display
    water_only : boolean
        If only water is displayed
    display_range : list [2]
            Range in ppm used for display
    """
    # dummy full>free>full conversion to apply master/slave rules
    params_free = params.toFreeParams()
    params_full = params.toFullParams(params_free)

    # the data
    ax.plot(data.frequency_axis_ppm(), data.spectrum().real, 'k-', linewidth=0.5, label='data')
    # the full model
    mod = seq._model(params_full)
    ax.plot(mod.frequency_axis_ppm(), mod.spectrum().real, 'k-', linewidth=2, label='model')
    # the residue
    diff = data - mod
    ax.plot(diff.frequency_axis_ppm(), diff.spectrum().real, 'g-', linewidth=0.25, label='residue')

    # trying to set ylim smart way
    ymax = np.max(mod.spectrum().real)
    # if not only water and water exists
    try:
        if(not water_only):
            # generate a model signal, water set to 0 to set later ylim
            params_full_nowater = params_full.copy()
            params_full_nowater[xxx.m_Water, xxx.p_cm] = 0.0
            mod_nowater = seq._model(params_full_nowater)
            ymax = np.max(mod_nowater.spectrum().real) * 1.1
    except NameError:
        pass

    # metabolite list
    if(mIndex_list is None):
        meta_ind_list = [*range(xxx.n_MBs)]
    else:
        meta_ind_list = mIndex_list

    meta_names = params.get_meta_names()

    # filter out fixed LL?
    if(LL_exluding):
        # remove fixed parameters (LL==1)
        meta_ind_list_without_fixedLLs = []
        for im in meta_ind_list:
            if(params_full.linklock[im, xxx.p_cm] < 1):
                meta_ind_list_without_fixedLLs.append(im)
        meta_ind_list = meta_ind_list_without_fixedLLs

    # calculate a nice gap
    n_meta_to_display = len(meta_ind_list) + 1  # 1 line for MMs
    if(LL_merging):
        n_meta_to_display - np.unique(params_full.linklock[:, xxx.p_cm])
    ygap = ymax * (60 - n_meta_to_display) * 0.005

    # init metabolite-per-metabolite plot
    ignore_those_mLLs = []
    p_allzeros = params_full.copy()
    p_allzeros[:, xxx.p_cm] = 0.0
    koffset = 0
    for im, mLL in zip(meta_ind_list, params_full.linklock[meta_ind_list, xxx.p_cm]):
        # check if we should plot this guy
        mLL_abs = np.abs(mLL)
        if(mLL_abs not in ignore_those_mLLs):
            if(LL_merging and mLL_abs > 1):
                # find which metabolites to merge (NAA_CH3 and NAA_CH2 for example)
                this_params_LL_abs = np.abs(params_full.linklock[:, xxx.p_cm])
                this_meta_mask = (this_params_LL_abs == mLL_abs)
                this_meta_mask = np.tile(this_meta_mask, (1, 1))
                this_meta_mask = np.transpose(this_meta_mask)
                this_meta_mask = np.tile(this_meta_mask, (1, 4))
                # remember for next time we meet one of those metabolites
                ignore_those_mLLs.append(mLL_abs)
            else:
                this_meta_mask = np.full(params_full.shape, False, dtype=bool)
                this_meta_mask[im, :] = True

            # prepare metabolite name
            this_meta_name = ""
            for mName, mLL in zip(meta_names, this_meta_mask[:, xxx.p_cm]):
                if(mLL):
                    this_meta_name = this_meta_name + mName + " & "
            this_meta_name = this_meta_name[:-3]

            # prepare parameters
            this_params = p_allzeros.copy()
            this_params[this_meta_mask] = params_full[this_meta_mask]

            # call model
            s_single_meta = seq._model(this_params)

            # display
            ax.plot(s_single_meta.frequency_axis_ppm(), s_single_meta.spectrum().real - ygap * (koffset + 1), 'r-')
            ax.text(display_range[1] - 0.5, -ygap * (koffset + 1) + ygap / 5, this_meta_name)
            koffset = koffset + 1

    # prepare macromolecules parameters
    this_params = p_allzeros.copy()
    this_params[xxx.m_All_MMs, :] = params_full[xxx.m_All_MMs, :]

    # call model
    s_MMs = seq._model(this_params)

    # plot the final MM baseline
    ax.plot(s_MMs.frequency_axis_ppm(), s_MMs.spectrum().real - ygap * (koffset + 2), 'r-')

    # and now break down the MM baseline into the gaussian components
    for im in xxx.m_All_MMs:
        # prepare parameters for this MM
        this_params = p_allzeros.copy()
        this_params[im, :] = params_full[im, :]

        # call model
        s_single_MM = seq._model(this_params)

        # plot the a single MM keeping the same offset
        ax.plot(s_single_MM.frequency_axis_ppm(), s_single_MM.spectrum().real - ygap * (koffset + 2), 'r--')

    ax.text(display_range[1] - 0.5, -ygap * (koffset + 2) + ygap / 5, 'MM baseline')

    # finalize the plot
    ax.set_xticks(np.arange(-1, 10, 0.5))
    ax.set_xlim(display_range[1], display_range[0])
    # ylim trick: tricky
    ymax = min(ymax, np.max(mod.spectrum().real) * 1.1)
    ax.set_ylim([-ygap * (koffset + 3), ymax])
    ax.set_yticks([])
    ax.set_xlabel('chemical shift (ppm)')
    ax.grid('on')
    ax.legend(loc="upper right")
