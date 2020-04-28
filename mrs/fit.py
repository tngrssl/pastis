#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A fit_tool class that deals with the quantification of MRS data based on the scipy's least_squares optimisation function.

@author: Tangi Roussel
"""

import numpy as np
import mrs.sim as sim
import mrs.log as log
import mrs.aliases as xxx
import matplotlib.pylab as plt
import scipy.optimize as optimize
from enum import Enum
import time

import pdb


class fit_plot_type(Enum):
    """The enum fit_plot_type describes the type of plots that can be displayed real-time during the fit."""

    BARGRAPH_CM = 1
    BARGRAPH_DD = 2
    BARGRAPH_DF = 3
    BARGRAPH_DP = 4
    COST_FUNCTION = 5
    CORR_MATRIX = 6
    TIME_DOMAIN = 7


class prefit_tool:
    """The prefit_tool class is used to perfom area integration of peaks in the spectrum."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("You are trying to dynamically create a new attribute (%s) to this object and that is not cool! I will not let you do that because I believe it is a bad habit and can lead to terrible bugs. A clean way of doing this is to initialize your attribute (%s) in the __init__ method of this class. Bisou, bye :)" % (key, key))
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
        self.area_integration_peaks = [xxx.m_Cho_CH3, xxx.m_Cr_CH3, xxx.m_NAA_CH3]
        # ppm range to look for peak maximum
        self.area_integration_peak_search_range = 0.4  # ppm
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
        self.display_range_ppm = [1, 6]  # ppm

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
        self._peak_base_width = []

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
                raise Exception("> one of the metabolites you chose for peak integration does not have a singlet, aborting area integration...")
            else:
                # measure linewidth of the peak
                aipsr = self.area_integration_peak_search_range / 2.0
                this_lw = self.data.analyze_linewidth([this_meta_ppm - aipsr / 2.0, this_meta_ppm + aipsr / 2.0], this_meta_name, False, True)
                # store
                self._peak_names.append(this_meta_name)
                self._peak_ppms.append(this_meta_ppm)
                self._peak_nprots.append(this_meta_nprots)
                # according to https://doi.org/10.1002/nbm.4257
                # peak integration area should be 2 * peak linewidth
                self._peak_base_width.append(this_lw * 2.0)

        # I am ready now
        self._ready = True

    def run(self):
        """Run the area integration pipeline."""
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
        sf_abs = np.abs(sf)

        if(self.display_enable):
            fig = plt.figure(self.display_fig_index)
            fig.clf()
            axs = fig.subplots()
            fig.canvas.set_window_title("mrs.fit.prefit_tool.run")
            axs.plot(self.data.frequency_axis_ppm(), self.data.spectrum().real, "k-", label="data")

        # now for each peak, integrate the area
        self._peak_areas = []
        self._peak_areas_norm = []
        p = sim.params(self.meta_bs)
        p[:] = 0.0

        log.info("integrating peak for...")
        log.info_line________________________()
        for (this_peak_index, this_peak_name, this_peak_ppm, this_peak_np, this_peak_basewidth) in zip(self.area_integration_peaks, self._peak_names, self._peak_ppms, self._peak_nprots, self._peak_base_width):
            log.info("[%s] theoretically at %.2fppm in theory" % (this_peak_name, this_peak_ppm))
            # first find closest peak in range
            peakwidth_ppm = this_peak_basewidth / self.data.f0
            ippm_peak_range = ((this_peak_ppm - self.area_integration_peak_search_range) > ppm) | (ppm > (this_peak_ppm + self.area_integration_peak_search_range))
            sf_masked = sf_abs.copy()
            sf_masked[ippm_peak_range] = 0
            ippm_peak = np.argmax(sf_masked)
            if(ippm_peak == 0):
                log.error("no peak found in specified ppm range or badly phased data!")
            log.debug("found it in the spectrum at %.2fppm" % ppm[ippm_peak])

            # integrate area
            ippm_peak_range = (ppm > (ppm[ippm_peak] - peakwidth_ppm / 2)) & (ppm < (ppm[ippm_peak] + peakwidth_ppm / 2))
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
            p[this_peak_index, 0] = this_peak_area / this_peak_np

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

        # return a params vector
        return(p)


class fit_tool:
    """The fit_tool class is used to adjust the MRS model implemented in the mrs.sim.metabolite_database class on some experimentally acquired data."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("You are trying to dynamically create a new attribute (%s) to this object and that is not cool! I will not let you do that because I believe it is a bad habit and can lead to terrible bugs. A clean way of doing this is to initialize your attribute (%s) in the __init__ method of this class. Bisou, bye :)" % (key, key))
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
        # should we nomalize the concentrations displayed during the fit?
        self.display_normalized = False
        # if yes, use those T2 and T1 values from metabolite database
        self.display_normalized_T2_key = "T2s_human_brain_7T"
        self.display_normalized_T1_key = "T1s_human_brain_7T"
        # if yes, use those parameters estimated on reference non water suppressed data
        self.display_normalized_abs_ref_params = sim.params(self._meta_bs)
        # and assume a water concentration in mmol/L
        self.display_normalized_abs_water_concentration = 55000.0

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

            # prepare bar values
            if(self.display_normalized):
                params_T2norm = params.correct_T2s(self.data.te, self.display_normalized_T2_key)
                params_T2norm_T1norm = params_T2norm.correct_T1s(self.data.tr, self.display_normalized_T1_key)
                params_T2norm_T1norm_abs = params_T2norm_T1norm.correct_absolute(self.display_normalized_abs_ref_params, self.display_normalized_abs_water_concentration)
                params_val_list = [self.params_min, self.params_max, params_T2norm_T1norm_abs]
            else:
                params_val_list = [self.params_min, self.params_max, params]

            # real-time CRBs
            if(self.display_CRBs):
                params_CRBs_abs, params_CRBs_rel = self.get_CRBs(params)
            else:
                params_CRBs_abs = params * 0.0

            params_std_list = [self.params_min * 0.0, self.params_max * 0.0, params_CRBs_abs]
            params_leg_list = ['lower bounds', 'upper bounds', 'current']

            # let's display the subplots as required by self.display_subplots_types
            axs_list = [axs[0, 2], axs[0, 3], axs[1, 2], axs[1, 3]]
            for ax, plt_type in zip(axs_list, self.display_subplots_types):
                if(plt_type == fit_plot_type.BARGRAPH_CM):
                    sim.disp_bargraph(ax, params_val_list, params_std_list, params_leg_list, True, True, self._water_only, True, xxx.p_cm)
                elif(plt_type == fit_plot_type.BARGRAPH_DD):
                    sim.disp_bargraph(ax, params_val_list, params_std_list, params_leg_list, True, True, self._water_only, True, xxx.p_dd)
                elif(plt_type == fit_plot_type.BARGRAPH_DF):
                    sim.disp_bargraph(ax, params_val_list, params_std_list, params_leg_list, True, True, self._water_only, True, xxx.p_df)
                elif(plt_type == fit_plot_type.BARGRAPH_DP):
                    sim.disp_bargraph(ax, params_val_list, params_std_list, params_leg_list, True, True, self._water_only, True, xxx.p_dp)
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
            sim.disp_fit(ax, self.data, params, self.seq, True, True, None, self.display_range_ppm)
            current_fit_time = time.time() - self._fit_time
            if(fit_terminated):
                ax.set_title("Fit terminated! %ds elapsed | iteration #%d | residue = %.2E" % (current_fit_time, self._model_call_count, current_err))
            else:
                ax.set_title("Running fit... %ds elapsed | iteration #%d | residue = %.2E" % (current_fit_time, self._model_call_count, current_err))
            fig.subplots_adjust()

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
                    raise Exception(" One lower bound (%f) is equal/greater than a upper bound (%f) for [%s] at index (%d,%d)!" % (
                        self.params_min[l, c], self.params_max[l, c], self.params_init.get_meta_names()[l], l, c))
                # min vs init
                if(self.params_min[l, c] >= self.params_init[l, c] and self.params_init.linklock[l, c] != 1):
                    raise Exception(" One initial value (%f) is equal/lower than the lower bound value (%f) for [%s] at index (%d,%d)!" % (
                        self.params_init[l, c], self.params_min[l, c], self.params_init.get_meta_names()[l], l, c))
                # max vs init
                if(self.params_max[l, c] <= self.params_init[l, c] and self.params_init.linklock[l, c] != 1):
                    raise Exception(" One initial value (%f) is equal/greater than the upper bound value (%f) for [%s] at index (%d,%d)!" % (
                        self.params_init[l, c], self.params_max[l, c], self.params_init.get_meta_names()[l], l, c))

        # checking if LL is not broken
        if(not self.params_init.check()):
            raise Exception("> the link-lock vector looks broken!")

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
        params_CRBs_abs : params object
            Estimated absolute Cramér-Rao Lower Bounds for each parameter
        params_CRBs_rel : params object
            Estimated relative Cramér-Rao Lower Bounds for each parameter
        optim_result : 'OptimizeResult' object returned by scipy.optimize.least_squares
            Numerical optimization parameters
        """
        log.info_line________________________()
        log.info("running fit...")
        log.info_line________________________()

        # ready or not, here I come
        if(not self.ready):
            raise Exception("> This fit_tool object was not initialized!")

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
        params_CRBs_abs, params_CRBs_rel = self.get_CRBs(params_fit)

        return(params_fit, params_CRBs_abs, params_CRBs_rel, optim_result)

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
            for j, p in enumerate(sim.PARS_SHORT_NAMES):
                if(mat_linklock_mask[i, j]):
                    cor_labels.append(m + "|" + p)

        if(np.all(mat_linklock_mask==False)):
            raise Exception(">> mrs.fit_tool.disp_CorrMat: nothing to display, please adjust the mIndex_list/pIndex_list parameters!")

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
        ax.set_xticklabels(cor_labels, rotation=90, fontsize=9)
        ax.set_yticklabels(cor_labels, fontsize=9)
        if(color_bar):
            plt.colorbar(im)

        return(params_cor)
