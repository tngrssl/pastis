#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A fit_tool class that deals with the quantification of MRS data based on the scipy's least_squares optimisation function.

@author: Tangi Roussel
"""

import numpy as np
import pandas as pd
import suspect
from mrs import sim
from mrs import log
from mrs import aliases as xxx
import matplotlib.pylab as plt
import scipy.optimize as optimize
import matplotlib._color_data as mcd
from enum import Enum
import os
import time
import hashlib
import copy
from parse import parse

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


class fit_adjust_metabolites_mode(Enum):
    """The enum fit_adjust_metabolite_mode describes how the algorithm should look for fittable metabolites. See the method fit_pastis.adjust_metabolites() for more info."""

    NONE = 0
    OPTIMISTIC = 1
    AVERAGE = 2
    PESSIMISTIC = 3


class fit_adjust_metabolites_when(Enum):
    """The enum fit_adjust_metabolite_when describes when the algorithm should look for fittable metabolites."""

    NEVER = 0
    BEFORE_SECOND_FIT_ONLY = 1
    BEFORE_EACH_FIT = 2


class fit_tool():
    """The fit_tool class is a virtual class, mother of other classes used for fitting."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, data=None):
        """Construct a fit_tool object.

        Parameters
        ----------
        data : MRSData2 objet
            Data to fit
        """
        # name of the fit (fit approach for example)
        self.name = "fit"
        # data signal to process
        self.data = data

    def copy(self):
        """Copy method."""
        obj = copy.copy(self)
        if(self.data is not None):
            obj.data = self.data.copy()

        return(obj)


class fit_lcmodel(fit_tool):
    """The fit_lcmodel class is used to fit acquired MRS data using the blackbox now opensource old FORTRAN LCModel reference algorithm which has been used for the past 30 years."""

    def __init__(self, data=None, meta_bs=None):
        """
        Construct a fit_lcmodel object.

        Parameters
        ----------
        data : MRSData2 objet
            Data to fit
        meta_bs : sim.metabolite_basis_set object
            Metabolite db to use to fit. If None, use meta_bs from sequence. If sequence is None, use default meta_bs.
        """
        super().__init__(data)

        # some LCModel parameters
        self.lcmodel_executable_fullpath = "/home/tangir/crmbm/soft/lcmodel/lcm-64/bin/lcmodel"
        self.lcmodel_rawfile_fullpath = "/home/tangir/crmbm/soft/lcmodel/data/data.raw"
        self.lcmodel_filbas = "/home/tangir/crmbm/soft/lcmodel/metabolite_basis_sets/7t/gamma_press_te{}_7t_v1.basis"
        self.lcmodel_lcsv = True

        # data
        self.data = data

        # default metabolite basis set used for simulation
        if(meta_bs is not None):
            self.meta_bs = meta_bs
        else:
            self.meta_bs = sim.metabolite_basis_set()

        # display stuff
        self.display_enable = True
        self.display_pdf_software = "evince"

        # final results parameters
        self.params_fit = None

        # freeze
        self.__isfrozen = True

    def copy(self):
        """Copy method."""
        obj = super().copy()

        return(obj)

    def _adjust_filbas(self):
        """Adjust the metabolite basis set according to te."""
        log.info("adjusting metabolite basis set TE...")

        # hold debugs logs during parse
        log.pause()

        # scan basis set files
        metabs_folder = "/".join(self.lcmodel_filbas.split("/")[0:-1])
        metabs_filename_mask = self.lcmodel_filbas.split("/")[-1]
        metabs_files_dict = {}
        for _, _, metabs_files in os.walk(metabs_folder):
            for this_metabs_filename in metabs_files:
                parse_res = parse(metabs_filename_mask, this_metabs_filename)
                if(parse_res is not None):
                    this_metabs_te = parse_res[0]
                    metabs_files_dict[this_metabs_te] = metabs_folder + "/" + this_metabs_filename

        # put back log
        log.resume()

        # find metabs with a TE close to the data's TE
        metabs_te_arr = np.array(list(metabs_files_dict.keys())).astype(float)
        te_diff = np.abs(metabs_te_arr - self.data.te)
        ite_mindiff = np.argmin(te_diff)
        optim_te = list(metabs_files_dict.keys())[ite_mindiff]

        # change metabolite basis set file
        self.lcmodel_filbas = metabs_files_dict[optim_te]
        log.debug("found that it is better to use the TE=%sms metabolite basis set file!" % optim_te)

    def run(self):
        """Call LCModel using suspect (https://suspect.readthedocs.io/en/latest/notebooks/tut04_quant.html)."""
        log.info("running LCModel fit...")

        # check that we have everything we need
        if(self.data is None):
            log.error("canceling fit, no data was given to this fit object!")

        # adjust metabolite basis set file
        self._adjust_filbas()

        # create a parameters dictionary to set the basis set to use
        params = {
            "FILBAS": self.lcmodel_filbas,
            "LCSV": self.lcmodel_lcsv
        }
        # call suspect to prepare LCModel files
        suspect.io.lcmodel.write_all_files(self.lcmodel_rawfile_fullpath, self.data, self.data.data_ref, params=params)

        # guess CONTROL file name
        fullpath_noext, ext = os.path.splitext(self.lcmodel_rawfile_fullpath)
        controlfile_fullpath = fullpath_noext + "_sl0.CONTROL"

        # LCModel is now free, replace KEY by 210387309
        with open(controlfile_fullpath) as f:
            controlfile_data = f.read().replace('KEY = 123456789', 'KEY = 210387309')

        with open(controlfile_fullpath, "w") as f:
            f.write(controlfile_data)

        # make system call to LCModel
        system_res = os.system(self.lcmodel_executable_fullpath + "<" + controlfile_fullpath)
        log.debug("LCmodel binary system call returned " + str(system_res))

        log.debug("extracting LCmodel results from CSV...")
        csvfile_fullpath = fullpath_noext + ".CSV"
        df = pd.read_csv(csvfile_fullpath)

        # build params vector from this dataframe
        self.params_fit = sim.params(self.meta_bs)
        meta_names = self.params_fit.get_meta_names()
        meta_names_lcmodel = self.params_fit.get_meta_names(True)

        # make everything lowercase
        meta_names = [m.lower() for m in meta_names]
        meta_names_lcmodel = [m.lower() if (m is not None) else None for m in meta_names_lcmodel]
        df.columns = df.columns.str.strip()
        df.columns = df.columns.str.lower()

        for im, (m, mlc) in enumerate(zip(meta_names, meta_names_lcmodel)):
            if(mlc is None):
                # concentration
                self.params_fit[im, xxx.p_cm] = 0.0
            else:
                # concentration
                self.params_fit[im, xxx.p_cm] = df[mlc].iloc[0]
                # abs CRB error
                self.params_fit._errors[im, xxx.p_cm] = df[mlc].iloc[0] * df[mlc + " %sd"].iloc[0] / 100

        # display if needed
        if(self.display_enable):
            log.debug("displaying LCModel results...")
            psfile_fullpath = fullpath_noext + ".PS"
            pdffile_fullpath = fullpath_noext + ".pdf"
            os.system("ps2pdf " + psfile_fullpath + " " + pdffile_fullpath)
            os.system(self.display_pdf_software + " " + pdffile_fullpath + " &")

    def get_hash(self):
        """
        Generate a hash of this fit object.

        Returns
        -------
        h : string
            Hash code of this strategy. Usefull later when dealing with databases...
        """
        bytes_to_hash = "lcmodel".encode()

        h = hashlib.md5(bytes_to_hash)

        return(h.hexdigest())

    def to_dataframe(self, prefix_str="fit_"):
        """
        Convert the object's attributes to dataframe. Can include the object itself.

        Parameters
        ----------
        prefix_str : string
            Prefix string to add to column names

        Returns
        -------
        df : Dataframe
            With all the datasets and parameters from this fit object
        """
        log.debug("converting to dataframe...")

        # use pandas json_normalize function to flatten all nested stuff (magic!)
        # this takes care of the job attribute too
        df_attr = pd.json_normalize(vars(self), sep='_')

        # and need a specific call for the MRSData2 data
        df_data = self.data.to_dataframe(True, "data_")
        del df_attr["data"]

        df_params_fit = self.params_fit.to_dataframe("params_fit_")
        df_attr = df_attr.rename(columns = {'params_fit': 'params_fit_obj'})

        # append columns
        df = pd.concat([df_attr.reset_index(),
                        df_data.reset_index(),
                        df_params_fit], axis=1)

        # add fit hash
        df["fit_hash"] = self.get_hash()

        # remove index column
        df = df.reset_index()
        del df["index"]

        # add prefix
        df = df.add_prefix(prefix_str)

        return(df)


class fit_pastis(fit_tool):
    """The fit_pastis class is used to fit acquired MRS data using a linear combination algorithm based on QUEST, see Ratiney et al. NMR Biomed, 2005."""

    def __init__(self, data=None, sequence=None, meta_bs=None):
        """
        Construct a fit_pastis object.

        Parameters
        ----------
        data : MRSData2 objet
            Data to fit
        sequence : sim.mrs_sequence object
            Sequence to use to fit. If None, use sequence from data
        meta_bs : sim.metabolite_basis_set object
            Metabolite db to use to fit. If None, use meta_bs from sequence. If sequence is None, use default meta_bs.
        """
        super().__init__(data)

        # default sequence used for simulation (if none, will take sequence from data)
        self.sequence = sequence

        # default metabolite basis set used for simulation (if none, will take sequence from data)
        if(meta_bs is not None):
            self.meta_bs = meta_bs
        else:
            if(self.sequence is not None):
                if(self.sequence.meta_bs is not None):
                    self.meta_bs = self.sequence.meta_bs
                else:
                    # default meta_bs
                    self.meta_bs = sim.metabolite_basis_set()
            else:
                # default meta_bs
                self.meta_bs = sim.metabolite_basis_set()

        # list of metabolites to fit
        self.metabolites = None
        # list of metabolites to integrate
        self.metabolites_area_integration = [xxx.m_Water]

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
        # with early minimal damping
        self.params_init[xxx.m_All_MBs, xxx.p_dd] = self.params_min[xxx.m_All_MBs, xxx.p_dd] * 1.1

        # linklock array
        self.params_linklock = self.params_init._linklock.copy()
        self._set_unique_linklock()

        # --- final results parameters ---
        self.params_fit = None
        self.params_area = None
        self.params_area_pnorm = None
        self.optim_results = None

        # private option to know if we are dealing with a water fit
        self._water_only = False

        # --- optimization options ---
        # should we use jacobin information during the fit or not?
        self.optim_jacobian = True
        # ppm range (By default, use metabolite_basis_set ppm range)
        self.optim_ppm_range = self.meta_bs.ppm_range
        # stop fit if parameter changes go below this tolerance
        self.optim_xtol = 1e-9
        # stop fit if (error?) function changes go below this tolerance
        self.optim_ftol = 1e-9
        # stop fit if gradient (?) changes go below this tolerance
        self.optim_gtol = 1e-9
        # least squares optimization method (check spicy.optimize.least_squares doc)
        self.optim_method = 'trf'
        # count number of calls to error function
        self._model_call_count = 0
        # record cost function
        self._cost_function = []
        # time
        self._fit_time = 0
        # number of successive fits
        self._fit_count = 0
        # total time
        self._fit_total_time = 0
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

        # --- metabolite list adjustment options (experimental) ---
        # SNR threshold to include or exclude a metabolite
        self.metabolites_auto_adjust_threshold = 1
        # see fit_adjust_metabolite_mode above
        self.metabolites_auto_adjust_mode = fit_adjust_metabolites_mode.NONE
        # when to apply this tweak?
        self.metabolites_auto_adjust_when = fit_adjust_metabolites_when.BEFORE_SECOND_FIT_ONLY

        # --- peak area integration options ---
        # ppm range to look for peak maximum
        self.area_integration_peak_ranges = [0.4]  # ppm

        # freeze
        self.__isfrozen = True

    def copy(self):
        """Copy method."""
        obj = super().copy()

        if(self.sequence is not None):
            if(isinstance(self.sequence, sim.mrs_sequence)):
                obj.sequence = self.sequence.copy()
            else:
                obj.sequence = self.sequence

        if(self.metabolites is not None):
            obj.metabolites = self.metabolites.copy()

        if(self.params_linklock is not None):
            obj.linklock = self.params_linklock.copy()

        # --- fit init/min/max vectors ---
        if(self.params_min is not None):
            obj.params_min = self.params_min.copy()

        if(self.params_max is not None):
            obj.params_max = self.params_max.copy()

        if(self.params_init is not None):
            obj.params_init = self.params_init.copy()

        # --- fit results vectors ---
        if(self.params_fit is not None):
            obj.params_fit = self.params_fit.copy()

        if(self.params_area is not None):
            obj.params_area = self.params_area.copy()

        if(self.params_area_pnorm is not None):
            obj.params_area_pnorm = self.params_area_pnorm.copy()

        # --- optimization options vectors ---
        if(self.optim_results is not None):
            obj.optim_results = self.optim_results.copy()

        if(self.optim_ppm_range is not None):
            obj.optim_ppm_range = self.optim_ppm_range.copy()

        # --- reinit internal paramters for fit ---
        # count number of calls to error function
        obj._model_call_count = 0
        # record cost function
        obj._cost_function = []
        # time
        obj._fit_time = None

        # --- fqn and display options ---
        if(self.fqn_noise_range is not None):
            obj.fqn_noise_range = self.fqn_noise_range.copy()

        if(self.display_range_ppm is not None):
            obj.display_range_ppm = self.display_range_ppm.copy()

        if(self.display_subplots_types is not None):
            obj.display_subplots_types = self.display_subplots_types.copy()

        # --- peak area integration stuff ---
        if(self.metabolites_area_integration is not None):
            obj.metabolites_area_integration = self.metabolites_area_integration.copy()

        if(self.area_integration_peak_ranges is not None):
            obj.area_integration_peak_ranges = self.area_integration_peak_ranges.copy()

        # reinit internal parameters for peak area integration
        self._peak_names = []
        self._peak_ppms = []
        self._peak_nprots = []
        self._peak_areas = []
        self._peak_areas_norm = []

        return(obj)

    def _set_unique_linklock(self):
        """Link the params_min, _max and _init linlock vectors using a reference/pointer."""
        # link-lock vectors should be the same for initial, minimum and maximum parameter sets: link them here
        self.params_init._linklock = self.params_linklock
        self.params_min._linklock = self.params_linklock
        self.params_max._linklock = self.params_linklock

    def _run_area_integration(self):
        """
        Run the area integration pipeline.

        Returns
        -------
        pars : params object
            Estimated relative metabolic concentration using area integration
        pars_norm : params object
            Estimated relative metabolic concentration using area integration, normalized by proton number
        """
        log.info("initializing peak area integration...")

        # clear internal parameters
        self._peak_names = []
        self._peak_ppms = []
        self._peak_nprots = []
        self._peak_areas = []
        self._peak_areas_norm = []

        # check and fix area_integration_peak_ranges
        if(len(self.area_integration_peak_ranges) == 1):
            self.area_integration_peak_ranges = self.area_integration_peak_ranges * len(self.metabolites_area_integration)

        # first find chemical shifts for those peaks
        meta_keys = list(self.meta_bs.keys())
        integration_metagroup_keys = [meta_keys[ind] for ind in self.metabolites_area_integration]
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

        log.info("running peak area integration...")

        # init
        s = self.data.copy()
        # zero-fill for improved resolution
        s = s.correct_zerofill_nd()
        # chemical shift axis
        ppm = s.frequency_axis_ppm()

        # find maximum peak in range and its chemical shift
        sf = s.spectrum()
        sf_real = np.real(sf)

        if(self.display_enable):
            fig = plt.figure(self.display_fig_index)
            fig.clf()
            axs = fig.subplots()
            fig.canvas.set_window_title("run (mrs.fit.fit_pastis)")
            axs.plot(self.data.frequency_axis_ppm(), self.data.spectrum().real, "k-", label="data")

        # now for each peak, integrate the area
        self._peak_areas = []
        self._peak_areas_norm = []
        pars = sim.params(self.meta_bs)
        pars[:] = 0.0
        pars._linklock[:] = 1
        pars_pnorm = sim.params(self.meta_bs)
        pars_pnorm[:] = 0.0
        pars_pnorm._linklock[:] = 1

        log.info("integrating peak for...")
        log.info_line________________________()
        for (this_peak_index, this_peak_range, this_peak_name, this_peak_ppm, this_peak_np) in zip(self.metabolites_area_integration, self.area_integration_peak_ranges, self._peak_names, self._peak_ppms, self._peak_nprots):
            log.info("[%s] theoretically at %.2fppm in theory" % (this_peak_name, this_peak_ppm))

            # first find closest peak in range
            this_peak_search_range = [this_peak_ppm - this_peak_range / 2.0, this_peak_ppm + this_peak_range / 2.0]
            this_peak_ppm, _, this_peak_lw_hz, _, _ = self.data._analyze_peak_1d(this_peak_search_range)
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
            pars[this_peak_index, xxx.p_cm] = this_peak_area
            pars_pnorm[this_peak_index, xxx.p_cm] = this_peak_area / this_peak_np
            # fix linklock vectors
            pars._linklock[this_peak_index, xxx.p_cm] = 0
            pars_pnorm._linklock[this_peak_index, xxx.p_cm] = 0

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

    def _run_linear_combination(self):
        """
        Run the linear combination fit.

        Returns
        -------
        params_fit : params object
            Estimated metabolic parameters using the fitting algorithm
        optim_result : 'OptimizeResult' object returned by scipy.optimize.least_squares
            Numerical optimization parameters, includes the R^2 coefficients [rsq_t and rsq_f] and the FQN coefficient [fqn]
        """
        # ckeck consistency lower<>init<>upper bounds
        log.info("checking consistency with parameter bounds...")
        for ll in range(self.params_init.shape[0]):
            for c in range(self.params_init.shape[1]):
                # min vs max
                if(self.params_min[ll, c] >= self.params_max[ll, c] and self.params_init.linklock[ll, c] != 1):
                    log.error(" One lower bound (%f) is equal/greater than a upper bound (%f) for [%s] at index (%d,%d)!" % (self.params_min[ll, c], self.params_max[ll, c], self.params_init.get_meta_names()[ll], ll, c))
                # min vs init
                if(self.params_min[ll, c] >= self.params_init[ll, c] and self.params_init.linklock[ll, c] != 1):
                    log.error(" One initial value (%f) is equal/lower than the lower bound value (%f) for [%s] at index (%d,%d)!" % (self.params_init[ll, c], self.params_min[ll, c], self.params_init.get_meta_names()[ll], ll, c))
                # max vs init
                if(self.params_max[ll, c] <= self.params_init[ll, c] and self.params_init.linklock[ll, c] != 1):
                    log.error(" One initial value (%f) is equal/greater than the upper bound value (%f) for [%s] at index (%d,%d)!" % (self.params_init[ll, c], self.params_max[ll, c], self.params_init.get_meta_names()[ll], ll, c))

        # checking if LL is not broken
        if(not self.params_init.check()):
            log.error("the link-lock vector looks broken!")

        # checking if we are dealing with a reference scan fit (water only)
        where_LL = np.where(self.params_init.linklock[:, 0] == 0)
        self._water_only = (len(where_LL[0]) == 1 and where_LL[0][0] == xxx.m_Water)
        if(self._water_only):
            log.info("mmmh, looks like you are fitting a water reference scan!")

        log.info("running fit...")
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
            optim_result = optimize.least_squares(self._minimizeThis, params_free, bounds=(params_min_free, params_max_free), jac=self._jac, xtol=self.optim_xtol, gtol=self.optim_gtol, ftol=self.optim_ftol, method=self.optim_method)
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
        params_CRBs_abs = self.get_CRBs(params_fit)
        params_fit._errors = params_CRBs_abs
        # final corr mat
        corr_mat, _ = self.get_corr_mat(params_fit)
        params_fit._corr_mat = corr_mat

        # add fit criteria to optim_result
        rsq_t, rsq_f = self.get_Rsq(params_fit)
        optim_result.rsq_t = rsq_t
        optim_result.rsq_f = rsq_f
        optim_result.fqn = self.get_FQN(params_fit)

        return(params_fit, optim_result)

    def run(self):
        """Run the peak area integration followed by the linear combination fit."""
        log.info("initializing...")

        # check that we have everything we need
        if(self.data is None):
            log.error("canceling fit, no data was given to this fit object!")

        if(self.metabolites is None):
            log.error("canceling fit, the list of metabolites to fit was not specified!")

        # sequence
        if(self.sequence is None):
            self.sequence = self.data.sequence
        else:
            self.sequence = self.sequence(self.data.sequence.te,
                                          self.data.sequence.tr,
                                          self.data.sequence.na,
                                          self.data.sequence.ds,
                                          self.data.sequence.nuclei,
                                          self.data.sequence.npts,
                                          self.data.sequence.voxel_size,
                                          self.data.sequence.fs,
                                          self.data.sequence.f0)

        # do we need to adapt our ppm range?
        if(self.optim_ppm_range != self.meta_bs.ppm_range):
            # we need to reinitialize the metabolite basis set
            self.meta_bs.ppm_range = self.optim_ppm_range
            # force reinit of sequence
            self.sequence._ready = False

        # initialize sequence if needed
        if(not self.sequence.ready):
            self.sequence.initialize(self.meta_bs)

        # running metabolite list optimization now if required and if we already have fit results
        if((self.metabolites_auto_adjust_mode is not fit_adjust_metabolites_mode.NONE) and
           (self.metabolites_auto_adjust_when is not fit_adjust_metabolites_when.NEVER) and
           (self.params_fit is not None)):
            if(((self.metabolites_auto_adjust_when == fit_adjust_metabolites_when.BEFORE_SECOND_FIT_ONLY) and (self._fit_count == 1)) or
               ((self.metabolites_auto_adjust_when == fit_adjust_metabolites_when.BEFORE_EACH_FIT))):

                # test metabolite basis set
                new_metabolites_list, metabolites_excluded_list = self.get_fittable_metabolites(self.metabolites_auto_adjust_threshold, self.metabolites_auto_adjust_mode)

                if(len(metabolites_excluded_list) > 0):
                    self.metabolites = new_metabolites_list
                    # fix params according to new metabolite basis (hoping it does not break anything)
                    # null excluded metabolites
                    self.params_init[metabolites_excluded_list, xxx.p_cm] = 0.0
                    # and lock them
                    self.params_linklock[metabolites_excluded_list, :] = 1

        # reset linklock
        self._set_unique_linklock()

        # run
        pars_area, pars_area_pnorm = self._run_area_integration()
        params_fit, optim_result = self._run_linear_combination()

        # number of successive fits
        self._fit_count += 1
        # total time
        self._fit_total_time += self._fit_time

        # store
        log.info("saving fit results...")
        self.params_fit = params_fit.copy()
        self.params_area = pars_area.copy()
        self.params_area_pnorm = pars_area_pnorm.copy()
        self.optim_results = optim_result

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
        j_full = self.sequence._jac(params)

        # linlock the jacobian vector
        # for each master parameter, sum the slave derivates
        LL_list = np.unique(params.linklock)
        LL_list = LL_list[LL_list >= +2]
        for this_LL in LL_list:
            j_full[params.linklock == -this_LL, :] += np.sum(j_full[params.linklock == +this_LL, :], 0)

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
        mod = self.sequence._model(params)
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
            fig.canvas.set_window_title("_minimizeThis (mrs.fit.fit_pastis)")

            # display real-time bargraphs
            # real-time CRBs
            if(self.display_CRBs):
                params._errors = self.get_CRBs(params)

            # let's display the subplots as required by self.display_subplots_types
            axs_list = [axs[0, 2], axs[0, 3], axs[1, 2], axs[1, 3]]
            for ax, plt_type in zip(axs_list, self.display_subplots_types):
                if(plt_type == fit_plot_type.BARGRAPH_CM):
                    self.disp_bargraph(ax, params, xxx.p_cm)
                elif(plt_type == fit_plot_type.BARGRAPH_DD):
                    self.disp_bargraph(ax, params, xxx.p_dd)
                elif(plt_type == fit_plot_type.BARGRAPH_DF):
                    self.disp_bargraph(ax, params, xxx.p_df)
                elif(plt_type == fit_plot_type.BARGRAPH_DP):
                    self.disp_bargraph(ax, params, xxx.p_dp)
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
            disp_fit(ax, self.data, params, self.sequence, True, True, None, self._water_only, self.display_range_ppm)
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

    def get_fittable_metabolites(self, snr_threshold=1, mode=fit_adjust_metabolites_mode.AVERAGE):
        """
        Adjust the list of metabolites to fit. This method is experimental. It runs an algorithm to try and test which metabolites are actually fittable. How? Based on their common concentration and simulated time-domain SNR. If the metabolite SNR is above the threshold (1, by default, meaning above noise level in short), the metabolite is kept in the fit list. The final metabolite basis set is printed on screen. This method should be run after the run() method in order to get some first fit results.

        Parameters
        ----------
        snr_threshold : float
            SNR threshold in order to keep the metabolite in the fit metabolite basis set. This is a time-domain SNR calculated using the maximum point of the real part of FID divided by the STD of noise at the end of the FID.
        mode : fit_adjust_metabolite_mode
            This describes how optimistic is the algorithm. When running on the OPTIMISTIC mode, the algorithm will assume a low concentration for the reference metabolite (Cr_CH3 for example) based on the concentration ranges in the xls metabolite basis set. It will also keep in the basis set metabolites giving a SNR above the threshold even with the lowest concentration. AVERAGE mode uses average concentrations while PESSIMISTIC uses maximum concentrations available in the xls metabolite basis set.

        Returns
        -------
        new_metabolites_list : list
            New list of metabolite indexes to fit
        metabolites_excluded_list : list
            List of metabolite indexes to exclude from fit
        """
        # first get original noise level
        noise_std = self.data.noise_level

        # get reference metabolite concentration
        # correct for T1/T2 effects
        params_fit_ref = self.params_fit.correct_T2s(self.data.te).correct_T1s(self.data.tr)
        # go and get concentration for this metabolites according to literature
        if(mode == fit_adjust_metabolites_mode.AVERAGE):
            cm_lit_ref = (self.meta_bs.get_literature_min() + self.meta_bs.get_literature_max()) / 2.0
            cm_lit_test = (self.meta_bs.get_literature_min() + self.meta_bs.get_literature_max()) / 2.0
        elif(mode == fit_adjust_metabolites_mode.OPTIMISTIC):
            cm_lit_ref = self.meta_bs.get_literature_min()
            cm_lit_test = self.meta_bs.get_literature_max()
        if(mode == fit_adjust_metabolites_mode.PESSIMISTIC):
            cm_lit_ref = self.meta_bs.get_literature_max()
            cm_lit_test = self.meta_bs.get_literature_min()

        # get our reference metabolite from metabolite basis set
        m_ref = xxx.m_Ref_MB

        # for each metabolite we attempt to fit here, we test its SNR
        snr_metabolites_list = []
        metabolites_kept_list = []
        metabolites_excluded_list = []

        # metabolites
        meta_names = params_fit_ref.get_meta_names()
        metabolites_to_test_list = [m for m in self.metabolites if(m in xxx.m_All_MBs)]
        macromolecules_to_keep_list = [m for m in self.metabolites if(m in xxx.m_All_MMs)]

        # init display
        if(self.display_enable):
            n_subplots = len(metabolites_to_test_list)
            n_subplots_side = int(np.ceil(np.sqrt(n_subplots)))
            fig = plt.figure(self.display_fig_index)
            fig.clf()
            axs = fig.subplots(n_subplots_side, n_subplots_side, sharex='all')
            fig.canvas.set_window_title("adjust_metabolites (mrs.fit.fit_pastis)")

        for plot_index, m in enumerate(metabolites_to_test_list):
            # simulate its signal with the noise level estimated on the data
            p = sim.params(self.meta_bs)
            p[:, xxx.p_cm] = 0.0
            # concentration parameter according to fit results and known concentrations
            p[m, xxx.p_cm] = params_fit_ref[m_ref, xxx.p_cm] * cm_lit_test[m] / cm_lit_ref[m_ref]
            # take damping/T2* parameter from fit results of ref metabolite
            p[:, xxx.p_dd] = params_fit_ref[m_ref, xxx.p_dd]
            # go simulation without noise
            s = self.sequence.simulate_signal(p, 0.0)

            # noise only
            p[:, xxx.p_cm] = 0.0
            # go simulation
            s_noise = self.sequence.simulate_signal(p, noise_std)

            # estimate time SNR
            this_snr_s = np.max(np.real(s))
            this_snr_n = noise_std
            this_snr = this_snr_s / this_snr_n
            snr_metabolites_list.append(this_snr)

            # check snr and store
            if(this_snr > snr_threshold):
                metabolites_kept_list.append(m)
            else:
                metabolites_excluded_list.append(m)

            # display
            if(self.display_enable):
                axs.flat[plot_index].plot(s_noise.time_axis() * 1000.0, np.real(s_noise), 'k')
                axs.flat[plot_index].axhline(this_snr_n, color='w', linestyle='--')
                if(this_snr > snr_threshold):
                    axs.flat[plot_index].plot(s.time_axis() * 1000.0, np.real(s), 'g')
                else:
                    axs.flat[plot_index].plot(s.time_axis() * 1000.0, np.real(s), 'r')
                axs.flat[plot_index].set(xlabel='FID time (ms)', title=meta_names[m])

        # display
        if(self.display_enable):
            for ax in fig.get_axes():
                ax.label_outer()

            fig.subplots_adjust(left=0.02, bottom=0.1, right=0.98, top=0.95, wspace=0.2, hspace=None)
            fig.show()
            plt.pause(0.5)

        # display results
        log.info("Initial basis set (n=%d): %s" % (len(self.metabolites), " ".join([meta_names[im] for im in self.metabolites])))
        log.info("Metabolite tested (n=%d): %s" % (len(metabolites_to_test_list), " ".join([meta_names[im] for im in metabolites_to_test_list])))
        log.info("Metabolite included (n=%d): %s" % (len(metabolites_kept_list), " ".join([meta_names[im] for im in metabolites_kept_list])))
        log.info("Metabolite excluded (n=%d): %s" % (len(metabolites_excluded_list), " ".join([meta_names[im] for im in metabolites_excluded_list])))

        # merge back with macromolecules
        new_metabolites_list = list(np.sort(list(set(metabolites_kept_list + macromolecules_to_keep_list))))
        log.info("Final basis set (n=%d): %s" % (len(new_metabolites_list), " ".join([meta_names[im] for im in new_metabolites_list])))

        # and return
        return(new_metabolites_list, metabolites_excluded_list)

    def get_hash(self):
        """
        Generate a hash of this strategy.

        Returns
        -------
        h : string
            Hash code of this strategy. Usefull later when dealing with databases...
        """
        bytes_to_hash = np.array(self.metabolites).tobytes() + self.params_linklock.tobytes() + str(self.sequence).encode()

        h = hashlib.md5(bytes_to_hash)

        return(h.hexdigest())

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
        """
        # calling jacobian model
        j_full = self.sequence._jac(params)

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

        # put back to full parameter form
        params_CRBs_abs = self.params_init.copy()
        params_CRBs_abs = params_CRBs_abs.toFullParams(CRBs_abs)

        return(params_CRBs_abs)

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
        mod = self.sequence._model(params)

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
        mod = self.sequence._model(params)
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

    def get_corr_mat(self, params, mIndex_list=None, pIndex_list=[xxx.p_cm]):
        """
        Return and display the parameter correlation matrix.

        Parameters
        ----------
        params : params object
            Array of simulation parameters
        mIndex_list : list of int
            Metabolite indexes to display in bargraph
        pIndex_list : list of int
            Parameter indexes to display in bargraph

        Returns
        -------
        corr_mat : numpy array
            Correlation matrix
        corr_mat_lbls : numpy array
            Correlation matrix labels for plot output
        """
        # calling jacobian model
        j_full = self.sequence._jac(params)

        # probably do not want to look at the full matrix but only some metabolites/parameters
        mat_linklock_mIndex = np.full(params.linklock.shape, False, dtype=bool)
        mat_linklock_mIndex[mIndex_list, :] = True
        mat_linklock_pIndex = np.full(params.linklock.shape, False, dtype=bool)
        mat_linklock_pIndex[:, pIndex_list] = True
        mat_linklock_mask = np.logical_and(np.logical_and(mat_linklock_mIndex, mat_linklock_pIndex), (params.linklock <= 0))

        # build labels vector the same way
        corr_mat_lbls = []
        meta_names = params.get_meta_names()
        for i, m in enumerate(meta_names):
            for j, p in enumerate(PARS_SHORT_NAMES):
                if(mat_linklock_mask[i, j]):
                    corr_mat_lbls.append(m + "|" + p)

        if(np.all(mat_linklock_mask == False)):
            log.error("empty fit correlation matrix, please adjust the mIndex_list/pIndex_list parameters!")

        # reducing it to metabolites, parameters of interest + free params
        j_reduced = j_full[mat_linklock_mask, :]
        j_reduced = np.real(j_reduced)

        # Fisher (number of free metabolites*number of params, time points)
        F = np.dot(j_reduced, np.transpose(j_reduced))

        # covariance
        covar = np.linalg.inv(np.real(F))

        # parameter correlation
        corr_mat = np.zeros(covar.shape)
        for i in range(corr_mat.shape[0]):
            for j in range(corr_mat.shape[1]):
                corr_mat[i, j] = covar[i, j] / (np.sqrt(covar[i, i]) * np.sqrt(covar[j, j]))

        return(corr_mat, corr_mat_lbls)

    def disp_corr_mat(self, ax, params, disp_color_bar=False):
        """
        Return and display the parameter correlation matrix.

        Parameters
        ----------
        ax : matplotlib axis
            Axis to use for plotting
        params : params object
            Array of simulation parameters
        disp_color_bar : boolean
            Add a color (True) or not (False)
        """
        # calling jacobian model
        corr_mat, corr_mat_lbls = self.get_corr_mat(params)

        ax.cla()
        im = ax.imshow(corr_mat, clim=(-1, +1))
        ax.set_xticks(np.arange(len(corr_mat_lbls)))
        ax.set_yticks(np.arange(len(corr_mat_lbls)))
        ax.set_xticklabels(corr_mat_lbls, rotation=90, fontsize=6)
        ax.set_yticklabels(corr_mat_lbls, fontsize=6)
        if(disp_color_bar):
            plt.colorbar(im)

    def disp_bargraph(self, ax, params, pIndex=xxx.p_cm):
        """
        Plot a bargraph of concentrations.

        Parameters
        ----------
        ax : matplotlib axis
            Axis to use for plotting
        params : sim.params
            Parameters to barplot
        pIndex : int
            Parameter index to display in bargraph
        """
        # LLcolor_names=[name for name in mcd.XKCD_COLORS]
        LLcolor_names = ['xkcd:blue', 'xkcd:red', 'xkcd:violet', 'xkcd:green',
                         'xkcd:goldenrod', 'xkcd:crimson', 'xkcd:salmon', 'xkcd:wheat', 'xkcd:lightblue']
        # fixed bar width, maybe an issue one day
        width = 0.3
        # label
        lbl = PARS_LONG_NAMES[pIndex]

        # build logical mask
        meta_mask = np.full(params.shape, True, dtype=bool)
        # filter out the fixed LLs
        meta_mask[params.linklock[:, xxx.p_cm] == 1, :] = False
        # filter out water?
        if(not self._water_only):
            meta_mask[xxx.m_Water, :] = False

        # parameter filter
        meta_mask_p = meta_mask.copy()
        meta_mask_p[:] = False
        meta_mask_p[:, pIndex] = True
        meta_mask = np.logical_and(meta_mask, meta_mask_p)

        # check that we still have something to display
        if(np.all(meta_mask == False)):
            log.error("nothing to display! :(")

        # prepare data
        params = params[meta_mask]
        params_min = self.params_min[meta_mask]
        params_max = self.params_max[meta_mask]
        params_range = params_max - params_min  # we want to plot the ranges
        params_LL = params.linklock[meta_mask]
        params_std = params.errors[meta_mask]
        meta_names = params.get_meta_names()
        meta_names = [meta_names[i] for i in range(len(meta_names)) if meta_mask[i, pIndex]]

        # prepare bars
        n = np.sum(meta_mask[:])
        pos_bars = np.arange(n)

        # draw bars
        ax.bar(pos_bars, params_range, width * 2, bottom=params_min, color="grey")
        bar_params = ax.bar(pos_bars, params, width, yerr=params_std, color="red")

        # set color according to linklock
        for this_bar, this_LL in zip(bar_params, params_LL):
            this_LLcolor_name = LLcolor_names[np.mod(int(np.abs(this_LL)), len(LLcolor_names))]
            this_LLcolor = mcd.XKCD_COLORS[this_LLcolor_name].upper()
            this_bar.set_color(this_LLcolor)

        # ylim
        if(pIndex in [xxx.p_cm, xxx.p_dd]):
            ax.set_ylim([0.0, params_max.max()])
        elif(pIndex in [xxx.p_df, xxx.p_dp]):
            minmax = max(np.max(np.abs(params_min)), np.max(np.abs(params_max)))
            ax.set_ylim([-minmax, +minmax])

        ax.set_ylabel(lbl)
        ax.set_xticks(pos_bars)
        ax.set_xticklabels(meta_names, rotation=90)
        ax.grid('on')

    def to_dataframe(self, prefix_str="fit_"):
        """
        Convert the object's attributes to dataframe. Can include the object itself.

        Parameters
        ----------
        prefix_str : string
            Prefix string to add to column names

        Returns
        -------
        df : Dataframe
            With all the datasets and parameters from this fit object
        """
        log.debug("converting to dataframe...")

        # use pandas json_normalize function to flatten all nested stuff (magic!)
        # this takes care of the job attribute too
        df_attr = pd.json_normalize(vars(self), sep='_')

        # and need a specific call for the MRSData2 data
        df_data = self.data.to_dataframe(True, "data_")
        del df_attr["data"]

        # params objects
        df_params_min = self.params_min.to_dataframe("params_min_")
        df_params_max = self.params_max.to_dataframe("params_max_")
        df_params_init = self.params_init.to_dataframe("params_init_")

        # keep objects
        df_attr = df_attr.rename(columns = {'params_min': 'params_min_obj',
                                            'params_max': 'params_max_obj',
                                            'params_init': 'params_init_obj'})

        # if any results, do them too
        if((self.params_fit is not None) and (self.params_area is not None) and (self.params_area_pnorm is not None)):
            df_params_fit = self.params_fit.to_dataframe("params_fit_")
            df_params_area = self.params_area.to_dataframe("params_area_")
            df_params_area_pnorm = self.params_area_pnorm.to_dataframe("params_area_pnorm_")

            # keep objects
            df_attr = df_attr.rename(columns = {'params_fit': 'params_fit_obj',
                                            'params_area': 'params_area_obj',
                                            'params_area_pnorm_init': 'params_area_pnorm_obj'})

            # append columns
            df = pd.concat([df_attr.reset_index(),
                            df_data.reset_index(),
                            df_params_min.reset_index(),
                            df_params_max.reset_index(),
                            df_params_init.reset_index(),
                            df_params_fit.reset_index(),
                            df_params_area.reset_index(),
                            df_params_area_pnorm.reset_index()], axis=1)
        else:
            # append columns
            df = pd.concat([df_attr.reset_index(),
                            df_data.reset_index(),
                            df_params_min.reset_index(),
                            df_params_max.reset_index(),
                            df_params_init.reset_index()], axis=1)


        # add fit hash
        df["fit_hash"] = self.get_hash()

        # remove index column
        df = df.reset_index()
        del df["index"]

        # remove this weird column
        del df["level_0"]

        # add prefix
        df = df.add_prefix(prefix_str)

        return(df)


def disp_fit(ax, data, params, seq, LL_exluding=True, LL_merging=False, mIndex_list=None, water_only=False, display_range=[1, 5], zerofilling_final_npts=16384):
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
    zerofilling_final_npts : int
        Final npts of points after zero-filling data to smooth display; if (None), no zero-filling
    """
    # dummy full>free>full conversion to apply master/slave rules
    params_free = params.toFreeParams()
    params_full = params.toFullParams(params_free)

    # the data
    if(zerofilling_final_npts is not None):
        data = data.correct_zerofill_nd(zerofilling_final_npts)
    ax.plot(data.frequency_axis_ppm(), data.spectrum().real, 'k-', linewidth=0.5, label='data')
    # the full model
    mod = seq._model(params_full)
    if(zerofilling_final_npts is not None):
        mod = mod.correct_zerofill_nd(zerofilling_final_npts)
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

    ax.text(display_range[1] - 0.5, -ygap * (koffset + 2) + ygap / 5, 'Lip baseline')

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
