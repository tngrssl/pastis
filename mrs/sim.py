#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Three classes' definition in here.

    * a params class which stores and manipulate the parameters of our MRS fitting/simulation model
    * a metabolite class which stores and can compute a MRS modeled signal for a single metabolite, based on the pyGAMMA library (for python 3!) using a specific MR sequence described by pulse flip angles and delays
    * a metabolite_group class which contains several metabolites
    * a metabolite_basis_set class which contains a whole database of metabolites chemical shifts, J-couplings, nucleis and computed signals. Usefull to simulate all kind of MRS data for various metabolites, concentrations acquired with various sequences...

@author: Tangi Roussel
"""

try:
    import pygamma as pg
    GAMMA_LIB_LOADED = True
except ImportError:
    GAMMA_LIB_LOADED = False

# GAMMA_LIB_LOADED forced to False for debug
# GAMMA_LIB_LOADED = False

import suspect
import numpy as np
import math as ma
import matplotlib.pylab as plt
import pickle
import warnings
import pathlib
from xlrd import open_workbook
from termcolor import cprint
from enum import Enum
from mrs import reco
from mrs import aliases as xxx
from mrs import log
from mrs import paths as default_paths
import copy as copy

import pdb


class sequence_exc_type(Enum):
    """The enum sequence_exc_type describes the type of excitation scheme of the sequence. Can be usefull when comparing sequences."""

    PULSE_ACQUIRE = 1
    STIMULATED_ECHO = 2
    SPIN_ECHO = 3


class gating_signal_source(Enum):
    """The enum gating_signal_source describes the type of gating used during the acquisition."""

    NO_GATING = 0
    CARDIAC_ECG = 2
    CARDIAC_GATING = 4
    RESP_GATING = 16


class params(np.ndarray):
    """A class that stores the parameters used to modelize a MR spectrum during simulation or fit."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, meta_bs):
        """
        Initialize a params object.

        Parameters
        ----------
        meta_bs: metabolite_basis_set
            A metabolite_basis_set object to which this params object is linked to
        """
        super().__init__()
        # those parameters are related to a metabolite database
        self._meta_bs = meta_bs
        # the link-lock vector used to control the model
        self._linklock = np.zeros(self.shape)
        # the error vector
        self._errors = np.zeros(self.shape)
        # the corr vector
        self._corr_mat = None

        # freeze
        self.__isfrozen = True

    def __new__(cls, meta_bs):
        """
        Construct a params object that inherits of numpy array's class. This class is used to deal with metabolite parameters.

        Parameters
        ----------
        meta_bs: metabolite_basis_set
            A metabolite_basis_set object to which this params object is linked to

        Returns
        -------
        obj : params numpy array [n,4]
            Resulting constructed params object
        """
        obj = super(params, cls).__new__(cls, [len(meta_bs), 4])
        obj[:] = 0.0
        return(obj)

    def __array_finalize__(self, obj):
        """
        Overload of special numpy array method called when playing around with stuff relative to object copy etc...

        Parameters
        ----------
        obj : params numpy array [n,4]
        """
        # to begin, I followed online help and wrote:
        #self._meta_bs = getattr(obj, 'meta_bs', None)
        #self._linklock = getattr(obj, 'linklock', None)
        #self._errors = getattr(obj, 'errors', None)
        #self._corr_mat = getattr(obj, 'corr_mat', None)

        # but that only works for some simple attribute types
        # if the attributes are nd arrays, only the pointers will be copied...
        # leading to terrible bugs

        # for now, I could only find this ugly way:

        self._meta_bs = getattr(obj, 'meta_bs', None)
        if(self.meta_bs is not None):
            self._meta_bs = obj.meta_bs.copy()

        self._linklock = getattr(obj, 'linklock', None)
        if(self.linklock is not None):
            self._linklock = obj.linklock.copy()

        self._errors = getattr(obj, 'errors', None)
        if(self.errors is not None):
            self._errors = obj.errors.copy()

        self._corr_mat = getattr(obj, 'corr_mat', None)
        if(self.corr_mat is not None):
            self._corr_mat = obj.corr_mat.copy()

    @property
    def meta_bs(self):
        """Property method for meta_bs."""
        return(self._meta_bs)

    @property
    def linklock(self):
        """Property method for linklock."""
        return(self._linklock)

    @property
    def errors(self):
        """Property method for errors."""
        return(self._errors)

    @property
    def corr_mat(self):
        """Property method for corr_mat."""
        return(self._corr_mat)

    def get_meta_names(self):
        """
        Return list of metabolite names controlled by this params object.

        Returns
        -------
        list(self._meta_bs.keys()) : list
            List of metabolite names
        """
        return(list(self._meta_bs.keys()))

    def get_errors_prct(self):
        """
        Return relative errors in percentage.

        Returns
        -------
        errors_prct : numpy array
            Errors in percentage (%)
        """
        # calculate percentages
        errors_prct = self.errors / self[:] * 100.0

        return(errors_prct)

    def check(self):
        """
        Check if the linklock vector is consistent and read to be used.

        Returns
        -------
        all_right : boolean
            True if eveything is ok, False if linklock vector if broken
        """
        # by default, eveything is ok
        all_right = True

        LL_list = np.unique(self.linklock)
        LL_list = LL_list[LL_list >= +2]
        for this_LL in LL_list:
            # count number of master
            tmp = (self.linklock == -this_LL)
            n_masters = np.sum(tmp[:])

            # count number of slaves
            tmp = (self.linklock == +this_LL)
            n_slaves = np.sum(tmp[:])

            # master ratio values should be all equal to one
            tmp = (self.linklock == -this_LL)
            if(np.any(tmp[:]) != 1.0):
                all_right = False

            # now, we should have one single master
            if(n_masters != 1):
                all_right = False

            # and several slaves
            if(n_slaves < 1):
                all_right = False

        return(all_right)

    def toFreeParams(self):
        """
        Return free parameters, in other words, parameters for which LL=0 (free) or LL<0 (free and master).

        Returns
        -------
        self[self.linklock<=0] : 1D numpy array of free parameters
            Array of free parameters
        """
        return(self[self.linklock <= 0])

    def toFullParams(self, pFree):
        """
        Convert an array of free paramters to a full array of parameters based on the link-lock array and the initial parameter values.

        Parameters
        ----------
        pFree : 1D numpy array of free parameters
            Array of free parameters

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        # init: we are working on a copy
        p = self.copy()

        # link-lock stuff
        p[p.linklock <= 0] = pFree

        # some of them are masters, copying values to slaves
        # find unique LL values
        LL_list = np.unique(p.linklock)
        LL_list = LL_list[LL_list >= +2]
        for this_LL in LL_list:
            p[p.linklock == +this_LL] = p[p.linklock == -this_LL]

        return(p)

    def set_T2_weighting(self, te):
        """
        Recalculate metabolite amplitude parameters by applying a T2 relaxation, usefull for simulations.

        Parameters
        ----------
        te : float
            Echo time in (ms)

        Returns
        -------
        self.copy() : params object
            Full array of parameters corrected in T2
        """
        # init: we are working on a copy
        p = self.copy()

        # browse though the database and find T2s
        params_T2s = []
        for this_metagroup_key, this_metagroup_entry in self._meta_bs.items():
            params_T2s.append(this_metagroup_entry["T2"])

        # convert to np
        params_T2s = np.array(params_T2s)

        # apply T2w for given TE
        multiplication_factor = np.exp(-te / params_T2s)
        p[:, xxx.p_cm] = p[:, xxx.p_cm] * multiplication_factor
        # to errors too
        p.errors[:, xxx.p_cm] = p.errors[:, xxx.p_cm] * multiplication_factor
        return(p)

    def correct_T2s(self, te):
        """
        Correct the concentration values of a parameter array depending on the TE and the common values of T2s for each metabolite.

        Parameters
        ----------
        te : float
            Echo time in (ms)

        Returns
        -------
        self.copy() : params object
            Full array of parameters corrected in T2
        """
        # init: we are working on a copy
        p = self.copy()

        # browse though the database and find T2s
        params_T2s = []
        for this_metagroup_key, this_metagroup_entry in self._meta_bs.items():
            params_T2s.append(this_metagroup_entry["T2"])

        # convert to np
        params_T2s = np.array(params_T2s)

        # finding real concentration values at TE=0ms
        multiplication_factor = 1 / np.exp(-te / params_T2s)
        p[:, xxx.p_cm] = p[:, xxx.p_cm] * multiplication_factor
        # to errors too
        p.errors[:, xxx.p_cm] = p.errors[:, xxx.p_cm] * multiplication_factor
        return(p)

    def correct_T1s(self, tr):
        """
        Correct the concentration values of a parameter array depending on the TR and the common values of T1s for each metabolite.

        Parameters
        ----------
        tr : float
            Repetition time in (ms)

        Returns
        -------
        self.copy() : params object
            Full array of parameters corrected in T1
        """
        # init: we are working on a copy
        p = self.copy()

        # browse though the database and find T1s
        params_T1s = []
        for this_metagroup_key, this_metagroup_entry in self._meta_bs.items():
            params_T1s.append(this_metagroup_entry["T1"])

        # convert to np
        params_T1s = np.array(params_T1s)

        # finding real concentration values at TE=0ms
        multiplication_factor = 1 / (1 - np.exp(-tr / params_T1s))
        p[:, xxx.p_cm] = p[:, xxx.p_cm] * multiplication_factor
        # to errors too
        p.errors[:, xxx.p_cm] = p.errors[:, xxx.p_cm] * multiplication_factor
        return(p)

    def get_relative_to_meta(self, mIndex=xxx.m_Water, m_concentration_mmolkg=55000.0, params_ref=None):
        """
        Calculate the metabolic concentration values relative to a metabolite. The metabolite relative concentration can be taken from the current params vector or another params vector (params_ref), assuming the absolute metabolite concentration value. Usefull to get concentrations relative to water (called absolute concentrations).

        Parameters
        ----------
        mIndex : int
            Index of metabolite of reference
        water_concentration : float
            Assumed water concentration used to calculate absolute concentration estimates (mmol/kg)
        params_ref : params object
            Array of parameters used to get the concentration estimate of the reference metabolite

        Returns
        -------
        self.copy() : params object
            Full array of parameters
        """
        # init: check ref params
        if(params_ref is None):
            params_ref = self.copy()

        # init: we are working on copies
        p1 = self.copy()
        p2 = self.copy()

        # do the maths
        multiplication_factor = m_concentration_mmolkg / params_ref[mIndex, xxx.p_cm]
        p2[:, xxx.p_cm] = p1[:, xxx.p_cm] * multiplication_factor

        # deal with errors too
        # see paper http://dx.doi.org/10.1007/s10334-005-0018-7
        # inspired from ratio, considering zero correlation between metabolite of interest and ref:
        # relCRB(1/ref) = sqrt(relCRB1^2 + relCRB2^2)
        # get numerator rel CRBs for cm
        relCRBs_num = p1.get_errors_prct()[:, xxx.p_cm]
        # get denominator rel CRBs for cm
        relCRBs_den = p2.get_errors_prct()[mIndex, xxx.p_cm]
        # calculate the final relCRB
        rel_CRBs_ratio = np.sqrt( relCRBs_num**2 + relCRBs_den**2 )
        # back to absCRB
        abs_CRBs_ratio = p2[:, xxx.p_cm] * rel_CRBs_ratio / 100.0

        p2._errors[:, xxx.p_cm] = abs_CRBs_ratio
        return(p2)

    def get_ratios(self, mIndex):
        """
        Calculate the metabolite ratios.

        Parameters
        ----------
        mIndex : int
            Index of metabolite of reference used as the denominator

        Returns
        -------
        self.copy() : params object
            Full array of parameters
        """
        # init: we are working on a copy
        p1 = self.copy()
        p2 = self.copy()

        # ratio
        p2[:, xxx.p_cm] = p1[:, xxx.p_cm] / p1[mIndex, xxx.p_cm]

        # deal with errors too
        # if correlation matrix available
        if(self.corr_mat is not None):
            # see paper http://dx.doi.org/10.1007/s10334-005-0018-7
            # relCRB(1/2) = sqrt(relCRB1^2 + relCRB2^2 - 2*corr1_2*relCRB1*relCRB2)
            # get numerator rel CRBs for cm
            relCRBs_num = p1.get_errors_prct()[:, xxx.p_cm]
            # get denominator rel CRBs for cm
            relCRBs_den = p1.get_errors_prct()[mIndex, xxx.p_cm]

            # get corr coeff between num and den cm
            # first, get free param corr mat
            free_param_corr_mat = self.corr_mat
            # convert mIndex to free params index
            p3 = self.copy()
            # replace cm by index
            p3[:, xxx.p_cm] = np.arange(0, p3.shape[0], 1)
            # lock all other pars
            p3.linklock[:, xxx.p_dd] = 1
            p3.linklock[:, xxx.p_df] = 1
            p3.linklock[:, xxx.p_dp] = 1
            # convert ot free params
            p3_free = p3.toFreeParams()
            # here we should get a list of indexes
            # find where is mIndex and we did it
            ind_free_pars_mIndex = np.where(p3_free == mIndex)[0][0]
            # extract corr vector
            free_param_corr_vec = free_param_corr_mat[:, ind_free_pars_mIndex]
            # convert it to full params
            p3 = self.copy()
            # lock all other pars
            p3.linklock[:, xxx.p_dd] = 1
            p3.linklock[:, xxx.p_df] = 1
            p3.linklock[:, xxx.p_dp] = 1
            # replace cm by cor coeffs
            p3 = p3.toFullParams(free_param_corr_vec)
            # get full param corr vector
            full_param_corr_vec = p3[:, xxx.p_cm]
            # make it absolute
            rc_num_den_abs = np.abs(full_param_corr_vec)

            # calculate the final relCRB
            rel_CRBs_ratio = np.sqrt( relCRBs_num**2 + relCRBs_den**2 - 2 * rc_num_den_abs * relCRBs_num * relCRBs_den )
            # back to absCRB
            abs_CRBs_ratio = p2[:, xxx.p_cm] * rel_CRBs_ratio / 100.0
            p2._errors[:, xxx.p_cm] = abs_CRBs_ratio

        return(p2)

    def _set_default_min(self):
        """
        Initialize the params object to the default minimum values, no macromolecules.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        # all to zero
        self[:] = 0.0
        # metabolites min values
        self[xxx.m_All_MBs, xxx.p_cm] = 0.0
        self[xxx.m_All_MBs, xxx.p_dd] = 5.0
        self[xxx.m_All_MBs, xxx.p_df] = -10.0
        self[xxx.m_All_MBs, xxx.p_dp] = -0.1

        # link all to the NAA singlet
        self.linklock[:] = np.tile([0, 2, 3, 4], (xxx.n_All, 1))
        self.linklock[xxx.m_Ref_MB, :] = [0, -2, -3, -4]
        self.linklock[xxx.m_Water, :] = 0

        # no MMs
        self.linklock[xxx.m_All_MMs, :] = 1

        return(self.copy())

    def _set_default_max(self):
        """
        Initialize the params object to the default maximum values, no macromolecules.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        # all to zero
        self[:] = 0.0
        # metabolites min values
        self[xxx.m_All_MBs, xxx.p_cm] = 50.0
        self[xxx.m_All_MBs, xxx.p_dd] = 100.0
        self[xxx.m_All_MBs, xxx.p_df] = +10.0
        self[xxx.m_All_MBs, xxx.p_dp] = +0.1

        # link all to the NAA singlet
        self.linklock[:] = np.tile([0, 2, 3, 4], (xxx.n_All, 1))
        self.linklock[xxx.m_Ref_MB, :] = [0, -2, -3, -4]
        self.linklock[xxx.m_Water, :] = 0

        # no MMs
        self.linklock[xxx.m_All_MMs, :] = 1

        return(self.copy())

    def set_default_water_min(self):
        """
        Initialize the params object to the default minimum values for a water only fit.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        self._set_default_min()
        # all concentrations to zero
        self[:, xxx.p_cm] = 0.0

        # water min values
        self[xxx.m_All_MBs, xxx.p_df] = -20.0

        # lock everything except water
        self.linklock[:] = 1
        self.linklock[xxx.m_Water, :] = 0
        return(self.copy())

    def set_default_water_max(self):
        """
        Initialize the params object to the default maximum values for a water only fit.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        self._set_default_max()
        # all concentrations to zero
        self[:, xxx.p_cm] = 0.0

        # water max values
        self[xxx.m_All_MBs, xxx.p_df] = +20.0
        self[xxx.m_Water, 0] = 100000.0

        # lock everything except water
        self.linklock[:] = 1
        self.linklock[xxx.m_Water, :] = 0
        return(self.copy())

    def set_default_min(self):
        """
        Initialize the params object to the minimum in vivo human values, no macromolecules.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        self._set_default_min()

        # browse though the database and find min values from literature
        cmin = []
        for this_metagroup_key, this_metagroup_entry in self._meta_bs.items():
            cmin.append(this_metagroup_entry["Concentration min"])

        # convert to np
        cmin = np.array(cmin)

        self[xxx.m_All_MBs, xxx.p_cm] = cmin[xxx.m_All_MBs]

        return(self.copy())

    def set_default_max(self):
        """
        Initialize the params object to the maximum in vivo values, no macromolecules.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        self._set_default_max()

        # browse though the database and find min values from literature
        cmax = []
        for this_metagroup_key, this_metagroup_entry in self._meta_bs.items():
            cmax.append(this_metagroup_entry["Concentration max"])

        # convert to np
        cmax = np.array(cmax)

        self[xxx.m_All_MBs, xxx.p_cm] = cmax[xxx.m_All_MBs]

        return(self.copy())

    def get_MBs_above_concentration(self, cm_threshold):
        """
        Return indexes for metabolites with average concentrations above a threshold.

        Parameters
        ----------
        cm_threshold : float
            Concentration threshold

        Returns
        -------
        im : list of integers
            Metabolite indexes
        """
        cm_mean = (self.set_default_min() + self.set_default_max()) / 2.0
        cm_mean = cm_mean[:, 0]
        metabolites_inds = np.where(cm_mean > cm_threshold)
        metabolites_inds = metabolites_inds[0]

        # display
        log.info("found %d metabolites above %f mmol/kg..." % (len(metabolites_inds), cm_threshold))
        meta_names = self.get_meta_names()
        found_meta_names = ""
        for im in metabolites_inds:
            found_meta_names += str(meta_names[im]) + " "
        log.info(found_meta_names)

        return(metabolites_inds)

    def add_macromolecules_min(self):
        """
        Enable the macromolecules modelization, minimum values.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        # macromolecules minimum values
        self[xxx.m_All_MMs, xxx.p_cm] = 0.001
        self[xxx.m_All_MMs, xxx.p_dd] = 100
        self[xxx.m_All_MMs, xxx.p_df] = -5
        self[xxx.m_All_MMs, xxx.p_dp] = -0.1

        # link all MM parameters to MM1
        self.linklock[xxx.m_All_MMs, :] = np.tile([0, 2000, 3000, 4000], (xxx.n_MMs, 1))
        self.linklock[xxx.m_Ref_MM, :] = [0, -2000, -3000, -4000]

        return(self.copy())

    def add_macromolecules_max(self):
        """
        Enable the macromolecules modelization, maximum values.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        self[xxx.m_All_MMs, xxx.p_cm] = 100.0
        self[xxx.m_All_MMs, xxx.p_dd] = 200.0
        self[xxx.m_All_MMs, xxx.p_df] = +5
        self[xxx.m_All_MMs, xxx.p_dp] = +0.1

        # link all MM parameters to MM1
        self.linklock[xxx.m_All_MMs, :] = np.tile([0, 2000, 3000, 4000], (xxx.n_MMs, 1))
        self.linklock[xxx.m_Ref_MM, :] = [0, -2000, -3000, -4000]
        return(self.copy())

    def print(self, bMM=False, bLL=True):
        """
        Display an array of parameters in the console.

        Parameters
        ----------
        bMM : boolean
            Includes macromolecular parameters (True) or not (False)
        bLL : boolean
            Displays link-lock status for each parameter (True) or not (False)
        """
        cell_nchar = 11

        log.info("displaying parameters...")
        log.info_line________________________()
        if(bMM):
            n = xxx.n_All
        else:
            n = xxx.n_MBs

        meta_names = self.get_meta_names()
        LLtransTab = str.maketrans("-+0123456789", "⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹")
        LLcolors = ['white', 'grey', 'yellow', 'blue', 'magenta', 'cyan', 'white', 'grey', 'yellow', 'blue', 'magenta', 'cyan', 'green']

        print("[#] " + "Metabolite".ljust(cell_nchar) + "[cm]".ljust(cell_nchar) + "[dd]".ljust(cell_nchar) + "[df]".ljust(cell_nchar) + "[dp]".ljust(cell_nchar))
        for k in range(n):
            print(("#%2d " % k) + meta_names[k].ljust(cell_nchar), end="")
            for kp in range(4):
                if(bLL):
                    thisLL = self.linklock[k, kp]
                    if(thisLL > 0):
                        thisLL_str = "+%1.0f" % thisLL
                    elif(thisLL == 0):
                        thisLL_str = "%1.0f " % thisLL
                    else:
                        thisLL_str = "%1.0f" % thisLL

                    thisLL_str = thisLL_str.translate(LLtransTab)
                    if(thisLL == +1.0):
                        thisLL_color = 'red'
                    else:
                        ind_color = -int(np.mod(np.abs(thisLL), len(LLcolors)))
                        thisLL_color = LLcolors[ind_color]

                    this_cell_str = ("(%4.1f)" % self[k, kp]) + thisLL_str
                    cprint(this_cell_str.ljust(cell_nchar), thisLL_color, attrs=['bold'], end="")
                else:
                    print(("(%4.1f)" % self[k, kp]).ljust(cell_nchar), end="")

            print("", flush=True)
        log.info_line________________________()

    def __reduce__(self):
        """Reduce internal pickling method used when dumping. Copied from MRSData2 class. Modified so that params attributes are not forgotten. See for more info: https://docs.python.org/3/library/pickle.html ."""
        # get numpy reduce tuple
        rd = super().__reduce__()
        # add params attributes
        rd2 = rd[2] + (self.__dict__,)
        # return the new reduce tuple version
        return(rd[0], rd[1], rd2)

    def __setstate__(self, d):
        """Set new state to object. Internal pickling method used when loading. Copied from MRSData2 class. Modified so that params attributes are not forgotten. See for more info: https://docs.python.org/3/library/pickle.html ."""
        # load params attributes
        self.__dict__ = d[-1]
        # load all the rest relative to numpy
        super().__setstate__(d[0:-1])
        return(self)


class mrs_sequence:
    """A class that stores a sequence and all its parameters used for simulation. This is a generic sequence class that you need to overload. By default, the simulated sequence is a simple pulse-acquire NMR experiment."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0):
        """
        Initialize the sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        """
        # all sequences have those parameters
        # name
        self.name = "fid"
        # type
        self.exc_type = sequence_exc_type.PULSE_ACQUIRE
        # echo time (ms)
        self.te = te
        # repetition time (ms)
        self.tr = tr
        # number of averages
        self.na = na
        # number of dummy scans
        self.ds = ds
        # which nuclei we are pulsing on/looking at (examples: '1H', '31P')
        self.nuclei = nuclei
        # number of acquired time points (int)
        self.npts = npts
        # voxel dimensions
        self.voxel_size = voxel_size
        # sampling frequency (Hz)
        self.fs = fs
        # larmor frequency of water (MHz)
        self.f0 = f0
        # reference voltage (V)
        self.vref = vref
        # shim vector
        self.shims = shims
        # start timestamp
        self.timestamp = timestamp
        # gating mode
        self.gating_mode = gating_mode
        # effective acquisition time
        self.eff_acquisition_time = eff_acquisition_time
        # kind of receiver gain
        self.scaling_factor = scaling_factor
        # ppm shift (ppm)
        self.ppm_water = 4.7
        # some 0th phase to add? (rd)
        self.additional_phi0 = 0.0

        # metabolite_basis_set object
        self._meta_bs = metabolite_basis_set()

        # try to simplify spin systems when possible to speed up simulations
        self.allow_spin_system_simplification = True
        # NMR simulation option: when hard zero-duration RF pulse are employed, should we take into account the duration of the real pulses in the evolution times or not? Experimentally, looks like yes.
        self.allow_evolution_during_hard_pulses = True

        # in case GAMMA could not be loaded, we can load the sequence from a stored file
        self.seqdb_file = None

        # pre-calculated stuff
        # 'metabase': set of numerically computed metabolite FID signals
        self._meta_signals = None
        # time vector
        self._t = []
        # last parameter call
        self._last_params = None
        self._last_model = None

        # initialized or not
        self._ready = False

    def copy(self):
        """Copy method"""
        obj = copy.copy(self)
        if(self._meta_signals is not None):
            obj._meta_signals = self._meta_signals.copy()
        if(self._t is not None):
            obj._t = self._t.copy()
        if(self._last_params is not None):
            obj._last_params = self._last_params.copy()
        if(self._last_model is not None):
            obj._last_model = self._last_model.copy()

        return(obj)

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
    def meta_bs(self):
        """
        Property get function for meta_bs.

        Returns
        -------
        self._meta_bs : metabolite_basis_set object
            Metabolite database to use for simulation
        """
        return(self._meta_bs)

    def _init_pulses(self):
        """Virtual method which initialize RF pulse waveforms if any."""

    def _prepare_spin_system(self, metabolite):
        """
        Return pyGAMMA spin system for a given metabolite, knowing all its properties. Simplify the system in simple cases like singulets (one single chemical shift, no J-couplings).

        Parameters
        ----------
        metabolite : dict
            metabolite_basis_set entry for one single metabolite

        Returns
        -------
        sys : pyGAMMA system object
            Spin system object used for simulation
        scaling_factor : float
            Scaling factor after system simplification
        """
        log.debug("preparing spin system")

        # init
        scaling_factor = 1.0

        # extract metabolite properties needed to create spin system
        ppm = metabolite["ppm"]
        iso = metabolite["iso"]
        j = metabolite["J"]

        # check if we can simplify
        if(self.allow_spin_system_simplification and len(ppm) > 1 and not np.any(j != 0.0) and len(np.unique(ppm)) == 1):
            log.debug("simplifying the spin system")
            # we have a non-coupled singulet with N spins here
            # let's simplify to one single spin + amplification factor
            scaling_factor = float(len(ppm))
            ppm = np.array([ppm[0]])
            iso = np.array([iso[0]])
            j = np.array([j[0, 0]])

        # init system
        nSpins = len(ppm)
        sys = pg.spin_system(nSpins)
        sys.Omega(self.f0)

        # for each spin
        for i_spin in range(nSpins):
            # set the nuclei
            if(iso[i_spin] == 1):
                sys.isotope(i_spin, '1H')
            elif(iso[i_spin] == 14):
                sys.isotope(i_spin, '14N')
            elif(iso[i_spin] == 31):
                sys.isotope(i_spin, '31P')
            else:
                log.error(str(iso[i_spin]) + ", that is weird nuclei!")

            # set the ppm
            sys.PPM(i_spin, ppm[i_spin])

            # set the couplings
            for icolJ in range(i_spin + 1, len(ppm)):
                sys.J(i_spin, icolJ, j[i_spin, icolJ])

        # water shift
        sys.offsetShifts(self.ppm_water * self.f0, self.nuclei)

        return(sys, scaling_factor)

    def _run_sequence(self, metabolite):
        """
        Simulate the NMR signal acquired using a single-pulse experiment using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_basis_set entry for one single metabolite

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        log.debug("acquiring pulse-acquire sequence: (90)...")
        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite)

        # amplitude normalization
        afactor = self.scaling_factor / sys.HS()

        # coupling
        H = pg.Hcs(sys) + pg.HJ(sys)

        # detection stuff
        D = pg.Fm(sys, self.nuclei)
        dt = np.double(1 / self.fs)

        # run the pulse-acquire experiment
        sigma0 = pg.sigma_eq(sys)
        te_real = 0.0
        # excitation: hard 90 pulse
        sigma1 = pg.Iypuls(sys, sigma0, self.nuclei, 90.0)
        te_real += 0.0
        # evolution
        sigma0 = pg.evolve(sigma1, pg.prop(H, self.te / 1000.0))  # TE evolution
        te_real += self.te / 1000.0
        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        log.debug("done. (real TE=%.2fms)" % (te_real * 1000.0))

        return(s_MRSData2)

    def _compute_metabolite_signals(self):
        """For each goup of metabolites, generate the simulated MRS signal."""
        # clear metabase
        self._meta_signals = []

        # zero signal vector
        s_full_meta = suspect.MRSData(np.zeros([self.npts, ]), 1 / self.fs, self.f0)
        s_full_meta = s_full_meta.view(reco.MRSData2)

        # browse though the database and display everything
        for this_metagroup_key, this_metagroup_entry in self._meta_bs.items():
            s = s_full_meta.copy()
            for this_meta_key, this_meta_entry in self._meta_bs[this_metagroup_key]["metabolites"].items():
                log.info("simulating MRS signal for metabolite [%s/%s]..." % (this_metagroup_key, this_meta_key))
                # build up the metabolite group
                s = s + self._run_sequence(this_meta_entry)
            # append this metabolite group to the metabase
            self._meta_signals.append(s)

    def _load_from_seqdb_file(self, te_tol=5.0, f0_tol=5.0):
        """
        Try to load the simulated metabolite signals from a stored PKL file.

        Parameters
        ----------
        te_tol : float
            Maximum TE difference tolerated when looking for a sequence (ms)
        f0_tol : float
            Maximum f0 difference tolerated when looking for a sequence (MHz)
        """
        if(self.seqdb_file is None):
            log.error("pyGAMMA library could not be loaded and no sequence database PKL file (seqdb_file) was specified :(")

        log.info("reading sequence database file...")
        # load pickle file
        with open(self.seqdb_file, 'rb') as f:
            [seqdb] = pickle.load(f)

        # convert to np
        seqdb = np.array(seqdb)

        # compare sequences
        log.info("trying to find the right simulated sequence for you...")

        # sequence excitation type
        seqdb_exc_type_mask = np.array([s.exc_type == self.exc_type for s in seqdb])
        seqdb_exc_type_n = np.sum(seqdb_exc_type_mask)
        log.info("sequence type match: n=%d (%s)" % (seqdb_exc_type_n, self.exc_type))

        # sequence name
        seqdb_name_mask = np.array([s.name == self.name for s in seqdb])
        seqdb_name_n = np.sum(seqdb_name_mask)
        log.info("sequence name match: n=%d (%s)" % (seqdb_name_n, self.name))

        # initialize final mask
        if(seqdb_name_n > 0):
            seqdb_final_mask = seqdb_name_mask
        else:
            seqdb_final_mask = seqdb_exc_type_mask

        # nuclei
        seqdb_nuclei_mask = np.array([s.nuclei == self.nuclei for s in seqdb])
        seqdb_nuclei_n = np.sum(seqdb_nuclei_mask)
        log.info("nuclei match: n=%d (%s)" % (seqdb_nuclei_n, self.nuclei))
        # adjust final mask
        seqdb_final_mask = np.logical_and(seqdb_final_mask, seqdb_nuclei_mask)

        # f0
        seqdb_f0_diff = np.abs([s.f0 - self.f0 for s in seqdb])
        seqdb_f0_diff_mask = (seqdb_f0_diff < f0_tol)
        seqdb_f0_diff_n = np.sum(seqdb_f0_diff_mask)
        log.info("f0 match: need %.3fMHz, found n=%d sequences simulated at +/-%.3fMHz" % (self.f0, seqdb_f0_diff_n, f0_tol))
        # adjust final mask
        seqdb_final_mask = np.logical_and(seqdb_final_mask, seqdb_f0_diff_mask)

        # ppm water
        seqdb_ppm_water_mask = np.array([s.ppm_water == self.ppm_water for s in seqdb])
        seqdb_ppm_water_n = np.sum(seqdb_ppm_water_mask)
        log.info("water ppm match: n=%d (%.2f)" % (seqdb_ppm_water_n, self.ppm_water))
        # adjust final mask
        seqdb_final_mask = np.logical_and(seqdb_final_mask, seqdb_ppm_water_mask)

        # additional_phi0
        seqdb_ph0_mask = np.array([s.additional_phi0 == self.additional_phi0 for s in seqdb])
        seqdb_ph0_n = np.sum(seqdb_ph0_mask)
        log.info("0th order phase match: n=%d (%.2f)" % (seqdb_ph0_n, self.additional_phi0))
        # adjust final mask
        seqdb_final_mask = np.logical_and(seqdb_final_mask, seqdb_ph0_mask)

        # allow_evolution_during_hard_pulses
        seqdb_evol_mask = np.array([s.allow_evolution_during_hard_pulses == self.allow_evolution_during_hard_pulses for s in seqdb])
        seqdb_evol_n = np.sum(seqdb_evol_mask)
        log.info("complicated stuff match: n=%d (%r)" % (seqdb_evol_n, self.allow_evolution_during_hard_pulses))
        # adjust final mask
        seqdb_final_mask = np.logical_and(seqdb_final_mask, seqdb_evol_mask)

        # metabolites
        seqdb_meta_bs_mask = np.array([s.meta_bs == self.meta_bs for s in seqdb])
        seqdb_meta_bs_n = np.sum(seqdb_meta_bs_mask)
        log.info("metabolites match: n=%d" % seqdb_meta_bs_n)
        # adjust final mask
        seqdb_final_mask = np.logical_and(seqdb_final_mask, seqdb_meta_bs_mask)

        # check mask
        if(not np.any(seqdb_final_mask)):
            log.error("sorry, there is no simulated sequence that matches with your acquisition criterias :(")

        # mask all the useless sequences and clean memory
        seqdb_match = seqdb[seqdb_final_mask]
        seqdb = []

        # te
        seqdb_te_diff = np.abs([s.te - self.te for s in seqdb_match])
        imin_te_diff = np.argmin(seqdb_te_diff)
        log.info("TE match: need %.2fms, found %.2fms" % (self.te, seqdb_match[imin_te_diff].te))
        if(seqdb_te_diff[imin_te_diff] > te_tol and self.meta_bs.non_coupled_only):
            log.info("TE match: This might seem a lot but I see you want to generate only singulet's metabolites so that is really not a big deal ;)")

        # optimal sequence and clean memory
        optim_seq = seqdb_match[imin_te_diff]
        seqdb_match = []

        # display
        log.info("comparing what you asked/what you got...")
        cell_nchar = 20
        log.info("parameter".ljust(cell_nchar) + " requested".ljust(cell_nchar) + " obtained".ljust(cell_nchar))

        unique_key_list = list(set(list(self.__dict__.keys()) + list(optim_seq.__dict__.keys())))
        for this_key in unique_key_list:
            my_val = self.__dict__[this_key]
            my_val_len = 1
            try:
                my_val_len = len(my_val)
            except:
                pass

            if(this_key[0] != '_' and my_val_len == 1):
                found_val = optim_seq.__dict__[this_key]
                # format floats
                if(type(my_val) == float):
                    my_val_str = "%.2f" % my_val
                    found_val_str = "%.2f" % found_val
                else:
                    my_val_str = str(my_val)
                    found_val_str = str(found_val)

                # crop strings
                this_key_str = this_key[0:cell_nchar]
                my_val_str = my_val_str[0:cell_nchar]
                found_val_str = found_val_str[0:cell_nchar]
                # pretty print
                log.info(this_key_str.ljust(cell_nchar) + " " + my_val_str.ljust(cell_nchar) + " " + found_val_str.ljust(cell_nchar))

        # ok well done. Now we maybe have to fix a few issues: number of samples, sampling frequency, amplification factor, this can be done with some signal processing stuff

        # resampling (even if not needed)
        log.info("resampling metabolite signals: %dpts/%.2fHz to %dpts/%.2fHz..." % (optim_seq.npts, optim_seq.fs, self.npts, self.fs))
        log.info("rescaling metabolite signals by a factor of %.2f..." % (self.scaling_factor / optim_seq.scaling_factor))
        old_dt = 1 / optim_seq.fs
        old_t = np.arange(0, optim_seq.npts * old_dt, old_dt)
        new_dt = 1 / self.fs
        new_t = np.arange(0, self.npts * new_dt, new_dt)

        # resample each metabolite signal
        new_meta_signals = []
        for s in optim_seq._meta_signals:
            s2_np = np.interp(new_t, old_t, s)
            # convert to suspect
            s_MRSData = suspect.MRSData(s2_np, new_dt, optim_seq.f0)
            s_MRSData2 = s_MRSData.view(reco.MRSData2)
            # rescale
            s_MRSData2 = s_MRSData2 * self.scaling_factor / optim_seq.scaling_factor
            # and rebuild metabase
            new_meta_signals.append(s_MRSData2)

        # final: carefully copy attributes except some
        keys_not_to_copy = ["_meta_bs", "_meta_signals", "_t", "_ready", "_last_params", "_last_model", "npts", "fs"]
        for this_key in list(self.__dict__.keys()):
            if(this_key not in keys_not_to_copy):
                self.__dict__[this_key] = optim_seq.__dict__[this_key]

        # apply changes
        self._t = new_t.copy()
        self._meta_signals = new_meta_signals.copy()

        log.info("successfully imported metabolite signal basis set! :)")

    def initialize(self, meta_bs=None):
        """
        Initialize the sequence object before using it to simulate MRS signals.

        Parameters
        ----------
        meta_bs : metabolite_basis_set
            Metabolite database to use for simulation. If none specified, use default metabolite_basis_set object.
        """
        # want to use a custom metabolite db?
        if(meta_bs is not None):
            self._meta_bs = meta_bs

        # now let's initialize the metabolite db if needed
        if(not self._meta_bs.ready):
            self._meta_bs.initialize()

        # wait, was the GAMMA library imported ok?
        if(GAMMA_LIB_LOADED):
            log.info("initializing sequence using pyGAMMA...")
            # oh ok, so let's run the simulations and stuff
            self._init_pulses()
            self._compute_metabolite_signals()
            dt = 1 / self.fs
            self._t = np.arange(0, self.npts * dt, dt)
        else:
            log.info("loading sequence from disk...")
            # ops, so let's try to load the simulations from a pkl file
            self._load_from_seqdb_file()

        # initialize the model persistent memory
        self._last_params = None
        self._last_model = None

        self._ready = True

    def _model(self, p):
        """
        Most important function in here... Simulates the MRS signal using a parametric model that allows the control in amplitude, linewidth, frequency and phase shift for each metabolite.

        Parameters
        ----------
        p : params object
            Array of simulation parameters

        Returns
        -------
        s : MRSData2 numpy array [timepoints]
            Modeled MRS data signal stored in a MRSData2 object
        """
        # ready or not, here I come
        if(not self.ready):
            log.error("this mrs_sequence object was not initialized!")

        # check last call
        if(self._last_params is not None and np.all(self._last_params == p)):
            # we are asking for the same model than last time
            # calling back persistent momery
            return(self._last_model)

        # time MRS model
        s_MRSData = suspect.MRSData(np.zeros([self.npts, ]), 1 / self.fs, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)
        for k, s_meta in enumerate(self._meta_signals):
            s_MRSData2 = s_MRSData2 + s_meta * p[k, 0] * np.exp((-p[k, 1] + 2.0 * np.pi * 1j * p[k, 2]) * self._t + 1j * (p[k, 3] + self.additional_phi0))

        # remember
        self._last_params = p
        self._last_model = s_MRSData2

        return(s_MRSData2)

    def _jac(self, p):
        """
        Return the jacobian vector: helps optimization, need for CRBs estimation and parameter correlation matrix.

        Parameters
        ----------
        p : params object
            Array of simulation parameters

        Returns
        -------
        j : numpy array [number of metabolites,number of parameters,number of time points]
            Jacobian vector, has to be reformatted before use with fit
        """
        # ready or not, here I come
        if(not self.ready):
            log.error("this mrs_sequence object was not initialized!")

        # for each metabolite, caculate derivative for each parameter
        j = np.full((len(self._meta_signals), 4, self.npts), 0.0, dtype=np.complex128)
        for k, s_meta in enumerate(self._meta_signals):
            # for each parameter, calculate derivative
            # pre-calculate exp term for this metabolite
            exp_term = np.exp((-p[k, 1] + 2.0 * np.pi * 1j * p[k, 2]) * self._t + 1j * (p[k, 3] + self.additional_phi0))
            # concentration
            j[k, 0, :] = s_meta * exp_term
            # damping
            j[k, 1, :] = -p[k, 0] * s_meta * self._t * exp_term
            # freq
            j[k, 2, :] = p[k, 0] * s_meta * 2.0 * np.pi * 1j * self._t * exp_term
            # phase
            j[k, 3, :] = p[k, 0] * s_meta * 1j * exp_term

        return(j)

    def simulate_signal(self, p, sigma_noise=0.0, na=1, lbl="simulated MRS signal"):
        """
        Print out the parameter values and returns the modeled signal using above member function.

        Parameters
        ----------
        p : params object
            Array of simulation parameters
        sigma_noise : float
            Noise level
        na : int
            Number of averages
        lbl : string
            Label describing simulated data

        Returns
        -------
        s : MRSData2 numpy array [timepoints]
            Simulated MRS data signal stored in a MRSData2 object
        """
        # ready or not, here I come
        if(not self.ready):
            log.error("this mrs_sequence object was not initialized!")

        # checking if LL is not broken
        if(not p.check()):
            log.error("the link-lock vector of this params object looks broken!")

        log.info("simulating signal...")
        s = self._model(p)

        if(sigma_noise > 0.0 and na == 1):
            # simple noisy 1 shot simulation
            log.info("adding complex gaussian noise, std. deviation = " + str(sigma_noise))
            s = s + np.random.normal(0.0, sigma_noise, s.shape[0] * 2).view(np.complex128)

        elif(sigma_noise > 0.0 and na > 1):
            # averaged noisy simulation
            b = log.progressbar("averaging", na)
            s_single_shot = s.copy()
            s_averaged = s.copy() * 0.0
            for a in range(na):
                s_averaged += s_single_shot + np.random.normal(0.0, sigma_noise, s.shape[0] * 2).view(np.complex128)
                b.update(a)
            b.finish("done")
            s = s_averaged
        else:
            log.warning("ambiguous arguments: simulating a single-shot signal without noise!")

        # adding a few attributes
        s._te = self.te
        s._sequence = self
        s._noise_level = s.analyze_noise_nd()
        s.set_display_label(lbl)
        return(s)


class mrs_seq_press(mrs_sequence):
    """A class that represents a general PRESS sequence."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, te1=np.nan, te2=np.nan):
        """
        Initialize a virtual PRESS sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        te1 : float
            First part of TE (ms)
        te2 : float
            Second part of TE (ms)
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor)
        # name of sequence
        self.name = "press (not specific)"
        # type
        self.exc_type = sequence_exc_type.SPIN_ECHO
        # TE timing
        self.te = te
        if(te1 is not np.nan and te2 is not np.nan):
            self.te1 = te1
            self.te2 = te2
        else:
            self.te1 = np.nan
            self.te2 = np.nan

        # flip phase for PRESS
        self.additional_phi0 = np.pi

        # freeze
        self.__isfrozen = True

    def _run_sequence(self, metabolite):
        """
        Simulate the NMR signal acquired using a PRESS sequence using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_basis_set entry for one single metabolite

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        log.debug("acquiring PRESS sequence: (90)-(180)-(180)...")
        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite)

        # amplitude normalization
        afactor = self.scaling_factor / sys.HS()

        # coupling
        H = pg.Hcs(sys) + pg.HJ(sys)

        # detection stuff
        D = pg.Fm(sys, self.nuclei)
        dt = np.double(1 / self.fs)

        # PRESS timing implementation
        # 90-a-180-b-c-180-d-FID
        # all in ms
        if(self.te1 is np.nan and self.te2 is np.nan):
            # not TE1/TE2 specified, assume a symmetric scheme
            a = self.te / 4.0
            bc = self.te / 2.0
            d = self.te / 4.0
        else:
            # TE1/TE2 timing
            ab = self.te1
            cd = self.te2
            a = ab / 2.0
            bc = ab / 2.0 + cd / 2.0
            d = cd / 2.0

        # delay list in seconds
        evol_delays_s = np.array([a, bc, d]) / 1000.0
        # flip angle list
        flip_angles_deg = np.array([90.0, 180.0, 180.0])

        # run the sequence
        sigma0 = pg.sigma_eq(sys)
        te_real = 0.0

        # excitation: hard 90 pulse
        for d, p in zip(evol_delays_s, flip_angles_deg):
            sigma1 = pg.Iypuls(sys, sigma0, self.nuclei, p)
            te_real += 0.0
            # evolution
            sigma0 = pg.evolve(sigma1, pg.prop(H, d))
            te_real += d

        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        log.debug("done. (real TE=%.2fms)" % (te_real * 1000.0))

        return(s_MRSData2)


class mrs_seq_steam(mrs_sequence):
    """A class that represents a general STEAM sequence."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor)
        # name of sequence
        self.name = "steam (not specific)"
        # type
        self.exc_type = sequence_exc_type.STIMULATED_ECHO
        # mixing time
        self.tm = tm
        # need some 180deg phase here
        self.additional_phi0 = np.pi

        # freeze
        self.__isfrozen = True

    def _run_sequence(self, metabolite):
        """
        Simulate the NMR signal acquired using a STEAM sequence using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_basis_set entry for one single metabolite

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        log.debug("acquiring PRESS sequence: (90)-(90)-(90)...")
        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite)

        # amplitude normalization
        afactor = self.scaling_factor / sys.HS()

        # coupling
        H = pg.Hcs(sys) + pg.HJ(sys)

        # detection stuff
        D = pg.Fm(sys, self.nuclei)
        dt = np.double(1 / self.fs)

        # STEAM timing implementation
        # 90-a-90-b-90-c-FID
        # all in ms
        a = self.te / 2.0
        b = self.tm
        c = self.te / 2.0

        # delay list in seconds
        evol_delays_s = np.array([a, b, c]) / 1000.0

        # run the sequence
        sigma0 = pg.sigma_eq(sys)
        te_real = 0.0

        # 1st 90 pulse
        sigma1 = pg.Iypuls(sys, sigma0, self.nuclei, 90.0)
        te_real += 0.0
        # evolution
        sigma0 = pg.evolve(sigma1, pg.prop(H, evol_delays_s[0]))
        te_real += evol_delays_s[0]

        # the following code is from Vespa's pulse_sequences.xml file
        # STEAM Ideal version 2 from Brian Soher
        # seems that whatever TM we have, it will give the same result

        # Now we need to create the effect of crushers around the 2nd and 3rd
        # 90 pulses. This is done by creating 4 copies of spin state and repeating
        # the rest of the sequence for four different rotations around z-axis
        dephase_ang = [0.0, 90.0, 180.0, 270.0]
        sigma_mult = []
        for i in dephase_ang:
            sigma_mult.append(pg.gen_op(sigma0))

        for i, angle in enumerate(dephase_ang):
            # calculate and apply rotation around z-axis
            riz = pg.gen_op(pg.Rz(sys, angle))
            sigma_mult[i] = pg.evolve(sigma_mult[i], riz)

            # second 90 pulse
            sigma_mult[i] = pg.Ixpuls(sys, sigma_mult[i], self.nuclei, 90.0)

            # this function removes all coherences still in transverse plane
            # this removes all stimulated echos from first and second 90 pulse
            pg.zero_mqc(sys, sigma_mult[i], 0, -1)

            # third 90 pulse
            sigma_mult[i] = pg.Ixpuls(sys, sigma_mult[i], self.nuclei, 90.0)
            # undo rotation around z-axis
            sigma_mult[i] = pg.evolve(sigma_mult[i], riz)
            # scale results based on the number of phase angles
            sigma_mult[i] *= 1.0 / float(len(dephase_ang))

            # sum up each rotated/unrotated results
            if i == 0:
                sigma_res = pg.gen_op(sigma_mult[i])
            else:
                sigma_res += sigma_mult[i]

        # last TE/2 nutation
        sigma0 = pg.evolve(sigma_res, pg.prop(H, evol_delays_s[2]))
        te_real += evol_delays_s[2]

        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        log.debug("done. (real TE=%.2fms)" % (te_real * 1000.0))

        return(s_MRSData2)


class mrs_seq_eja_svs_slaser(mrs_sequence):
    """A class that represents the semi-LASER sequence from CMRR."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, exc_pulse_duration=2.0, exc_pulse_voltage=350.0, rfc_pulse_duration=9.0, rfc_pulse_fa=180.0, rfc_pulse_r=20.0, rfc_pulse_n=1.0, rfc_pulse_voltage=350.0, ref_pulse_voltage=300.0, spoiler_duration=1.0):
        """
        Initialize a virtual semi-LASER sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        exc_pulse_duration : float
            Excitation pulse duration (ms)
        exc_pulse_voltage : float
            Excitation pulse voltage (V)
        rfc_pulse_duration : float
            Refocussing pulse duration (ms)
        rfc_pulse_fa : float
            Refocussing pulse flip angle, which does not mean much when dealing with HSn pulses (deg)
        rfc_pulse_r : float
            AHP HSn refocussing pulse R (Hz.s)
        rfc_pulse_n : float
            AHP HSn refocussing pulse N
        rfc_pulse_voltage : float
            AHP HSn refocussing pulse voltage (V)
        spoiler_duration : float
            Spoiler duration (ms)
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor)
        # name of sequence
        self.name = "eja_svs_slaser"
        # type
        self.exc_type = sequence_exc_type.SPIN_ECHO
        # excitation and AFP pulse properties
        self.pulse_exc_duration = exc_pulse_duration
        self.pulse_rfc_duration = rfc_pulse_duration
        self.pulse_rfc_flipangle = rfc_pulse_fa
        self.pulse_rfc_r = rfc_pulse_r
        self.pulse_rfc_n = rfc_pulse_n
        self.pulse_rfc_voltage = rfc_pulse_voltage
        self.pulse_rfc_npts = None  # no default number of points
        # spoiler properties
        self.spoiler_duration = spoiler_duration

        # should we use real shaped RF pulses or not
        self.pulse_rfc_real_shape_enable = False

        # sLASER needs some rephasing
        self.additional_phi0 = np.pi

        # --- RF power adjustement ---
        # adjust RF by sweeping B1 power
        self.pulse_rfc_optim_power_enable = False
        # voltage range to sweep (V)
        self.pulse_rfc_optim_power_voltage_range = [0, 500]
        # number of measurements in the previous range
        self.pulse_rfc_optim_power_voltage_n = 50
        # which metabolites should we test (separately)
        self.pulse_rfc_optim_power_metabolites = [xxx.m_Water, xxx.m_Eth]
        # margins/threshold parameters to estimate optimal RF power. By default:
        # [0] 10% last points (for higher power) used for plateau estimation
        # [1] 10% change allowed
        # [2] 100% power increase to be sure to be in adiabatic regime
        self.pulse_rfc_optim_power_margins = [10.0, 10.0, 100.0]  # %
        # display all this in a nice fig
        self.pulse_rfc_optim_power_display = True
        # final pulse voltage found by optimization
        self.pulse_rfc_optim_power_voltage = None
        # final pulse flip angle found by optimization
        self.pulse_rfc_optim_power_flipangle = None

        # RF power (w1 in Hz) used for AFP pulses during simulation
        self._pulse_rfc_w1max = None

        # interpulse delay lists
        self._subdelay_c2c_ms_list = []
        self._subdelay_evol_ms_list = []
        self.tcp_ms_list = []

        # freeze
        self.__isfrozen = True

    def get_csde(self, voxel_size_mm_list=None):
        """
        Calculate chemical shift displacement in three directions.

        Returns
        -------
        csd : numpy array of floats
            Estimated CS displacement in %/ppm
        """
        # init
        log.debug("estimating chemical shift displacement error...")
        if(voxel_size_mm_list is None):
            # TODO: I am not smart enough to make this work properly. Voxel orientation is sometimes modified meaning that RF pulses #1, #2 and #3 do not always select along RL, AP, HF or X, Y, Z. I am not very good with transformation matrix and stuff... :( So for now, I give the possibility to the user to input his voxel sizes in the order of pulsing.
            voxel_size_mm_list = self.voxel_size

        csd = [None, None, None]
        df_abs_Hz = self.f0  # yes, 1ppm==f0[MHz]/1e6=297Hz at 7T

        # assuming X-Y-Z is done with 90-180-180
        log.info("estimating CS displacement for semiLASER: assuming (90x)-(180y)-(180z)!...")

        # X selection done with asymmetric 90° pulse
        # we do not know much about this pulse. We can only say it is 3.4kHz large if the duration is 2ms
        log.debug("estimating CS displacement for (90x): this pulse is the weird asymmetric one, we have no idea what it is exactly!")
        if(self.pulse_laser_exc_length == 2000.0):
            log.debug("since its duration is 2ms here, we assume, according to Oz & Tkac, MRM 65:901-910 (2011), that its bandwidth is 3.4kHz.")
            bw_x_Hz = 3400.0
            grad_x_Hz_m = bw_x_Hz / (voxel_size_mm_list[0] * 0.001)
            d_x_mm = 1000.0 * df_abs_Hz / grad_x_Hz_m
            d_x_prct = d_x_mm / voxel_size_mm_list[0] * 100.0
        else:
            log.debug("since its duration is not 2ms here, we do not know its bandwith. Therefore, no way to calculate the CS displacement for this axis, sorry ;)")
            d_x_m = None

        # Y selection done with 180°
        bw_y_Hz = self.pulse_laser_rfc_r / (self.pulse_laser_rfc_length / 1000000.0)
        grad_y_Hz_m = bw_y_Hz / (voxel_size_mm_list[1] * 0.001)
        d_y_mm = 1000.0 * df_abs_Hz / grad_y_Hz_m
        d_y_prct = d_y_mm / voxel_size_mm_list[1] * 100.0

        # Z selection done with 180°
        bw_z_Hz = self.pulse_laser_rfc_r / (self.pulse_laser_rfc_length / 1000000.0)
        grad_z_Hz_m = bw_z_Hz / (voxel_size_mm_list[2] * 0.001)
        d_z_mm = 1000.0 * df_abs_Hz / grad_z_Hz_m
        d_z_prct = d_z_mm / voxel_size_mm_list[2] * 100.0

        csd = np.array([d_x_prct, d_y_prct, d_z_prct])

        return(csd)

    def _sech(self, x):
        return(1.0 / np.cosh(x))

    def _asech(self, x):
        return(ma.acosh(1 / x))

    def _rf_pulse_hsn(self, N: int, n: int, tbw, Tp, trunc=0.01, display=False):
        """
        Return a n-order hyperbolic secant (HSn) RF pulse. This function was 100% taken from the MATLAB FID-A (Jamie Near's code) toolbox and tranlated to Python. It is based on the method of gradient modulated offset-independent adiabaticity as first described by Tannus and Garwood in NMR Biomed 1997; 10:423-434.

        Parameters
        ----------
        N : int
            Number of points in RF waveform
        n : int
            Order of the HS pulse
        tbw : float
            Time bandwidth product (Hz*s)
        Tp : float
            Duration of the RF pulse (ms)
        trunc : float
            Truncation of the amplitude modulation function
        display : bool
            Display resulting time and frequency pulse profiles

        Returns
        -------
        pulse_amp  : float list
            RF amplitude waveform for a HSn pulse (normalized 0-1)
        pulse_phi  : float list
            RF phase waveform for a HSn pulse (deg)
        """
        # make sure N is even
        if(np.mod(N, 2) != 0):
            N = N + 1

        Tp_s = Tp / 1000.0

        # initialize the time vectors ta has N steps from 0 to Tp.
        t = np.linspace(0, Tp_s, int(N))

        # tau has N steps from -1 to 1. (useful for defining our AM and GM functions.)
        tau = t * 2 / Tp_s - 1

        # create truncation factor.
        B = self._asech(trunc)

        # Find Bandwith Factor A
        bw = tbw / Tp_s
        A = bw / 2

        # find time step size
        dt = t[1]

        # First define the AM function:
        F1 = self._sech(B * (np.power(tau, n)))

        # now calculate the FM function based on the assumption of a constant gradient by integrating the AM function:
        F2 = np.zeros(F1.shape)
        F2[int(F1.shape[0] / 2 + 1):F1.shape[0]] = np.cumsum(F1[int(F1.shape[0] / 2 + 1):])
        F2[0:int(F1.shape[0] / 2)] = -F2[F1.shape[0]:int(F1.shape[0] / 2 - 1):-1]
        F2 = F2 / np.max(F2)

        pulse_amp = F1
        FM = A * F2

        # create phase modulation function
        pulse_phi = np.cumsum(FM) * dt * 360.0

        # display
        if(display):
            # display
            fig = plt.figure(11)
            fig.clf()
            axs = fig.subplots(3, 1)
            fig.canvas.set_window_title("_rf_pulse_hsn (mrs.sim.mrs_seq_eja_svs_slaser)")
            fig.suptitle("HSn pulse generation | %.2f ms, BW = %.2f Hz, TBW/R = %.2f, n-order = %d, %d pts" % (Tp, bw, tbw, n, N))

            axs[0].plot(t * 1000.0, pulse_amp, 'k-', linewidth=1)
            axs[0].set_xlabel('time (ms)')
            axs[0].set_ylabel('pulse amplitude')
            axs[0].grid('on')

            axs[1].plot(t * 1000.0, pulse_phi, 'k-', linewidth=1)
            axs[1].set_ylabel('phase profile')
            axs[1].grid('on')

            # frequency profile
            f = np.fft.fftshift(np.fft.fftfreq(int(N), d=(t[1] - t[0])))
            pulse_cmplx = pulse_amp * np.exp(pulse_phi * 1j)
            pulse_freq_abs = np.abs(np.fft.fftshift(np.fft.fft(pulse_cmplx)))

            axs[2].plot(f / 1000.0, pulse_freq_abs, 'k-', linewidth=1)
            axs[2].set_xlabel('frequency (kHz)')
            axs[2].set_ylabel('frequency profile (magnitude)')
            axs[2].grid('on')

            fig.subplots_adjust()
            fig.show()

        return(pulse_amp, pulse_phi)

    def _init_pulses(self):
        """Calculate HSn pulse profiles for later use during simulation. This could look like a detail but the semi-LASER behaves in particular way regarding effective TE and J-coupling BECAUSE of those pulses and their long durations. Depending on TE and pulse properties, MR peak lineshapes can be severly affected."""
        if(self.pulse_rfc_real_shape_enable):
            # pulse bandwidth
            pulse_bw_hz = self.pulse_rfc_r / (self.pulse_rfc_duration / 1000.0)

            if(self.pulse_rfc_npts is None):
                # number of points: critical because pulse waveforms containing many points will slow down simulation! Let's check Nyquist criteria and choose a minium number of points...
                # nyquist sampling
                pulse_min_fs = pulse_bw_hz * 2.0
                pulse_n = pulse_min_fs * self.pulse_rfc_duration / 1000.0
                # multiply this by 2 (empiric finding)
                pulse_n = pulse_n * 2
            else:
                pulse_n = self.pulse_rfc_npts

            # call HSn generation
            log.debug("generating a %.1f-ms HS%d pulse: R=%.0f (BW=%.1fkHz) using %d points..." % (self.pulse_rfc_duration, self.pulse_rfc_n, self.pulse_rfc_r, pulse_bw_hz / 1000.0, pulse_n))
            pulse_amp, pulse_phi = self._rf_pulse_hsn(pulse_n, self.pulse_rfc_n, self.pulse_rfc_r, self.pulse_rfc_duration)

            # store
            self.pulses_rfc_shape_amp_norm = pulse_amp
            self.pulses_rfc_shape_phi_deg = pulse_phi

            if(self.pulse_rfc_optim_power_enable):
                self._pulse_rfc_w1max = self._get_rfc_pulse_w1max_by_optim()
            else:
                self._pulse_rfc_w1max = self._get_rfc_pulse_w1max_using_ref_pulse_voltage()
        else:
            log.debug("not using real shaped pulses so nothing to do here!")

    def _get_rfc_pulse_w1max_using_ref_pulse_voltage(self):
        """
        Return the RF power (w1max) of the HSn pulses using the reference voltage.

        Returns
        -------
        w1_rfc_pulse_hz  : float
            RF power or w1 max for HSn pulse (Hz)
        """
        # optimize RF power
        log.debug("setting HSn pulse RF power using reference voltage...")

        # pulse area
        pulse_t = np.linspace(0.0, 1.0 / 1000.0, len(self.pulses_rfc_shape_amp_norm))
        pulse_rfc_area = np.trapz(self.pulses_rfc_shape_amp_norm, pulse_t, pulse_t[1])

        # find real RF power that was used
        w1_ref_pulse_degs = 180 / (360 * 1e-3)  # deg/s
        w1_rfc_pulse_degs = self.pulse_rfc_voltage * w1_ref_pulse_degs / self.vref  # deg/s
        w1_rfc_pulse_hz = w1_rfc_pulse_degs / (360.0 * pulse_rfc_area)

        log.debug("HSn pulse voltage is set to %.2fV" % self.pulse_rfc_voltage)
        log.debug("HSn pulse flip angle is set to %.2fdeg" % self.pulse_rfc_flipangle)
        log.debug("reference voltage is %.2fV" % self.vref)
        log.debug("this means the HSn pulse w1max is %.2fHz" % w1_rfc_pulse_hz)

        return(w1_rfc_pulse_hz)

    def _get_rfc_pulse_w1max_by_optim(self):
        """
        Optimize the RF power (w1max) of the HSn pulses. The current power used during acquisition, if available, will be displayed.

        Returns
        -------
        w1_rfc_pulse_hz_optim : float
            RF power or w1 max for HSn pulse (Hz)
        """
        # optimize RF power
        log.debug("optimizing RF power...")

        # pulse area
        pulse_t = np.linspace(0.0, 1.0 / 1000.0, len(self.pulses_rfc_shape_amp_norm))
        pulse_rfc_area = np.trapz(self.pulses_rfc_shape_amp_norm, pulse_t, pulse_t[1])

        # find real RF power that was used
        w1_ref_pulse_degs = 180 / (360 * 1e-3)  # deg/s
        w1_rfc_pulse_degs = self.pulse_rfc_voltage * w1_ref_pulse_degs / self.vref  # deg/s
        w1_rfc_pulse_hz = w1_rfc_pulse_degs / (360.0 * pulse_rfc_area)

        # prepare RF power triple axis
        w1_range_V = np.linspace(self.pulse_rfc_optim_power_voltage_range[0], self.pulse_rfc_optim_power_voltage_range[1], self.pulse_rfc_optim_power_voltage_n)
        w1_range_Hz = (w1_range_V * w1_ref_pulse_degs / self.vref) / (360.0 * pulse_rfc_area)

        # in uT?
        # gamma = GAMMA_DICT[self.nuclei]
        # w1_range_uT = w1_range_Hz / gamma

        # run RF adjustment
        meta_bs_keys = list(self.meta_bs.keys())
        acquired_signals = []
        peak_max_ppm_index = []
        peak_intensity_abs = np.zeros([len(self.pulse_rfc_optim_power_metabolites), len(w1_range_Hz)])
        peak_intensity_real = np.zeros([len(self.pulse_rfc_optim_power_metabolites), len(w1_range_Hz)])

        # for each RF power (w1)
        for kw, w in enumerate(w1_range_Hz):
            self._pulse_rfc_w1max = w
            # for each metabolite
            for km, im in enumerate(self.pulse_rfc_optim_power_metabolites):
                meta_key = meta_bs_keys[im]
                meta_dict_entry = self.meta_bs[meta_key]["metabolites"][meta_key]
                # acquire
                s = self._run_sequence(meta_dict_entry)
                s = s.correct_apodization_nd(10.0)  # silent apodization
                acquired_signals.append(s)
                # analyze
                sf = s.spectrum()
                if(kw == 0):
                    # first power, let's peak pick
                    peak_max_ppm_index.append(np.argmax(np.abs(sf)))
                # store peak intensities
                peak_intensity_abs[km, kw] = np.abs(sf[peak_max_ppm_index[km]])
                peak_intensity_real[km, kw] = np.real(sf[peak_max_ppm_index[km]])

        # find optimal power
        # first, find plateau value by taking
        i_power = int(self.pulse_rfc_optim_power_voltage_n - self.pulse_rfc_optim_power_voltage_n * self.pulse_rfc_optim_power_margins[0] / 100.0)
        plateau_power = np.average(peak_intensity_abs[:, i_power:], axis=1).reshape([len(self.pulse_rfc_optim_power_metabolites), 1])
        peak_intensity_abs_rel = np.abs(peak_intensity_abs - plateau_power) / plateau_power * 100.0
        peak_intensity_abs_rel = np.mean(peak_intensity_abs_rel, axis=0)
        # second, find power for which variation is above 10%
        peak_intensity_abs_rel_above_threshold_mask = (peak_intensity_abs_rel > self.pulse_rfc_optim_power_margins[1])
        w1_range_Hz_optim = w1_range_Hz[peak_intensity_abs_rel_above_threshold_mask]
        w1_range_Hz_optim = w1_range_Hz_optim[-1]
        # third, increase optimal power of 10% for security (to be sure to be in adiabatic regime)
        w1_range_Hz_optim = w1_range_Hz_optim * (100 + self.pulse_rfc_optim_power_margins[2]) / 100.0

        # convert back to voltage (usefull for virtual pulse calibration)
        w1_range_degs_optim = w1_range_Hz_optim * 360.0 * pulse_rfc_area  # deg/s
        self.pulse_rfc_optim_power_voltage = w1_range_degs_optim * self.vref / w1_ref_pulse_degs  # V

        # what would that be in deg?
        # flip angles are linear with voltages
        self.pulse_rfc_optim_power_flipangle = self.pulse_rfc_flipangle * self.pulse_rfc_optim_power_voltage / self.pulse_rfc_voltage

        # peak ppm
        ppm = s.frequency_axis_ppm()
        peak_max_ppm_index = np.array(peak_max_ppm_index)
        peak_max_ppm = ppm[peak_max_ppm_index]

        # display
        fig = plt.figure(10)
        fig.clf()
        axs = fig.subplots(1, 2)
        fig.canvas.set_window_title("_get_rfc_pulse_w1max_by_optim (mrs.sim.mrs_seq_eja_svs_slaser)")
        fig.suptitle("RF power optimization results for sLASER")

        afactor = (w1_range_Hz[2] - w1_range_Hz[1]) / np.max(peak_intensity_abs)
        i = 0
        # for each RF power (w1)
        for kw, w in enumerate(w1_range_Hz):
            # for each metabolite
            for km, im in enumerate(self.pulse_rfc_optim_power_metabolites):
                s = acquired_signals[i]
                axs[0].plot(s.frequency_axis_ppm(), w + afactor * s.spectrum().real, 'k', linewidth=1)
                i = i + 1

        # for each metabolite, draw ppm measurement line
        for km, im in enumerate(self.pulse_rfc_optim_power_metabolites):
            axs[0].plot([peak_max_ppm[km], peak_max_ppm[km]], [np.min(w1_range_Hz), np.max(w1_range_Hz)], 'r', linewidth=1)

        axs[0].set_xlabel('chemical shift (ppm)')
        axs[0].set_ylabel('w1 (Hz)')
        axs[0].grid('on')
        axs[0].set_xlim(5, 1)

        # for metabolite
        for km, im in enumerate(self.pulse_rfc_optim_power_metabolites):
            meta_key = meta_bs_keys[im]
            axs[1].plot(peak_intensity_abs[km, :], w1_range_Hz, 'o-', linewidth=1, label=meta_key + " peak intensity at %.2fppm (magnitude)" % peak_max_ppm[km])
            axs[1].plot(peak_intensity_real[km, :], w1_range_Hz, 'o-', linewidth=1, label=meta_key + " peak intensity at %.2fppm (real)" % peak_max_ppm[km])

        # add real RF power as a line
        axs[1].plot([np.min(peak_intensity_real), np.max(peak_intensity_abs)], [w1_rfc_pulse_hz, w1_rfc_pulse_hz], linewidth=3, label="Current RF power")

        # add optimal RF power as a line
        axs[1].plot([np.min(peak_intensity_real), np.max(peak_intensity_abs)], [w1_range_Hz_optim, w1_range_Hz_optim], linewidth=3, label="Optimal RF power")

        axs[1].set_xlabel('peak intensity')
        axs[1].set_ylabel('w1 (Hz)')
        axs[1].grid('on')
        axs[1].legend(loc='upper center')
        axs[1].set_ylim(np.min(w1_range_Hz), np.max(w1_range_Hz))

        # dealing with axes Hz/V
        ax2 = axs[1].twinx()
        ax2.set_ylabel('Pulse voltage (V)')
        ax2.set_ylim(np.min(w1_range_V), np.max(w1_range_V))

        log.debug("HSn pulse voltage is set to %.2fV" % self.pulse_rfc_voltage)
        log.debug("HSn pulse flip angle is set to %.2fdeg" % self.pulse_rfc_flipangle)
        log.debug("reference voltage is %.2fV" % self.vref)
        log.debug("which is equivalent to a w1max of %.2fHz" % w1_rfc_pulse_hz)
        log.debug("the optimization process gave however something different...")
        log.debug("the optimal w1max found is %.2fHz" % w1_range_Hz_optim)
        log.debug("which corresponds to a voltage of %.2fV" % self.pulse_rfc_optim_power_voltage)
        log.debug("or a flip angle of %.2fdeg" % self.pulse_rfc_optim_power_flipangle)

        return(w1_range_Hz_optim)

    def _run_sequence(self, metabolite):
        """
        Simulate the NMR signal acquired using a sLASER sequence using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_basis_set entry for one single metabolite

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        log.debug("acquiring sLASER sequence: (90)-(180)-(180)-(180)-(180)...")

        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite)

        if(self.pulse_rfc_real_shape_enable):
            # convert pulse amplitude currently normalized (0-1) to Hz (compatible with GAMMA)
            pulse_rfc_amp_np = self.pulses_rfc_shape_amp_norm * self._pulse_rfc_w1max
            pulse_rfc_phi_np = self.pulses_rfc_shape_phi_deg

            # prepare pyGAMMA RF pulse vectors (from Vespa)
            pg_pulse_step_s = self.pulse_rfc_duration / (1000.0 * len(pulse_rfc_amp_np))
            pg_pulse_rfc = pg.row_vector(len(pulse_rfc_amp_np))
            pg_pulse_rfc_time = pg.row_vector(len(pulse_rfc_amp_np))
            for j, val in enumerate(zip(pulse_rfc_amp_np, pulse_rfc_phi_np)):
                pg_pulse_rfc.put(pg.complex(val[0], val[1]), j)
                pg_pulse_rfc_time.put(pg.complex(pg_pulse_step_s, 0), j)

            # create the pulse waveform and composite pulse objects  (from Vespa)
            pwf_hsn = pg.PulWaveform(pg_pulse_rfc, pg_pulse_rfc_time, "HSn LASER pulse")
            pulc180 = pg.PulComposite(pwf_hsn, sys, self.nuclei)
            Upulc180 = pulc180.GetUsum(-1)

        # amplitude normalization
        afactor = self.scaling_factor / sys.HS()

        # coupling
        H = pg.Hcs(sys) + pg.HJ(sys)

        # detection stuff
        D = pg.Fm(sys, self.nuclei)
        dt = np.double(1 / self.fs)

        # timing implementation close to what the CMRR's sLASER is doing
        # those delays are between pulses' centers
        # 90-a-180-b-c-180-d-e-180-f-g-180-h-FID
        # all in ms
        a = self.pulse_exc_duration / 2.0 + self.pulse_rfc_duration / 2.0 + self.spoiler_duration + 0.4e-3  # timing found experimentally on the 7T
        bc = a + self.pulse_rfc_duration / 2.0 + self.spoiler_duration  # SE timing
        de = self.pulse_rfc_duration / 2.0 + self.spoiler_duration + self.spoiler_duration + self.pulse_rfc_duration / 2.0  # assuming shortest
        f = self.pulse_rfc_duration / 2.0 + self.spoiler_duration
        gh_min = self.pulse_rfc_duration + 2.0 * self.spoiler_duration + 0.0  # according to measurements, TE padding happens on the last SE timing
        # let's calculate min TE to find TE filling value
        te_min = a + bc + de + f + gh_min
        te_fill = self.te - te_min
        # check for impossible TE
        if(te_fill < 0.0):
            log.error("your echo time is too short!")
        gh = self.pulse_rfc_duration + 2.0 * self.spoiler_duration + te_fill
        # final SE timing
        g = gh / 2.0
        h = gh / 2.0
        fg = f + g

        # final delay list: center to center of pulses
        self._subdelay_c2c_ms_list = np.array([a, bc, de, fg, h])

        # now, calculate evolution delays knowing the c2c delays and the user parameters
        # 4 cases possible
        # those delays are the actual time when spins will evolve
        # 90-a-180-b-c-180-d-e-180-f-g-180-h-FID
        if(self.pulse_rfc_real_shape_enable and self.allow_evolution_during_hard_pulses):
            # we are using real shaped rfc pulses and a hard exc pulse
            # we want to keep time evolution during hard pulsing
            # evolution already happens during shaped pulses
            a += -self.pulse_rfc_duration / 2.0
            bc += -self.pulse_rfc_duration
            de += -self.pulse_rfc_duration
            fg += -self.pulse_rfc_duration
            h += -self.pulse_rfc_duration / 2.0
        elif(self.pulse_rfc_real_shape_enable and not self.allow_evolution_during_hard_pulses):
            # we are using real shaped rfc pulses and a hard exc pulse
            # we do not want to keep time evolution during hard pulsing
            # evolution already happens during shaped pulses
            a += -self.pulse_rfc_duration / 2.0 - self.pulse_exc_duration / 2.0
            bc += -self.pulse_rfc_duration
            de += -self.pulse_rfc_duration
            fg += -self.pulse_rfc_duration
            h += -self.pulse_rfc_duration / 2.0
        elif(not self.pulse_rfc_real_shape_enable and self.allow_evolution_during_hard_pulses):
            # we are not using real shaped rfc pulses, only a hard exc and rfc pulses
            # we want to keep time evolution during hard pulsing
            pass
        elif(not self.pulse_rfc_real_shape_enable and not self.allow_evolution_during_hard_pulses):
            # we are not using real shaped rfc pulses, only a hard exc and rfc pulses
            # we do not want to keep time evolution during hard pulsing
            a += -self.pulse_rfc_duration / 2.0 - self.pulse_exc_duration / 2.0
            bc += -self.pulse_rfc_duration
            de += -self.pulse_rfc_duration
            fg += -self.pulse_rfc_duration
            h += -self.pulse_rfc_duration / 2.0

        # final delay list: evolution delays
        self._subdelay_evol_ms_list = np.array([a, bc, de, fg, h])

        # final tcp delay list (for CP effects)
        self.tcp_ms_list = self._subdelay_c2c_ms_list.copy()

        # run the sequence
        sigma0 = pg.sigma_eq(sys)
        evol_delays_s = self._subdelay_evol_ms_list / 1000.0  # seconds
        te_real_ms = 0.0

        # excitation: hard 90 pulse to simulate asymmetric sinc pulse
        sigma1 = pg.Iypuls(sys, sigma0, self.nuclei, 90.0)
        if(not self.allow_evolution_during_hard_pulses):
            te_real_ms += self.pulse_exc_duration / 2.0
        # evolution
        sigma0 = pg.evolve(sigma1, pg.prop(H, evol_delays_s[0]))
        te_real_ms += evol_delays_s[0] * 1000.0

        # LASER: 2 x (pair of 180 HSn pulses)
        for d in evol_delays_s[1:]:
            if(self.pulse_rfc_real_shape_enable):
                sigma1 = Upulc180.evolve(sigma0)
                te_real_ms += self.pulse_rfc_duration
            else:
                sigma1 = pg.Iypuls(sys, sigma0, self.nuclei, 180.0)
                if(not self.allow_evolution_during_hard_pulses):
                    te_real_ms += self.pulse_rfc_duration

            # evolution
            sigma0 = pg.evolve(sigma1, pg.prop(H, d))
            te_real_ms += d * 1000.0

        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        log.debug("done. (real TE=%.2fms)" % te_real_ms)

        return(s_MRSData2)


class mrs_seq_eja_svs_press(mrs_seq_press):
    """A class that represents the PRESS sequence from CMRR."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, te1=np.nan, te2=np.nan):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        te1 : float
            First part of TE (ms)
        te2 : float
            Second part of TE (ms)
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor, te1, te2)
        # name of sequence
        self.name = "eja_svs_press"


class mrs_seq_eja_svs_steam(mrs_seq_steam):
    """A class that represents the STEAM sequence from CMRR."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor, tm)
        # name of sequence
        self.name = "eja_svs_steam"


class mrs_seq_fid(mrs_sequence):
    """A class that represents the pulse-acquire sequence (fid), which is actually a clone of the super class."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0):
        """
        Initialize a virtual FID sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        scaling_factor : float
            Scaling FID intensity factor
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor)
        # name of sequence
        self.name = "fid"


class mrs_seq_svs_se(mrs_seq_press):
    """A class that represents the PRESS sequence from SIEMENS."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, te1=np.nan, te2=np.nan):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        te1 : float
            First part of TE (ms)
        te2 : float
            Second part of TE (ms)
        """
        super().__init__(te, tr, na, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor, te1, te2)
        # name of sequence
        self.name = "svs_se"


class mrs_seq_svs_st(mrs_seq_steam):
    """A class that represents the STEAM sequence from SIEMENS."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        eff_acquisition_time : float
            Effective acquisition time (s)
        scaling_factor : float
            Scaling FID intensity factor
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor, tm)
        # name of sequence
        self.name = "svs_st"


class mrs_seq_svs_st_vapor_643(mrs_seq_steam):
    """A class that represents the STEAM sequence from SIEMENS WiP 643."""

    def __init__(self, te, tr=3500.0, na=128, ds=4, nuclei="1H", npts=4096 * 4, voxel_size=[10, 10, 10], fs=5000.0, f0=297.2062580, vref=250.0, shims=[], timestamp=np.nan, gating_mode=gating_signal_source.NO_GATING, eff_acquisition_time=np.nan, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        na : int
            Number of averages/excitations
        ds : int
            Number of summy scans
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        voxel_size : list
            Dimensions of voxel (mm)
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        vref : float
            Reference voltage (V)
        shims : list of floats
            List of shim voltages in volts
        timestamp : float
            Timestamp in ms
        gating_mode : gating_signal_source
            Acquisition triggering mode
        scaling_factor : float
            Scaling FID intensity factor
        eff_acquisition_time : float
            Effective acquisition time (s)
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, na, ds, nuclei, npts, voxel_size, fs, f0, vref, shims, timestamp, gating_mode, eff_acquisition_time, scaling_factor, tm)
        # name of sequence
        self.name = "svs_st_vapor_643"


class metabolite_basis_set(dict):
    """The metabolite_basis_set class is a big dictionnary of metabolites with their respective chemical shift and J coupling information."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self):
        """Construct a metabolite_basis_set object."""
        super(metabolite_basis_set, self).__init__()
        # xls file that contains metabolites properties (you should not change that)
        self._database_xls_file = default_paths.DEFAULT_META_DB_FILE
        # xls file that contains your metabolite basis set
        self.basis_set_xls_file = default_paths.DEFAULT_META_BASIS_SET_FILE
        # include peaks only in this range of chemical shifts
        self.ppm_range = [0, 6]
        # include only non-coupled peaks
        self.non_coupled_only = False
        # set all metabolites to have only one peak (first chemical shift) and one proton
        # one proton only for each metabolite
        self.one_proton_mode = False
        # to know if the object is initialized
        self._ready = False

        # freeze
        self.__isfrozen = True

    def __eq_dict(self, d1, d2):
        """
        Compare and return equality check for nested dictionnaries d1 and d2.

        Parameters
        ----------
        d1 : dict
            Dictionnary 1
        d2 : dict
            Dictionnary 2

        Returns
        -------
        r : bool
            Self if equal to other
        """
        # check keys
        if(d1.keys() != d2.keys()):
            return(False)
        # check key values
        for k in list(d1.keys()):
            # check key values type
            if(type(d1[k]) == dict and type(d2[k]) == dict):
                r = self.__eq_dict(d1[k], d2[k])
            elif(type(d1[k]) == np.ndarray and type(d2[k]) == np.ndarray):
                r = (d1[k] == d2[k]).all()
            else:
                r = (d1[k] == d2[k])

            if(not r):
                return(False)

        return(True)

    def __eq__(self, other):
        """
        Overload of "equal to" operator so that we can compare two metabolite basis set (usefull when looking for a adequat sequence/basis set).

        Returns
        -------
        r : bool
            Self if equal to other
        """
        # check attributes
        if(self.__dict__ != other.__dict__):
            return(False)
        # check dict content
        r = self.__eq_dict(self, other)
        return(r)

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

    def _read_xls_file(self):
        """Read the xls file specified by "self.xls_file" and stores all the information (chemical shifts, nuclei, J couplings)."""
        # TODO: need to work with gitable xls format, FODS!
        # remove annoying warning
        warnings.simplefilter(action='ignore', category=FutureWarning)

        # first, read xls metabolite basis set file and parse it
        log.info("reading metabolite basis set from XLS file...")
        book_db = open_workbook(self._database_xls_file)
        book = open_workbook(self.basis_set_xls_file)

        # get sheets from database file
        sheet_names = book_db.sheet_by_name('Metabolites')
        sheet_ppm = book_db.sheet_by_name('PPM')
        sheet_iso = book_db.sheet_by_name('Nuclei')

        # get sheets from basis set file
        sheet_groups = book.sheet_by_name('Metabolites')
        sheet_ref = book.sheet_by_name('Reference_metabolite')
        sheet_MMs = book.sheet_by_name('Macromolecules')
        sheet_refMM = book.sheet_by_name('Reference_macromolecule')
        sheet_cmin = book.sheet_by_name('Concentrations_min')
        sheet_cmax = book.sheet_by_name('Concentrations_max')
        sheet_T1s = book.sheet_by_name('T1s')
        sheet_T2s = book.sheet_by_name('T2s')

        # metabolite used as a reference, a robust peak
        ref_meta_name = sheet_ref.col_values(0)[0]

        # MMs list and reference MM
        MM_list = sheet_MMs.col_values(0)
        if(len(sheet_refMM.col_values(0)) > 0):
            ref_MM_name = sheet_refMM.col_values(0)[0]
        else:
            ref_MM_name = ""

        # other sheets
        c_min = np.array(sheet_cmin.col_values(0))
        c_max = np.array(sheet_cmax.col_values(0))
        T1s = np.array(sheet_T1s.col_values(0))
        T2s = np.array(sheet_T2s.col_values(0))

        # browse metabolite groups and MMs
        for this_sheet in [sheet_groups, sheet_MMs]:
            # metagroup list
            for this_metagroup_name in this_sheet.col_values(0):
                # index
                this_metagroup_ind = this_sheet.col_values(0).index(this_metagroup_name)

                # prepare group entry
                this_metagroup_entry = {"metabolites": {}}

                # look for the members of the group
                for this_meta_name in this_sheet.row_values(this_metagroup_ind):
                    # is this metagroup name is a metabolite or not ?
                    if(this_meta_name in sheet_names.col_values(0)):
                        # its index
                        this_meta_ind = sheet_names.col_values(0).index(this_meta_name)

                        # find info on this meta
                        # ppm
                        this_meta_ppm_line = np.array(sheet_ppm.row_values(this_meta_ind))
                        this_meta_ppm_line[this_meta_ppm_line == ''] = np.nan
                        this_ppm_list = this_meta_ppm_line.astype(np.float64)
                        this_meta_mask = ~np.isnan(this_ppm_list)
                        this_ppm_list = this_ppm_list[this_meta_mask]

                        # ppm filtering: hope that works
                        this_ppm_list[(this_ppm_list < self.ppm_range[0]) | (this_ppm_list > self.ppm_range[1])] = 100.0

                        # iso
                        this_meta_iso_line = np.array(sheet_iso.row_values(this_meta_ind))
                        this_meta_iso_line[this_meta_iso_line == ''] = np.nan
                        this_iso_list = this_meta_iso_line.astype(np.float64)
                        this_iso_list = this_iso_list[this_meta_mask]

                        # J couplings
                        # check if there is a sheet for this metabolite
                        if(this_meta_name in book_db.sheet_names()):
                            this_sheet_Jcouplings = book_db.sheet_by_name(this_meta_name)
                            this_j_list = np.full([len(this_ppm_list), len(this_ppm_list)], np.nan)
                            for irowJ in range(len(this_ppm_list)):
                                # ppm
                                this_meta_j_line = np.array(this_sheet_Jcouplings.row_values(irowJ))
                                this_meta_j_line[this_meta_j_line == ''] = np.nan
                                this_meta_j_line = this_meta_j_line.astype(np.float64)
                                this_meta_j_line = this_meta_j_line[this_meta_mask]
                                this_j_list[irowJ][:] = this_meta_j_line
                        else:
                            # if no sheet, put all J-couplings to zeroes
                            this_j_list = np.zeros([len(this_ppm_list), len(this_ppm_list)])

                        # filter coupled and non-coupled peaks
                        if(self.non_coupled_only):
                            # sum up one axis
                            this_j_list_summed = np.sum(np.abs(this_j_list), axis=0)
                            this_j_list_summed_mask = (this_j_list_summed != 0.0)
                            this_ppm_list[this_j_list_summed_mask] = 100.0

                        # append metabolite entry for this metabolite group
                        this_metagroup_entry["metabolites"][this_meta_name] = {"ppm": this_ppm_list, "iso": this_iso_list, "J": this_j_list}

                # add extra infos to metabolite group
                # is it a MM?
                this_MM_bool = (this_metagroup_name in MM_list)
                this_metagroup_entry["Macromecule"] = this_MM_bool
                # omg, is it the reference MM?
                this_ref_MM_bool = (this_metagroup_name == ref_MM_name)
                this_metagroup_entry["Reference macromolecule"] = this_ref_MM_bool
                # oh maybe, is it the reference metabolite?
                this_ref_meta_bool = (this_metagroup_name == ref_meta_name)
                this_metagroup_entry["Reference metabolite"] = this_ref_meta_bool
                # min/max concentration
                this_metagroup_entry["Concentration min"] = c_min[this_metagroup_ind]
                this_metagroup_entry["Concentration max"] = c_max[this_metagroup_ind]
                # T1s/T2s
                this_metagroup_entry["T1"] = T1s[this_metagroup_ind]
                this_metagroup_entry["T2"] = T2s[this_metagroup_ind]

                # and add the group to the dict
                self[this_metagroup_name] = this_metagroup_entry

    def _write_header_file(self):
        """Interesting method here... It generates a python .py and writes very usefull aliases to access quickly to a specific metabolite or parameter. Since metabolite indexes depend on the metabolite database, this python header file is regenerated each time the metabolite basis set is initialized."""
        log.debug("generating metabolite and parameter aliases...")

        # find location of package
        # this file needs to be specifically in the package folder!
        # because we will later import it
        # (I know it is a bit ugly but so practical)
        pkg_folder = str(pathlib.Path(__file__).parent.absolute())

        with open(pkg_folder + "/aliases.py", 'w') as f:
            f.write("#!/usr/bin/env python3")
            f.write("# -*- coding: utf-8 -*-\n")
            f.write('\n')
            f.write("# This file is automatically generated on the fly by sim.py! Do not try to modify it please.\n")
            f.write("\n")

            # metabolite aliases
            metagroup_names = list(self.keys())
            for k, this_meta_name in enumerate(metagroup_names):
                f.write("m_" + this_meta_name + " = " + str(k) + "\n")
            f.write("\n")

            # parameter aliases
            f.write("p_cm = 0\n")
            f.write("p_dd = 1\n")
            f.write("p_df = 2\n")
            f.write("p_dp = 3\n")
            f.write("\n")

            # all metabolite alias
            ind_ref_meta = None
            n_meta = 0
            f.write("m_All_MBs = [")
            for k, this_meta_name in enumerate(metagroup_names):
                if(not self[this_meta_name]["Macromecule"]):
                    f.write(str(k) + ", ")
                    n_meta += 1
                if(self[this_meta_name]["Reference metabolite"]):
                    ind_ref_meta = k
            f.write("]\n")

            # ref metabolite alias
            if(ind_ref_meta is None):
                log.error("weird, could not find reference metabolite!?")
            else:
                f.write("m_Ref_MB = %d\n" % ind_ref_meta)
            f.write("\n")

            # all MMs alias
            ind_ref_MM = None
            n_MM = 0
            f.write("m_All_MMs = [")
            for k, this_meta_name in enumerate(metagroup_names):
                if(self[this_meta_name]["Macromecule"]):
                    f.write(str(k) + ", ")
                    n_MM += 1
                if(self[this_meta_name]["Reference macromolecule"]):
                    ind_ref_MM = k
            f.write("]\n")

            # ref MM alias
            if(ind_ref_MM is None):
                f.write("m_Ref_MM = None\n")
            else:
                f.write("m_Ref_MM = %d\n" % ind_ref_MM)
            f.write("\n")

            f.write("n_All = " + str(n_meta + n_MM) + "\n")
            f.write("n_MBs = " + str(n_meta) + "\n")
            f.write("n_MMs = " + str(n_MM) + "\n")
            f.write("\n")

    def print(self):
        """Print the metabolite database with all the information (chemical shifts, nuclei, J couplings)."""
        cell_nchar = 8

        # and now browse though the database and display everything
        for this_metagroup_key, this_metagroup_entry in self.items():
            for this_meta_key, this_meta_entry in self[this_metagroup_key]["metabolites"].items():
                print("")

                print("> metabolite [", end="")
                cprint(this_metagroup_key, 'green', attrs=['bold'], end="")
                print("/", end="")
                cprint(this_meta_key, 'green', attrs=['bold'], end="")
                print("]")

                # first nuclei line
                print("nuclei".ljust(cell_nchar) + "".ljust(cell_nchar), end="")
                for this_iso in this_meta_entry["iso"]:
                    print(("(%2d)" % this_iso).ljust(cell_nchar), end="")
                print()

                print("".ljust(cell_nchar) + "ppm".ljust(cell_nchar), end="")
                for this_ppm in this_meta_entry["ppm"]:
                    print(("%.1f" % this_ppm).ljust(cell_nchar), end="")
                print()

                for (this_iso, this_ppm, irowJ) in zip(this_meta_entry["iso"], this_meta_entry["ppm"], range(len(this_meta_entry["ppm"]))):
                    # recall nuclei and chemical shift for 1st and 2nd col
                    print(("(%2d)" % this_iso).ljust(cell_nchar), end="")
                    print(("%.1f" % this_ppm).ljust(cell_nchar), end="")

                    for icolJ in range(len(this_meta_entry["J"])):
                        print(("%.1f" % this_meta_entry["J"][irowJ][icolJ]).ljust(cell_nchar), end="")
                    print("", flush=True)

    def initialize(self):
        """Initialize metabolite database: run the two previous methods."""
        log.info("initializing metabolite database...")
        self.clear()
        self._read_xls_file()
        self._write_header_file()

        # single proton AMARES style?
        if(self.one_proton_mode):
            for m in list(self.keys()):
                for sm in list(self[m]["metabolites"].keys()):
                    self[m]["metabolites"][sm]["ppm"] = np.array([self[m]["metabolites"][sm]["ppm"][0]])
                    self[m]["metabolites"][sm]["iso"] = np.array([self[m]["metabolites"][sm]["iso"][0]])
                    self[m]["metabolites"][sm]["J"] = np.array([[self[m]["metabolites"][sm]["J"][0][0]]])

        self._ready = True
