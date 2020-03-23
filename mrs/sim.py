#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Three classes' definition in here.

    * a params class which stores and manipulate the parameters of our MRS fitting/simulation model
    * a metabolite class which stores and can compute a MRS modeled signal for a single metabolite, based on the pyGAMMA library (for python 3!) using a specific MR sequence described by pulse flip angles and delays
    * a metabolite_group class which contains several metabolites
    * a metabolite_db class which contains a whole database of metabolites chemical shifts, J-couplings, nucleis and computed signals. Usefull to simulate all kind of MRS data for various metabolites, concentrations acquired with various sequences...

@author: Tangi Roussel
"""

try:
    import pygamma as pg
    GAMMA_LIB_LOADED = True
except ImportError:
    GAMMA_LIB_LOADED = False

import suspect
import numpy as np
import math as ma
import mrs.reco as reco
import mrs.metabase as xxx
import matplotlib.pylab as plt
import os
import warnings
from xlrd import open_workbook
from termcolor import cprint
import matplotlib._color_data as mcd

import pdb

GAMMA_DICT = {"1H": 42.577478518, "13C": 10.7084, "19F": 40.052, "23Na": 11.262, "31P": 17.235}
PARS_LONG_NAMES = ["Concentration [mmol/kg]", "Linewidth factor [Hz]", "Frequency shift [Hz]", "Phase [rd]"]
PARS_SHORT_NAMES = ["cm", "dd", "df", "dp"]


class params(np.ndarray):
    """A class that stores the parameters used to modelize a MR spectrum during simulation or fit."""

    def __init__(self, metadb):
        """
        Initialize a params object.

        Parameters
        ----------
        metadb: metabolite_db
            A metabolite_db object to which this params object is linked to
        """
        super().__init__()
        # those parameters are related to a metabolite database
        self._metadb = metadb
        # the link-lock vector used to control the model
        self._linklock = np.zeros(self.shape)
        # the ratio vector
        self._ratios = np.ones(self.shape)

    def __new__(cls, metadb):
        """
        Construct a params object that inherits of numpy array's class. This class is used to deal with metabolite parameters.

        Parameters
        ----------
        metadb: metabolite_db
            A metabolite_db object to which this params object is linked to

        Returns
        -------
        obj : params numpy array [n,4]
            Resulting constructed params object
        """
        obj = super(params, cls).__new__(cls, [len(metadb), 4])
        obj[:] = 0.0
        return(obj)

    def __array_finalize__(self, obj):
        """
        Overload of special numpy array method called when playing around with stuff relative to object copy etc...

        Parameters
        ----------
        obj : params numpy array [n,4]
        """
        self._metadb = getattr(obj, 'metadb', None)
        self._linklock = getattr(obj, 'linklock', None)
        self._ratios = getattr(obj, 'ratios', None)

    @property
    def metadb(self):
        """Property method for metadb."""
        return(self._metadb)

    @property
    def linklock(self):
        """Property method for linklock."""
        return(self._linklock)

    @property
    def ratios(self):
        """Property method for ratios."""
        return(self._ratios)

    def get_meta_names(self):
        """
        Return list of metabolite names controlled by this params object.

        Returns
        -------
        list(self._metadb.keys()) : list
            List of metabolite names
        """
        return(list(self._metadb.keys()))

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
        for l in LL_list:
            # count number of master
            tmp = (self.linklock == -l)
            n_masters = np.sum(tmp[:])

            # count number of slaves
            tmp = (self.linklock == +l)
            n_slaves = np.sum(tmp[:])

            # master ratio values should be all equal to one
            tmp = (self.linklock == -l)
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
        for l in LL_list:
            p[p.linklock == +l] = p[p.linklock == -l] * p.ratios[p.linklock == +l]

        return(p)

    def set_T2_weighting(self, te, extra_info_key="T2s_human_brain_7T"):
        """
        Recalculate metabolite amplitude parameters by applying a T2 relaxation, usefull for simulations.

        Parameters
        ----------
        te : float
            Echo time in (ms)
        extra_info_key : dictionnary key / string
            Tab name in the metabase XLS file containing T2 values

        Returns
        -------
        self.copy() : params object
            Full array of parameters corrected in T2
        """
        # init: we are working on a copy
        p = self.copy()

        # browse though the database and find T2s
        params_T2s = []
        for this_metagroup_key, this_metagroup_entry in self._metadb.items():
            params_T2s.append(this_metagroup_entry[extra_info_key])

        # convert to np
        params_T2s = np.array(params_T2s)

        # apply T2w for given TE
        p[:, xxx.p_cm] = p[:, xxx.p_cm] * np.exp(-te / params_T2s)
        return(p)

    def correct_T2s(self, te, extra_info_key="T2s_human_brain_7T"):
        """
        Correct the concentration values of a parameter array depending on the TE and the common values of T2s for each metabolite.

        Parameters
        ----------
        te : float
            Echo time in (ms)
        extra_info_key : dictionnary key / string
            Tab name in the metabase XLS file containing T2 values

        Returns
        -------
        self.copy() : params object
            Full array of parameters corrected in T2
        """
        # init: we are working on a copy
        p = self.copy()

        # browse though the database and find T2s
        params_T2s = []
        for this_metagroup_key, this_metagroup_entry in self._metadb.items():
            params_T2s.append(this_metagroup_entry[extra_info_key])

        # convert to np
        params_T2s = np.array(params_T2s)

        # finding real concentration values at TE=0ms
        p[:, xxx.p_cm] = p[:, xxx.p_cm] / np.exp(-te / params_T2s)
        return(p)

    def correct_T1s(self, tr, extra_info_key="T1s_human_brain_7T"):
        """
        Correct the concentration values of a parameter array depending on the TR and the common values of T1s for each metabolite.

        Parameters
        ----------
        tr : float
            Repetition time in (ms)
        metabase_extra_infos_key : dictionnary key / string
            Tab name in the metabase XLS file containing T1 values

        Returns
        -------
        self.copy() : params object
            Full array of parameters corrected in T1
        """
        # init: we are working on a copy
        p = self.copy()

        # browse though the database and find T1s
        params_T1s = []
        for this_metagroup_key, this_metagroup_entry in self._metadb.items():
            params_T1s.append(this_metagroup_entry[extra_info_key])

        # convert to np
        params_T1s = np.array(params_T1s)

        # finding real concentration values at TE=0ms
        p[:, xxx.p_cm] = p[:, xxx.p_cm] / (1 - np.exp(-tr / params_T1s))
        return(p)

    def _correct_normalize_concentrations(self, multiplication_factor):
        """
        Multiply the concentrations by a factor.

        Parameters
        ----------
        multiplication_factor : float
            Multiplication factor

        Returns
        -------
        self.copy() : params object
            Full array of parameters
        """
        # init: we are working on a copy
        p = self.copy()
        p[:, xxx.p_cm] = p[:, xxx.p_cm] * multiplication_factor
        return(p)

    def correct_absolute(self, params_water, water_concentration=55000.0):
        """
        Calculate the absolute metabolic concentration values knowing the relative water concentration value and its assumed concentration.

        Parameters
        ----------
        params_water : params object
            Array of parameters obtained from the fit a non water suppressed MRS data (reference scan)
        water_concentration : float
            Assumed water concentration used to calculate absolute concentration estimates (mmol/kg)

        Returns
        -------
        self.copy() : params object
            Full array of parameters
        """
        return(self._correct_normalize_concentrations(water_concentration / params_water[xxx.m_Water, xxx.p_cm]))

    def correct_relative(self, params_water, water_concentration=55000.0):
        """
        Calculate the relative metabolic concentration values knowing the relative water concentration value and its assumed concentration.

        Parameters
        ----------
        params_water : params object
            Array of parameters obtained from the fit a non water suppressed MRS data (reference scan)
        water_concentration : float
            Assumed water concentration used to calculate absolute concentration estimates (mmol/kg)

        Returns
        -------
        self.copy() : params object
            Full array of parameters
        """
        return(self._correct_normalize_concentrations(params_water[xxx.m_Water, xxx.p_cm] / water_concentration))

    def set_default_min(self):
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
        self.linklock[:] = np.tile([0, 2, 3, 4], (xxx.n_ALL, 1))
        self.linklock[xxx.m_NAA_CH3, :] = [0, -2, -3, -4]
        self.linklock[xxx.m_Water, :] = 0

        # no MMs
        self.linklock[xxx.m_All_MMs, :] = 1

        return(self.copy())

    def set_default_max(self):
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
        self.linklock[:] = np.tile([0, 2, 3, 4], (xxx.n_ALL, 1))
        self.linklock[xxx.m_NAA_CH3, :] = [0, -2, -3, -4]
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
        self.set_default_min()
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
        self.set_default_max()
        # all concentrations to zero
        self[:, xxx.p_cm] = 0.0

        # water max values
        self[xxx.m_All_MBs, xxx.p_df] = +20.0
        self[xxx.m_Water, 0] = 100000.0

        # lock everything except water
        self.linklock[:] = 1
        self.linklock[xxx.m_Water, :] = 0
        return(self.copy())

    def set_default_human_brain_min(self):
        """
        Initialize the params object to the minimum in vivo human values, no macromolecules.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        self.set_default_min()

        # browse though the database and find human brain min values from literature
        cmin = []
        for this_metagroup_key, this_metagroup_entry in self._metadb.items():
            cmin.append(this_metagroup_entry["Concentrations_human_brain_min"])

        # convert to np
        cmin = np.array(cmin)

        self[xxx.m_All_MBs, xxx.p_cm] = cmin[xxx.m_All_MBs]

        return(self.copy())

    def set_default_human_brain_max(self):
        """
        Initialize the params object to the maximum in vivo human values, no macromolecules.

        Returns
        -------
        self.copy() : params object
            Copy of this current object with the applied modification
        """
        self.set_default_max()

        # browse though the database and find human brain min values from literature
        cmax = []
        for this_metagroup_key, this_metagroup_entry in self._metadb.items():
            cmax.append(this_metagroup_entry["Concentrations_human_brain_max"])

        # convert to np
        cmax = np.array(cmax)

        self[xxx.m_All_MBs, xxx.p_cm] = cmax[xxx.m_All_MBs]

        return(self.copy())

    def get_MBs_human_above(self, cm_threshold):
        """
        Return indexes for metabolites with human brain average concentrations above a threshold.

        Parameters
        ----------
        cm_threshold : float
            Concentration threshold

        Returns
        -------
        im : list of integers
            Metabolite indexes
        """
        cm_mean = (self.set_default_human_brain_min() + self.set_default_human_brain_max()) / 2.0
        cm_mean = cm_mean[:, 0]
        metabolites_inds = np.where(cm_mean > cm_threshold)
        metabolites_inds = metabolites_inds[0]

        # display
        print(">> mrs.sim.params.get_MBs_human_above: found %d metabolites above %f mmol/kg..." % (len(metabolites_inds), cm_threshold))
        print(" > ", end="", flush=True)
        meta_names = self.get_meta_names()
        for im in metabolites_inds:
            print("%s " % str(meta_names[im]), end="", flush=True)
        print()

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
        self.linklock[xxx.m_All_MMs, :] = np.tile(
            [0, 200, 300, 400], (xxx.n_MMs, 1))
        self.linklock[xxx.m_MM1, :] = [0, -200, -300, -400]

        return(self.copy())

    def add_macromolecules_max(self):
        """
        Initialize the params object to the maximum in vivo human values.

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
        self.linklock[xxx.m_All_MMs, :] = np.tile(
            [0, 200, 300, 400], (xxx.n_MMs, 1))
        self.linklock[xxx.m_MM1, :] = [0, -200, -300, -400]
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
        print(">> mrs.sim.params.display_parameters: parameters...")
        if(bMM):
            n = xxx.n_ALL
        else:
            n = xxx.n_MBs

        meta_names = self.get_meta_names()
        LLtransTab = str.maketrans("-+0123456789", "⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹")
        LLcolors = ['white', 'grey', 'yellow', 'blue', 'magenta', 'cyan', 'white', 'grey', 'yellow', 'blue', 'magenta', 'cyan', 'green']

        print("[#] Metabolite name          \t [cm] \t\t[dd] \t\t[df] \t\t [dp]")
        for k in range(n):
            print("#%2d %-20s" %
                  (k, meta_names[k]), end="", flush=True)
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
                    cprint("\t(%4.1f)" % self[k, kp], thisLL_color, attrs=['bold'], end="", flush=True)
                    cprint(thisLL_str, thisLL_color, attrs=['bold'], end="", flush=True)
                else:
                    print("\t(%4.1f)" % self[k, kp], end="", flush=True)

            print()
        print()


class mrs_sequence:
    """A class that stores a sequence and all its parameters used for simulation. This is a generic sequence class that you need to overload. By default, the simulated sequence is a simple pulse-acquire NMR experiment."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0):
        """
        Initialize the sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        """
        # all sequences have those parameters
        # name
        self._name = "fid"
        # echo time (ms)
        self.te = te
        # repetition time (ms)
        self.tr = tr
        # which nuclei we are pulsing on/looking at (examples: '1H', '31P')
        self.nuclei = nuclei
        # number of acquired time points (int)
        self.npts = npts
        # sampling frequency (Hz)
        self.fs = fs
        # larmor frequency of water (MHz)
        self.f0 = f0
        # kind of receiver gain
        self.scaling_factor = scaling_factor
        # ppm shift (ppm)
        self.ppm_water = 4.7
        # some 0th phase to add? (rd)
        self.additional_phi0 = 0.0

        # metabolite_db object
        self._meta_db = metabolite_db()

        # try to simplify spin systems when possible to speed up simulations
        self.allow_spin_system_simplification = True
        # NMR simulation option: when hard zero-duration RF pulse are employed, should we take into account the duration of the real pulses in the evolution times or not? Experimentally, looks like yes.
        self.allow_evolution_during_hard_pulses = True

        # pre-calculated stuff
        # 'metabase': set of numerically computed metabolite FID signals
        self._meta_signals = None
        # time vector
        self._t = []
        # initialized or not
        self._ready = False

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
    def name(self):
        """
        Property get function for name.

        Returns
        -------
        self._name : string
            name of sequence
        """
        return(self._name)

    @property
    def meta_db(self):
        """
        Property get function for meta_db.

        Returns
        -------
        self._meta_db : metabolite_db object
            Metabolite database to use for simulation
        """
        return(self._meta_db)

    def _init_pulses(self):
        """Virtual method which initialize RF pulse waveforms if any."""

    def _prepare_spin_system(self, metabolite, bSilent=False):
        """
        Return pyGAMMA spin system for a given metabolite, knowing all its properties. Simplify the system in simple cases like singulets (one single chemical shift, no J-couplings).

        Parameters
        ----------
        metabolite : dict
            metabolite_db entry for one single metabolite
        beSilent : boolean
            No output in console (True)

        Returns
        -------
        sys : pyGAMMA system object
            Spin system object used for simulation
        scaling_factor : float
            Scaling factor after system simplification
        """
        # hello
        if(not bSilent):
            cprint(">>>> mrs.sim.mrs_sequence._prepare_spin_system:", 'green', attrs=['bold'])
            print("  > preparing spin system", end="", flush=True)

        # init
        scaling_factor = 1.0

        # extract metabolite properties needed to create spin system
        ppm = metabolite["ppm"]
        iso = metabolite["iso"]
        j = metabolite["J"]

        # check if we can simplify
        if(self.allow_spin_system_simplification and len(ppm) > 1 and not np.any(j != 0.0) and len(np.unique(ppm)) == 1):
            # we have a non-coupled singulet with N spins here
            # let's simplify to one single spin + amplification factor
            scaling_factor = float(len(ppm))
            ppm = np.array([ppm[0]])
            iso = np.array([iso[0]])
            j = np.array([j[0, 0]])

            if(not bSilent):
                print(", simplified it.", end="", flush=True)

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
                raise Exception("  > mrs.sim.mrs_sequence._run_sequence: " + str(iso[i_spin]) + ", that is weird nuclei!")

            # set the ppm
            sys.PPM(i_spin, ppm[i_spin])

            # set the couplings
            for icolJ in range(i_spin + 1, len(ppm)):
                sys.J(i_spin, icolJ, j[i_spin, icolJ])

                if(not bSilent):
                    print(".", end="", flush=True)

            if(not bSilent):
                print(".", end="", flush=True)

        # water shift
        sys.offsetShifts(self.ppm_water * self.f0, self.nuclei)

        if(not bSilent):
            print("done.")

        return(sys, scaling_factor)

    def _run_sequence(self, metabolite):
        """
        Simulate the NMR signal acquired using a single-pulse experiment using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_db entry for one single metabolite

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        # hello
        # cprint(">>> mrs.sim.mrs_sequence._run_sequence:", 'green', attrs=['bold'])

        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite)

        print("  > acquiring.", end="", flush=True)

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
        print(".90.", end="", flush=True)
        # evolution
        sigma0 = pg.evolve(sigma1, pg.prop(H, self.te / 1000.0))  # TE evolution
        te_real += self.te / 1000.0
        print(".", end="", flush=True)
        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition
        print(".", end="", flush=True)

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        print("done. (real TE=%.2fms)" % (te_real * 1000.0))

        return(s_MRSData2)

    def _compute_metabolite_signals(self):
        """For each goup of metabolites, generate the simulated MRS signal."""
        # hello
        cprint(">>> mrs.sim.mrs_sequence._compute_metabolite_signals:", 'green', attrs=['bold'])

        # clear metabase
        self._meta_signals = []

        # zero signal vector
        s_full_meta = suspect.MRSData(np.zeros([self.npts, ]), 1 / self.fs, self.f0)
        s_full_meta = s_full_meta.view(reco.MRSData2)

        # browse though the database and display everything
        for this_metagroup_key, this_metagroup_entry in self._meta_db.items():
            s = s_full_meta.copy()
            for this_meta_key, this_meta_entry in self._meta_db[this_metagroup_key]["metabolites"].items():
                print(">>> simulating MRS signal for metabolite [", end="", flush=True)
                cprint(this_metagroup_key, 'green', attrs=['bold'], end="", flush=True)
                print("/", end="", flush=True)
                cprint(this_meta_key, 'green', attrs=['bold'], end="", flush=True)
                print("]")
                # build up the metabolite group
                s = s + self._run_sequence(this_meta_entry)
            # append this metabolite group to the metabase
            self._meta_signals.append(s)

    def initialize(self, metadb=None):
        """
        Initialize the sequence object before using it to simulate MRS signals.

        Parameters
        ----------
        metadb : metabolite_db
            Metabolite database to use for simulation. If none specified, use default metabolite_db object.
        """
        # hello
        cprint(">> mrs.sim.mrs_sequence.initialize:", 'green', attrs=['bold'])
        print(" > initializing sequence...")

        # want to use a custom metabolite db?
        if(metadb is not None):
            self._meta_db = metadb

        # now let's initialize the metabolite db if needed
        if(not self._meta_db.ready):
            self._meta_db.initialize()

        self._init_pulses()
        self._compute_metabolite_signals()
        dt = 1 / self.fs
        self._t = np.arange(0, self.npts * dt, dt)
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
            raise Exception("> This mrs_sequence object was not initialized!")

        # time MRS model
        s_MRSData = suspect.MRSData(np.zeros([self.npts, ]), 1 / self.fs, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)
        for k, s_meta in enumerate(self._meta_signals):
            s_MRSData2 = s_MRSData2 + s_meta * p[k, 0] * np.exp((-p[k, 1] + 2.0 * np.pi * 1j * p[k, 2]) * self._t + 1j * (p[k, 3] + self.additional_phi0))

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
            raise Exception("> This mrs_sequence object was not initialized!")

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

    def simulate_signal(self, p, sigma_noise=0.0):
        """
        Print out the parameter values and returns the modeled signal using above member function.

        Parameters
        ----------
        p : params object
            Array of simulation parameters

        Returns
        -------
        s : MRSData2 numpy array [timepoints]
            Simulated MRS data signal stored in a MRSData2 object
        """
        # hello
        cprint(">> mrs.sim.mrs_sequence.simulate_signal:",
               'green', attrs=['bold'])

        # ready or not, here I come
        if(not self.ready):
            raise Exception("> This mrs_sequence object was not initialized!")

        # checking if LL is not broken
        if(not p.check()):
            raise Exception(" > mrs.sim.mrs_sequence.simulate_signal the link-lock vector of this params object looks broken!")

        print(" > simulating signal")
        s = self._model(p)

        if(sigma_noise > 0.0):
            print(" > adding complex gaussian noise, std. deviation = " + str(sigma_noise))
            s = s + np.random.normal(0.0, sigma_noise, s.shape[0] * 2).view(np.complex128)

        # adding te (fix)
        s._te = self.te
        return(s)


class mrs_seq_press(mrs_sequence):
    """A class that represents a general PRESS sequence."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0, te1=None, te2=None):
        """
        Initialize a virtual PRESS sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        te1 : float
            First part of TE (ms)
        te2 : float
            Second part of TE (ms)
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor)
        # name of sequence
        self._name = "press (not specific)"
        # TE timing
        if(te1 is None or te2 is None):
            # symmetric TE by default
            self.te1 = te / 2.0
            self.te2 = te / 2.0
        else:
            self.te1 = te1 / 2.0
            self.te2 = te2 / 2.0
        # flip phase for PRESS
        self.additional_phi0 = np.pi

    def _run_sequence(self, metabolite):
        """
        Simulate the NMR signal acquired using a PRESS sequence using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_db entry for one single metabolite

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        # hello
        # cprint(">>> mrs.sim.mrs_seq_laser._run_sequence:", 'green', attrs=['bold'])

        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite)

        print("  > acquiring.", end="", flush=True)

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
        if(self.te1 is None and self.te2 is None):
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
            print("." + str(p) + ".", end="", flush=True)
            # evolution
            sigma0 = pg.evolve(sigma1, pg.prop(H, d))
            te_real += d
            print(".", end="", flush=True)

        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition
        print(".", end="", flush=True)

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        print("done. (real TE=%.2fms)" % (te_real * 1000.0))

        return(s_MRSData2)


class mrs_seq_steam(mrs_sequence):
    """A class that represents a general STEAM sequence."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor)
        # name of sequence
        self._name = "steam (not specific)"
        # mixing time
        self.tm = tm
        # need some 180deg phase here
        self.additional_phi0 = np.pi

    def _run_sequence(self, metabolite):
        """
        Simulate the NMR signal acquired using a STEAM sequence using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_db entry for one single metabolite

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        # hello
        # cprint(">>> mrs.sim.mrs_seq_laser._run_sequence:", 'green', attrs=['bold'])

        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite)

        print("  > acquiring.", end="", flush=True)

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
        print(".90.", end="", flush=True)
        # evolution
        sigma0 = pg.evolve(sigma1, pg.prop(H, evol_delays_s[0]))
        te_real += evol_delays_s[0]
        print(".", end="", flush=True)

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

        print(".90...90.", end="", flush=True)
        # last TE/2 nutation
        sigma0 = pg.evolve(sigma_res, pg.prop(H, evol_delays_s[2]))
        te_real += evol_delays_s[2]
        print(".", end="", flush=True)

        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition
        print(".", end="", flush=True)

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        print("done. (real TE=%.2fms)" % (te_real * 1000.0))

        return(s_MRSData2)


class mrs_seq_eja_svs_slaser(mrs_sequence):
    """A class that represents the semi-LASER sequence from CMRR."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0, exc_pulse_duration=2.0, exc_pulse_voltage=350.0, rfc_pulse_duration=9.0, rfc_pulse_fa=180.0, rfc_pulse_r=20.0, rfc_pulse_n=1.0, rfc_pulse_voltage=350.0, ref_pulse_voltage=300.0, spoiler_duration=1.0):
        """
        Initialize a virtual semi-LASER sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
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
        ref_pulse_voltage : float
            Reference voltage (V)
        spoiler_duration : float
            Spoiler duration (ms)
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor)
        # name of sequence
        self._name = "eja_svs_slaser"
        # excitation and AFP pulse properties
        self.pulse_exc_duration = exc_pulse_duration
        self.pulse_rfc_duration = rfc_pulse_duration
        self.pulse_rfc_flipangle = rfc_pulse_fa
        self.pulse_rfc_r = rfc_pulse_r
        self.pulse_rfc_n = rfc_pulse_n
        self.pulse_rfc_voltage = rfc_pulse_voltage
        self.pulse_rfc_npts = None  # no default number of points
        self.pulse_ref_voltage = ref_pulse_voltage
        # spoiler properties
        self.spoiler_duration = spoiler_duration

        # should we use real shaped RF pulses or not
        self.pulse_rfc_real_shape_enable = True

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
        # [2] 30% power increase to be sure to be in adiabatic regime
        self.pulse_rfc_optim_power_margins = [10.0, 10.0, 30.0]  # %
        # display all this in a nice fig
        self.pulse_rfc_optim_power_display = True
        # final pulse voltage found by optimization
        self.pulse_rfc_optim_power_voltage = None
        # final pulse flip angle found by optimization
        self.pulse_rfc_optim_power_flipangle = None

        # RF power (w1 in Hz) used for AFP pulses during simulation
        self._pulse_rfc_w1max = None

    def _sech(self, x):
        return(1.0 / np.cosh(x))

    def _asech(self, x):
        return(ma.acosh(1 / x))

    def _rf_pulse_hsn(self, N: int, n: int, tbw, Tp, trunc=0.01):
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

        return(pulse_amp, pulse_phi)

    def _init_pulses(self):
        """Calculate HSn pulse profiles for later use during simulation. This could look like a detail but the semi-LASER behaves in particular way regarding effective TE and J-coupling BECAUSE of those pulses and their long durations. Depending on TE and pulse properties, MR peak lineshapes can be severly affected."""
        # hello
        cprint(">>> mrs.sim.mrs_seq_laser._init_pulses:", 'green', attrs=['bold'])

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
            print("  > generating a %.1f-ms HS%d pulse: R=%.0f (BW=%.1fkHz) using %d points..." % (self.pulse_rfc_duration, self.pulse_rfc_n, self.pulse_rfc_r, pulse_bw_hz / 1000.0, pulse_n), end="", flush=True)
            pulse_amp, pulse_phi = self._rf_pulse_hsn(pulse_n, self.pulse_rfc_n, self.pulse_rfc_r, self.pulse_rfc_duration)

            # store
            self.pulses_rfc_shape_amp_norm = pulse_amp
            self.pulses_rfc_shape_phi_deg = pulse_phi
            print("done.")

            if(self.pulse_rfc_optim_power_enable):
                self._pulse_rfc_w1max = self._get_rfc_pulse_w1max_by_optim()
            else:
                self._pulse_rfc_w1max = self._get_rfc_pulse_w1max_using_ref_pulse_voltage()
        else:
            print("  > not using real shaped pulses so nothing to here!")

    def _get_rfc_pulse_w1max_using_ref_pulse_voltage(self):
        """
        Return the RF power (w1max) of the HSn pulses using the reference voltage.

        Returns
        -------
        w1_rfc_pulse_hz  : float
            RF power or w1 max for HSn pulse (Hz)
        """
        # hello
        cprint(">>> mrs.sim.mrs_seq_laser._get_rfc_pulse_w1max_using_ref_pulse_voltage:", 'green', attrs=['bold'])

        # optimize RF power
        print("  > setting HSn pulse RF power using reference voltage...")

        # pulse area
        pulse_t = np.linspace(0.0, 1.0 / 1000.0, len(self.pulses_rfc_shape_amp_norm))
        pulse_rfc_area = np.trapz(self.pulses_rfc_shape_amp_norm, pulse_t, pulse_t[1])

        # find real RF power that was used
        w1_ref_pulse_degs = 180 / (360 * 1e-3)  # deg/s
        w1_rfc_pulse_degs = self.pulse_rfc_voltage * w1_ref_pulse_degs / self.pulse_ref_voltage  # deg/s
        w1_rfc_pulse_hz = w1_rfc_pulse_degs / (360.0 * pulse_rfc_area)

        print("  > HSn pulse voltage is set to %.2fV" % self.pulse_rfc_voltage)
        print("  > HSn pulse flip angle is set to %.2fdeg" % self.pulse_rfc_flipangle)
        print("  > Reference voltage is %.2fV" % self.pulse_ref_voltage)
        print("  > This means the HSn pulse w1max is %.2fHz" % w1_rfc_pulse_hz)

        return(w1_rfc_pulse_hz)

    def _get_rfc_pulse_w1max_by_optim(self):
        """
        Optimize the RF power (w1max) of the HSn pulses. The current power used during acquisition, if available, will be displayed.

        Returns
        -------
        w1_rfc_pulse_hz_optim : float
            RF power or w1 max for HSn pulse (Hz)
        """
        # hello
        cprint(">>> mrs.sim.mrs_seq_laser._get_rfc_pulse_w1max_by_optim:", 'green', attrs=['bold'])

        # optimize RF power
        print("  > optimizing RF power...", end="", flush=True)

        # pulse area
        pulse_t = np.linspace(0.0, 1.0 / 1000.0, len(self.pulses_rfc_shape_amp_norm))
        pulse_rfc_area = np.trapz(self.pulses_rfc_shape_amp_norm, pulse_t, pulse_t[1])

        # find real RF power that was used
        w1_ref_pulse_degs = 180 / (360 * 1e-3)  # deg/s
        w1_rfc_pulse_degs = self.pulse_rfc_voltage * w1_ref_pulse_degs / self.pulse_ref_voltage  # deg/s
        w1_rfc_pulse_hz = w1_rfc_pulse_degs / (360.0 * pulse_rfc_area)

        # prepare RF power triple axis
        w1_range_V = np.linspace(self.pulse_rfc_optim_power_voltage_range[0], self.pulse_rfc_optim_power_voltage_range[1], self.pulse_rfc_optim_power_voltage_n)
        w1_range_Hz = (w1_range_V * w1_ref_pulse_degs / self.pulse_ref_voltage) / (360.0 * pulse_rfc_area)

        # in uT?
        # gamma = GAMMA_DICT[self.nuclei]
        # w1_range_uT = w1_range_Hz / gamma

        # run RF adjustment
        meta_db_keys = list(self.meta_db.keys())
        acquired_signals = []
        peak_max_ppm_index = []
        peak_intensity_abs = np.zeros([len(self.pulse_rfc_optim_power_metabolites), len(w1_range_Hz)])
        peak_intensity_real = np.zeros([len(self.pulse_rfc_optim_power_metabolites), len(w1_range_Hz)])

        # for each RF power (w1)
        for kw, w in enumerate(w1_range_Hz):
            self._pulse_rfc_w1max = w
            # for each metabolite
            for km, im in enumerate(self.pulse_rfc_optim_power_metabolites):
                meta_key = meta_db_keys[im]
                meta_dict_entry = self.meta_db[meta_key]["metabolites"][meta_key]
                # acquire
                s = self._run_sequence(meta_dict_entry, True)
                s = s._correct_apodization(10.0)  # silent apodization
                acquired_signals.append(s)
                # analyse
                sf = s.spectrum()
                if(kw == 0):
                    # first power, let's peak pick
                    peak_max_ppm_index.append(np.argmax(np.abs(sf)))
                # store peak intensities
                peak_intensity_abs[km, kw] = np.abs(sf[peak_max_ppm_index[km]])
                peak_intensity_real[km, kw] = np.real(sf[peak_max_ppm_index[km]])

            print(".", end="", flush=True)

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
        self.pulse_rfc_optim_power_voltage = w1_range_degs_optim * self.pulse_ref_voltage / w1_ref_pulse_degs  # V

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
        fig.canvas.set_window_title("mrs.sim.mrs_seq_eja_svs_slaser._get_rfc_pulse_w1max_by_optim")

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
            meta_key = meta_db_keys[im]
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

        print("done.")

        print("  > HSn pulse voltage is set to %.2fV" % self.pulse_rfc_voltage)
        print("  > HSn pulse flip angle is set to %.2fdeg" % self.pulse_rfc_flipangle)
        print("  > Reference voltage is %.2fV" % self.pulse_ref_voltage)
        print("  > Which is equivalent to a w1max of %.2fHz" % w1_rfc_pulse_hz)
        print("  > The optimization process gave however something different...")
        print("  > The optimal w1max found is %.2fHz" % w1_range_Hz_optim)
        print("  > Which corresponds to a voltage of %.2fV" % self.pulse_rfc_optim_power_voltage)
        print("  > Or a flip angle of %.2fdeg" % self.pulse_rfc_optim_power_flipangle)

        return(w1_range_Hz_optim)

    def _run_sequence(self, metabolite, bSilent=False):
        """
        Simulate the NMR signal acquired using a sLASER sequence using the GAMMA library via the pyGAMMA python wrapper for one metabolite (=one spin system).

        Parameters
        ----------
        metabolite : dict
            metabolite_db entry for one single metabolite
        beSilent : boolean
            No output in console (True)

        Returns
        -------
        s_MRSData2 : MRSData2 object
            Containing the MRS signal simulated for this metabolite
        """
        # hello
        if(not bSilent):
            cprint(">>> mrs.sim.mrs_seq_laser._run_sequence:", 'green', attrs=['bold'])

        # create spin system
        sys, amp_factor_spins = self._prepare_spin_system(metabolite, bSilent)

        if(not bSilent):
            print("  > acquiring.", end="", flush=True)

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
            raise Exception("  > mrs.sim.mrs_seq_eja_svs_slaser._run_sequence: your echo time is too short!")
        gh = self.pulse_rfc_duration + 2.0 * self.spoiler_duration + te_fill
        # final SE timing
        g = gh / 2.0
        h = gh / 2.0
        fg = f + g

        # now important: since rfc pulses here have a real durations, we need to convert the later a, bc, ... delays to time padding delays: 3 cases possible
        if(self.pulse_rfc_real_shape_enable and self.allow_evolution_during_hard_pulses):
            # we are using real shaped rfc pulses and a fake exc pulse
            # we want to keep  time evolution during hard pulsing (and keep TE correct)
            # remove all rfc pulse durations from delays because they are taken into account during simulation
            # BUT do not remove exc pulse duration because it is a hard pulse
            a = a - self.pulse_rfc_duration / 2.0
            bc = bc - self.pulse_rfc_duration
            de = de - self.pulse_rfc_duration
            fg = fg - self.pulse_rfc_duration
            h = h - self.pulse_rfc_duration / 2.0
        elif(not self.pulse_rfc_real_shape_enable and self.allow_evolution_during_hard_pulses):
            # we are using zero-duration exc and rfc hard pulses
            # we want to keep  time evolution during hard pulsing (and keep TE correct)
            # do not remove any pulse duration
            a = a
            bc = bc
            de = de
            fg = fg
            h = h
        elif(not self.pulse_rfc_real_shape_enable and not self.allow_evolution_during_hard_pulses):
            # we are using zero-duration exc and rfc hard pulses
            # we do not want to keep time evolution during hard pulsing (TE will be virtually shortened)
            # remove all pulse duration
            a = a - self.pulse_exc_duration / 2.0 - self.pulse_rfc_duration / 2.0
            bc = bc - self.pulse_rfc_duration
            de = de - self.pulse_rfc_duration
            fg = fg - self.pulse_rfc_duration
            h = h - self.pulse_rfc_duration / 2.0

        # final inter pulse delay list in seconds
        evol_delays_s = np.array([a, bc, de, fg, h]) / 1000.0

        # run the sequence
        sigma0 = pg.sigma_eq(sys)
        te_real = 0.0

        # excitation: hard 90 pulse to simulate asymmetric sinc pulse
        sigma1 = pg.Iypuls(sys, sigma0, self.nuclei, 90.0)
        te_real += 0.0  # fake pulse
        if(not bSilent):
            print(".90.", end="", flush=True)
        # evolution
        sigma0 = pg.evolve(sigma1, pg.prop(H, evol_delays_s[0]))
        te_real += evol_delays_s[0]
        if(not bSilent):
            print(".", end="", flush=True)

        # LASER: 2 x (pair of 180 HSn pulses)
        for d in evol_delays_s[1:]:
            if(self.pulse_rfc_real_shape_enable):
                sigma1 = Upulc180.evolve(sigma0)
                te_real += self.pulse_rfc_duration / 1000.0
            else:
                sigma1 = pg.Iypuls(sys, sigma0, self.nuclei, 180.0)
                te_real += 0.0  # fake pulse

            if(not bSilent):
                print(".180.", end="", flush=True)

            # evolution
            sigma0 = pg.evolve(sigma1, pg.prop(H, d))
            te_real += d
            if(not bSilent):
                print(".", end="", flush=True)

        data = pg.FID(sigma0, pg.gen_op(D), H, dt, self.npts)  # acquisition

        if(not bSilent):
            print(".", end="", flush=True)

        # extract complex time data points
        s = np.full([self.npts, ], np.nan, dtype=np.complex)
        for i in range(self.npts):
            s[i] = -afactor * amp_factor_spins * (data.getRe(i) - 1j * data.getIm(i))

        # convert to suspect
        s_MRSData = suspect.MRSData(s, dt, self.f0)
        s_MRSData2 = s_MRSData.view(reco.MRSData2)

        if(not bSilent):
            print("done. (real TE=%.2fms)" % (te_real * 1000.0))

        return(s_MRSData2)


class mrs_seq_eja_svs_press(mrs_seq_press):
    """A class that represents the PRESS sequence from CMRR."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0, te1=None, te2=None):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        te1 : float
            First part of TE (ms)
        te2 : float
            Second part of TE (ms)
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor, te1, te2)
        # name of sequence
        self._name = "eja_svs_press"


class mrs_seq_eja_svs_steam(mrs_seq_steam):
    """A class that represents the STEAM sequence from CMRR."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor, tm)
        # name of sequence
        self._name = "eja_svs_steam"


class mrs_seq_fid(mrs_sequence):
    """A class that represents the pulse-acquire sequence (fid), which is actually a clone of the super class."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0):
        """
        Initialize a virtual FID sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor)
        # name of sequence
        self._name = "fid"


class mrs_seq_svs_st(mrs_seq_steam):
    """A class that represents the STEAM sequence from SIEMENS."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor, tm)
        # name of sequence
        self._name = "svs_st"


class mrs_seq_svs_st_vapor_643(mrs_seq_steam):
    """A class that represents the STEAM sequence from SIEMENS WiP 643."""

    def __init__(self, te, tr=3500.0, nuclei="1H", npts=4096 * 4, fs=5000.0, f0=297.2062580, scaling_factor=1.0, tm=20.0):
        """
        Initialize a virtual STEAM sequence.

        Parameters
        ----------
        te : float
            Echo time (ms)
        tr : float
            Repetition time (ms)
        nuclei : string
            Observed nuclei. Examples: "1H", "31P", etc.
        npts : int
            Number of acquisition points
        fs : float
            Acquisition bandwidth (Hz)
        f0 : float
            Water Larmor frequency (MHz)
        scaling_factor : float
            Scaling FID intensity factor
        tm : float
            Mixing time (ms)
        """
        super().__init__(te, tr, nuclei, npts, fs, f0, scaling_factor, tm)
        # name of sequence
        self._name = "svs_st_vapor_643"


class metabolite_db(dict):
    """The metabolite_db class is a big dictionnary of metabolites with their respective chemical shift and J coupling information."""

    def __init__(self):
        """Construct a metabolite object."""
        super(metabolite_db, self).__init__()
        # xls file that contains all the info
        self.xls_file = "./mrs/metabase.xls"
        # include peaks only in this range of chemical shifts
        self.ppm_range = [0, 6]
        # include only non-coupled peaks
        self.non_coupled_only = False
        # this is the name of the first MM
        self._first_MM = "MM1"
        # to know if the object is initialized
        self._ready = False

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
        # hello
        cprint(">> mrs.sim.metabolite_db._read_xls_file:", 'green', attrs=['bold'])

        # remove annoying warning
        warnings.simplefilter(action='ignore', category=FutureWarning)

        # first, read xls metabolite basis set file and parse it
        print(" > reading metabolite basis set from XLS file...")
        book = open_workbook(self.xls_file)

        # get sheets
        all_sheets = book.sheet_names()
        sheet_names = book.sheet_by_name('Names')
        sheet_groups = book.sheet_by_name('Groups')
        sheet_ppm = book.sheet_by_name('PPM')
        sheet_iso = book.sheet_by_name('Nuclei')

        # find extra sheet names
        extra_sheet_names = [s for s in all_sheets if "Extra" in s]

        # metagroup list
        for this_metagroup_name in sheet_groups.col_values(0):
            # index
            this_metagroup_ind = sheet_groups.col_values(0).index(this_metagroup_name)

            # prepare group entry
            this_metagroup_entry = {"metabolites": {}}

            # look for the members of the group
            for this_meta_name in sheet_groups.row_values(this_metagroup_ind):
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
                    if(this_meta_name in all_sheets):
                        this_sheet_Jcouplings = book.sheet_by_name(this_meta_name)
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
                    this_metagroup_entry["metabolites"].update({this_meta_name: {"ppm": this_ppm_list, "iso": this_iso_list, "J": this_j_list}})

            # add extra infos to metabolite group from extra sheets (T2, T1, etc.)
            for this_extra_sheet_name in extra_sheet_names:
                # extract string from sheet name
                this_extra_sheet = book.sheet_by_name(this_extra_sheet_name)
                this_col = np.array(this_extra_sheet.col_values(0))
                this_extra_info_val = this_col[this_metagroup_ind]
                this_extra_info_key = this_extra_sheet_name[6:]  # remove Extra_
                this_metagroup_entry.update({this_extra_info_key: this_extra_info_val})

            # and add the group to the dict
            self.update({this_metagroup_name: this_metagroup_entry})

    def _write_header_file(self):
        """Interesting method here... It generates a python .py and writes very usefull aliases to access quickly to a specific metabolite or parameter. Since metabolite indexes depend on the metabolite database, this python header file is regenerated each time the metabolite basis set is initialized."""
        # hello
        cprint(">> mrs.sim.metabolite_db._write_header_file:", 'green', attrs=['bold'])
        print(" > generating metabolite and parameter aliases...")

        with open("./mrs/metabase.py", 'w') as f:
            f.write("#!/usr/bin/env python3")
            f.write("# -*- coding: utf-8 -*-\n")
            f.write('\n')
            f.write("# This file is automatically generated on the fly by sim.py!\n")
            f.write("\n")
            metagroup_names = list(self.keys())
            for k, thisMetaName in enumerate(metagroup_names):
                f.write("m_" + thisMetaName + " = " + str(k) + "\n")
            f.write("\n")
            f.write("p_cm = 0\n")
            f.write("p_dd = 1\n")
            f.write("p_df = 2\n")
            f.write("p_dp = 3\n")
            f.write("\n")

            f.write("m_All_MBs = [0")
            for im in range(1, metagroup_names.index(self._first_MM)):
                f.write(", " + str(im))
            f.write("]\n")

            f.write("m_All_MMs = [" + str(metagroup_names.index(self._first_MM)))
            for im in range(metagroup_names.index(self._first_MM) + 1, len(metagroup_names)):
                f.write(", " + str(im))
            f.write("]\n")
            f.write("\n")

            f.write("n_ALL = " + str(len(metagroup_names)) + "\n")
            f.write("n_MBs = " + str(metagroup_names.index(self._first_MM)) + "\n")
            f.write("n_MMs = " + str(len(metagroup_names) - metagroup_names.index(self._first_MM)) + "\n")
            f.write("\n")

    def print_out(self):
        """Print the metabolite database with all the information (chemical shifts, nuclei, J couplings)."""
        # hello
        cprint(">> mrs.sim.metabolite_db.print_out:", 'green', attrs=['bold'])

        # and now browse though the database and display everything
        for this_metagroup_key, this_metagroup_entry in self.items():
            print("")

            print(" > metabolite [", end="", flush=True)
            cprint(this_metagroup_key, 'green', attrs=['bold'], end="", flush=True)
            print("/", end="", flush=True)

            for this_meta_key, this_meta_entry in self[this_metagroup_key]["metabolites"].items():
                cprint(this_meta_key, 'green', attrs=['bold'], end="", flush=True)
                print("]")

                # first nuclei line
                print("nuclei  \t", end="", flush=True)
                for this_iso in this_meta_entry["iso"]:
                    print("\t(%2d)" % this_iso, end="", flush=True)
                print()

                print("\t     ppm", end="", flush=True)
                for this_ppm in this_meta_entry["ppm"]:
                    print("\t%4.1f" % this_ppm, end="", flush=True)
                print()

                for (this_iso, this_ppm, irowJ) in zip(this_meta_entry["iso"], this_meta_entry["ppm"], range(len(this_meta_entry["ppm"]))):
                    # recall nuclei and chemical shift for 1st and 2nd col
                    print("(%2d)" % this_iso, end="", flush=True)
                    print("\t%4.1f" % this_ppm, end="", flush=True)

                    for icolJ in range(len(this_meta_entry["J"])):
                        print("\t%4.1f" % this_meta_entry["J"][irowJ][icolJ], end="", flush=True)
                    print()

    def initialize(self):
        """Initialize metabolite database: run the two previous methods."""
        # hello
        cprint(">> mrs.sim.metabolite_db.initialize:", 'green', attrs=['bold'])
        print(" > initializing metabolite database...")
        self.clear()
        self._read_xls_file()
        self._write_header_file()
        # self.print_out()
        self._ready = False


def disp_bargraph(ax, params_val_list, params_std_list, params_leg_list, colored_LLs=True, bMM=False, bWater=False, bFitMode=False, pIndex=xxx.p_cm, mIndex_list=None, width=0.3):
    """
    Plot a bargraph of concentrations.

    Parameters
    ----------
    ax : matplotlib axis
        Axis to use for plotting
    params_val_list : list of params objects
        List of parameter arrays to display as bars
    params_std_list : list of params objects
        List of parameter arrays to display as error bars
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
        meta_mask[xxx.m_MM1:, :] = False

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
        raise Exception("  > mrs.sim.disp_bargraph: nothing to display! Please check disp_bargraph paramters like pIndex, mIndex_list, bWater and bMM!")

    # filter data now
    params_val_list = [p[meta_mask] for p in params_val_list]
    params_LL_list = [p.linklock[meta_mask] for p in params_val_list]
    params_std_list = [p[meta_mask] for p in params_std_list]
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
    ax.set_xticklabels(meta_names, rotation=90, fontsize=9)
    ax.grid('on')
    if(not bFitMode):
        ax.legend()


def disp_fit(ax, data, params, seq, LL_exluding=True, LL_merging=False, mIndex_list=None, display_range=[1, 5]):
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
    ygap = np.max(mod.spectrum().real) * (60 - n_meta_to_display) * 0.005

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
    ymax = np.max(mod.spectrum().real) * 1.1  # max(np.max(mod.spectrum().real), np.max(data.spectrum().real))
    ax.set_ylim([-ygap * (koffset + 3), ymax])
    ax.set_yticks([])
    ax.set_xlabel('chemical shift (ppm)')
    ax.grid('on')
    ax.legend(loc="upper right")
    plt.tight_layout()
