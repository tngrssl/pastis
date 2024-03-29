#! / usr / bin / env python3
# -*- coding: utf-8 -*-
"""
Three classes' definition in here.

    * a MRSData2 class with a bunch of methods based on the suspect module to deal with SIEMENS 7T MRS data
    * a pipeline class to run the reconstruction process on a bunch of acquired data

@author: Tangi Roussel
"""

import suspect
import numpy as np
import pandas as pd
from scipy import signal
from scipy import interpolate
import matplotlib.pylab as plt
import matplotlib._pylab_helpers
from datetime import datetime
import pickle
import os
import re
from enum import Enum
from pastis import io
from pastis import log

import pdb

# max number of datasets in a pipeline
MAX_NUM_DATASETS = 1000

# constants used during data rejection
# minimum step while setting amplitude threshold (% of signal relative change)
DATA_REJECTION_AMPLITUDE_STEP = 2.0
# minimum step while setting linewidth threshold (Hz)
DATA_REJECTION_LINEWIDTH_STEP = 1.0
# minimum step while setting chemical shift threshold (ppm)
DATA_REJECTION_FREQUENCY_STEP = 0.001
# minimum step while setting phase threshold (rad)
DATA_REJECTION_PHASE_STEP = 0.1

# spectral resolution for peak realignment using inter-correlation mode (experimental) (ppm)
RECO_CORRECT_REALIGN_INTER_CORR_MODE_DF = 0.1


class suspect_phasing_method(Enum):
    """The enum suspect_phasing_method describes the type of phasing method to use (phasing method from the suspect package."""

    MATCH_MAGNITUDE_REAL = 1
    MIN_IMAG_INTEGRAL = 2
    ACME = 3


class data_rejection_method(Enum):
    """The enum data_rejection_method describes the type of data rejection method used."""

    AUTO_AMPLITUDE = 0
    AUTO_LINEWIDTH = 1
    AUTO_FREQUENCY = 2
    AUTO_PHASE = 3


class MRSData2(suspect.mrsobjects.MRSData):
    """A class based on suspect's MRSData to store MRS data."""

    def __new__(cls, data_filepath, physio_log_file=None, anatomy_folderpath=None, obj=None, dt=None, f0=None, te=None, tr=None, ppm0=None, voxel_dimensions=None, transform=None, metadata=None, data_ref=None, label="", offset_display=0.0, patient={}, sequence_obj=None, noise_level=None, data_rejection=None, data_file_hash=None, is_concatenated=None, is_rawdata=None):
        """
        Construct a MRSData2 object that inherits of Suspect's MRSData class. In short, the MRSData2 class is a copy of MRSData + my custom methods for post-processing. To create a MRSData2 object, you need give a path that points to a SIEMENS DICOM or a SIEMENS TWIX file.

        Parameters
        ----------
        data_filepath : string
            Full absolute file path pointing to the stored signal (DCM or TWIX file) or the folder assuming that a dcm file named "original-primary_e09_0001.dcm" is stored inside.
        physio_log_file : string
            Full absolute file path pointing to a IDEA VB17 respiratory log file
        anatomy_folderpath : string
            Full absolute file path pointing a folder containing dicom DCM files contaning the anatomical images on which you want to dispaly the VOI
        obj,dt,f0,te,tr,ppm0,voxel_dimensions,transform,metadata
            Please check suspect's MRSData class for those arguments
        data_ref : MRSData2 object
            Reference data acquired for this signal
        label : string
            Label for this signal
        offset_display : float
            Y-axis offset display
        patient : dict
            Patient meta data
        sequence_obj : sim.mrs_sequence object
            Sequence object
        noise_level : float
            Noise level measured on real FID
        data_rejection : list of dicts
            Data rejection results (NA, SNR, LW, etc.)
        data_file_hash : string
            MD5 hash code of data file content
        is_concatenated : boolean
            Was this signal the result of a concatenation?
        is_rawdata : boolean
            Was this signal read from a DICOM file?

        Returns
        -------
        obj : MRSData2 numpy array [averages,channels,timepoints]
            Resulting constructed MRSData2 object
        """
        if(data_filepath == []):
            # calling the parent class' constructor
            obj = super(suspect.mrsobjects.MRSData, cls).__new__(cls, obj, dt, f0, te, tr, ppm0, voxel_dimensions, transform, metadata)
            # adding attributes
            obj.data_ref = data_ref
            obj._display_label = label
            obj._display_offset = offset_display
            obj._tr = tr
            obj._patient = patient
            obj._sequence = sequence_obj
            obj._noise_level = noise_level
            obj._data_rejection = data_rejection
            obj._data_file_hash = data_file_hash
            obj._is_concatenated = is_concatenated
            obj._is_rawdata = is_rawdata
            # bye
            return(obj)

        # hello
        log.debug("creating object...")

        # open data file
        log.info("reading data file...")
        log.info(data_filepath)
        # read data and header
        mfr = io.get_data_file_reader(data_filepath)
        # read data and get a suspect MRSData object
        MRSData_obj = mfr.data
        # get hash code
        hc = mfr.get_md5_hash()

        # --- reshape data ---

        # add dimensions if needed
        if(MRSData_obj.ndim == 1):
            # this could be a single-shot single-rx signal in twix
            # or an averaged single-rx signal in dicom...

            # add 2 dimensions
            MRSData_obj = MRSData_obj.reshape((1, 1,) + MRSData_obj.shape)

        if(MRSData_obj.ndim == 2):
            # this could be a single-shot multi-rx signal
            # or a averaged single-rx signal...
            coil_nChannels = mfr.get_number_rx_channels()

            if(coil_nChannels == MRSData_obj.shape[0]):
                # beware, same number of averages / coil elements, which is which?
                log.warning("ambiguous data dimensions: " + str(MRSData_obj.shape) + "-> Assuming (" + str(MRSData_obj.shape[0]) + ") to be the coil channels!")
                MRSData_obj = MRSData_obj.reshape((1,) + MRSData_obj.shape)
            elif(coil_nChannels == 1):
                # ok the user said it is a single channel coil, so the 2nd dimension here is the number of averages!
                # adding 1 dimension for number of channels
                MRSData_obj = MRSData_obj.reshape((1,) + MRSData_obj.shape)
                MRSData_obj = np.transpose(MRSData_obj, (1, 0, 2))
            else:
                # ok that is a single-shot multi-channel signal
                # adding 1 dimension for averages
                MRSData_obj = MRSData_obj.reshape((1,) + MRSData_obj.shape)

        log.info("read a " + str(MRSData_obj.shape) + " vector")

        # --- build MRSData object ---
        # calling the parent class' constructor
        obj = super(suspect.mrsobjects.MRSData, cls).__new__(cls, MRSData_obj, MRSData_obj.dt, MRSData_obj.f0, MRSData_obj.te, MRSData_obj.tr, MRSData_obj.ppm0, MRSData_obj.voxel_dimensions, MRSData_obj.transform, MRSData_obj.metadata)

        # --- get extra parameters ---

        # patient birthyear
        patient_birthday_datetime = mfr.get_patient_birthday()
        log.debug("extracted patient birthyear (" + str(patient_birthday_datetime) + ")")

        # patient sex
        patient_sex_str = mfr.get_patient_sex()
        log.debug("extracted patient sex (%s)" % patient_sex_str)

        # patient name
        patient_name_str = mfr.get_patient_name()
        log.debug("extracted patient name (%s)" % patient_name_str)

        # patient weight
        patient_weight_kgs = mfr.get_patient_weight()
        log.debug("extracted patient weight (%.2fkg)" % patient_weight_kgs)

        # patient height
        patient_height_m = mfr.get_patient_height()
        log.debug("extracted patient height (%.2fm)" % patient_height_m)

        # extract all the info to build the sequence object
        sequence_obj = mfr.get_sequence()

        # --- build MRSData2 object ---
        # adding MRSData2 attributes
        obj.data_ref = data_ref
        obj._patient = {"name": patient_name_str,
                        "birthday": patient_birthday_datetime,
                        "sex": patient_sex_str,
                        "weight": patient_weight_kgs,
                        "height": patient_height_m}

        obj._sequence = sequence_obj
        obj._noise_level = 0.0
        obj._data_rejection = None
        obj._data_file_hash = hc
        obj._is_concatenated = False
        obj._is_rawdata = mfr.is_rawdata()

        # those need to be called now, because they the attributes above
        obj.set_display_label()
        obj.set_display_offset()

        # respiratory trace if any
        if(physio_log_file is None):
            obj._physio_file = None
        else:
            # save this
            obj._physio_file = physio_log_file

        # anatomical images if any
        if(anatomy_folderpath is None):
            obj._anatomy_folderpath = None
        else:
            # save this
            obj._anatomy_folderpath = anatomy_folderpath

        return(obj)

    def __array_finalize__(self, obj):
        """
        Overload of special numpy array method called when playing around with stuff relative to object copy etc...

        Parameters
        ----------
        obj : MRSData2 numpy array [1,channels,timepoints]
        """
        super().__array_finalize__(obj)

        self.data_ref = getattr(obj, 'data_ref', None)
        # replace by a copy
        if(self.data_ref is not None):
            self.data_ref = obj.data_ref.copy()
        else:
            self.data_ref = None

        self._display_label = getattr(obj, 'display_label', None)
        self._display_offset = getattr(obj, 'display_offset', 0.0)
        self._patient = getattr(obj, 'patient', None)
        self._physio_file = getattr(obj, 'physio_file', None)
        self._anatomy_folderpath = getattr(obj, 'anatomy_folderpath', None)

        self._sequence = getattr(obj, 'sequence', None)
        # replace by a copy
        if(self.sequence is not None):
            self._sequence = obj.sequence.copy()

        self._noise_level = getattr(obj, 'noise_level', None)
        self._data_rejection = getattr(obj, 'data_rejection', None)
        self._data_file_hash = getattr(obj, 'data_file_hash', None)
        self._is_concatenated = getattr(obj, 'is_concatenated', None)
        self._is_rawdata = getattr(obj, 'is_rawdata', None)

    def inherit(self, obj):
        """
        Overload of special suspect's MRSData method to import a numpy array in here.

        Parameters
        ----------
        obj : MRSData2 numpy array [1,channels,timepoints]
            Multi-channel reference signal
        """
        obj2 = super().inherit(obj)

        obj2.data_ref = getattr(self, 'data_ref', None)
        # replace by a copy
        if(obj2.data_ref is not None):
            obj2.data_ref = self.data_ref.copy()
        else:
            self.data_ref = None

        obj2._display_label = getattr(self, 'display_label', None)
        obj2._display_offset = getattr(self, 'display_offset', 0.0)
        obj2._patient = getattr(self, 'patient', None)
        obj2._physio_file = getattr(self, 'physio_file', None)
        obj2._anatomy_folderpath = getattr(self, 'anatomy_folderpath', None)

        obj2._sequence = getattr(self, 'sequence', None)
        # replace by a copy
        if(obj2.sequence is not None):
            obj2._sequence = self.sequence.copy()

        obj2._noise_level = getattr(self, 'noise_level', 0.0)
        obj2._data_rejection = getattr(self, 'data_rejection', None)
        obj2._data_file_hash = getattr(self, 'data_file_hash', None)
        obj2._is_concatenated = getattr(self, 'is_concatenated', 0.0)
        obj2._is_rawdata = getattr(self, 'is_rawdata', 0.0)
        return(obj2)

    @property
    def display_label(self):
        """
        Property get function for display_label.

        Returns
        -------
        self._display_label : string
            Display label used in display_spectrum method
        """
        return(self._display_label)

    def set_display_label(self, lbl=""):
        """
        Set the display label.

        Parameters
        ----------
        lbl: string
            New display label for this signal
        """
        if((lbl == "" or lbl == []) and self.patient is not None):
            # create a usefull label based on patient name, sequence and object id
            new_lbl = ""
            if(self.patient["name"] is not None):
                new_lbl = new_lbl + self.patient["name"] + " | "
            if(self.sequence is not None):
                new_lbl = new_lbl + self.sequence.name + " | "
            new_lbl = new_lbl + str(id(self))
        else:
            self._display_label = lbl

    @property
    def display_offset(self):
        """
        Property get function for display_offset.

        Returns
        -------
        self._display_offset : float
            Display offset used in display_spectrum method
        """
        return(self._display_offset)

    def set_display_offset(self, ofs=0.0):
        """
        Set the display offset.

        Parameters
        ----------
        ofs: float
            New display offset for this signal
        """
        self._display_offset = ofs

    @property
    def patient(self):
        """
        Property get function for patient.

        Returns
        -------
        self._patient : dict
            Patient meta data
        """
        return(self._patient)

    @property
    def physio_file(self):
        """
        Property get function for physio_file.

        Returns
        -------
        self._physio_file : string
            Path to physio recording file
        """
        return(self._physio_file)

    @property
    def anatomy_folderpath(self):
        """
        Property get function for anatomy_folderpath.

        Returns
        -------
        self._anatomy_folderpath : string
            Full absolute file path pointing a folder containing dicom DCM files contaning the anatomical images on which you want to dispaly the VOI
        """
        return(self._anatomy_folderpath)

    @property
    def sequence(self):
        """
        Property get function for sequence.

        Returns
        -------
        self._sequence : sim.mrs_sequence object
            Sequence
        """
        return(self._sequence)

    @property
    def noise_level(self):
        """
        Property get function for noise_level.

        Returns
        -------
        self._noise_level : float
            Noise level in time-domain
        """
        return(self._noise_level)

    @property
    def data_rejection(self):
        """
        Property get function for data_rejection.

        Returns
        -------
        self._data_rejection : list of dict
            Data rejection results (NA, SNR, LW, etc.)
        """
        return(self._data_rejection)

    @property
    def data_file_hash(self):
        """
        Property get function for data_file_hash.

        Returns
        -------
        self._data_file_hash : string
            MD5 hash code of data file content
        """
        return(self._data_file_hash)

    @property
    def is_concatenated(self):
        """
        Property get function for is_concatenated.

        Returns
        -------
        self._is_concatenated : boolean
            True if this current signal is a result of a concatenation
        """
        return(self._is_concatenated)

    @property
    def is_rawdata(self):
        """
        Property get function for is_rawdata.

        Returns
        -------
        self._is_rawdata : boolean
            True if this current signal was originally read from a raw data file
        """
        return(self._is_rawdata)

    def _analyze_peak_1d(self, ppm_range, allowed_apodization=1.0):
        """
        Find peak in specific ppm range using magnitude mode and return stuff.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        ppm_range : list [2]
            Range in ppm used for peak searching
        allowed_apodization : float/boolean
            If >0 or !=False, apodize signal during peak analysis

        Returns
        -------
        peak_ppm : float
            Position in PPM of the peak
        peak_val : np.complex128
            Peak value
        peak_lw : float
            Peak linewidth in Hz
        peak_seg_ppm : numpy array
            Peak segment ppm axis
        peak_seg_val : numpy complex array
            Peak segment value
        """
        # silently zero-fill and apodize if needed
        log.pause()
        s = self.correct_zerofill_nd().correct_apodization_nd(allowed_apodization)
        log.resume()

        # init
        ppm = s.frequency_axis_ppm()
        sf = s.spectrum()
        sf_abs = np.abs(sf)

        # mask outside range
        ippm_peak_outside_range = (ppm_range[0] > ppm) | (ppm > ppm_range[1])
        sf_abs[ippm_peak_outside_range] = 0

        # max
        peak_index = np.argmax(sf_abs)

        # check
        if(peak_index == 0):
            log.error("no peak found in specified ppm range or badly phased data!")

        # ppm
        peak_ppm = ppm[peak_index]

        # complex value
        peak_val = sf[peak_index]

        # lw
        sf_real = np.real(sf)
        amp_peak = np.real(peak_val)
        ippm_max = peak_index + np.argmax(sf_real[peak_index:] < amp_peak / 2)
        ippm_min = peak_index - np.argmax(sf_real[peak_index::-1] < amp_peak / 2)
        dppm = np.abs(ppm[ippm_max] - ppm[ippm_min])
        peak_lw = dppm * s.f0

        # peak segment
        ippm_half_peak = np.arange(ippm_min, ippm_max)
        ppm_seg = ppm[ippm_half_peak]
        peak_seg = sf[ippm_half_peak]

        return(peak_ppm, peak_val, peak_lw, ppm_seg, peak_seg)

    def _analyze_peak_2d(self, peak_range=[4.5, 5], allowed_apodization=1.0):
        """
        Analyze a peak in the spectrum by estimating its amplitude, linewidth, frequency shift and phase for each average.

        * Works only with a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : array [2]
            Range in ppm used to analyze peak phase when no reference signal is specified
        allowed_apodization : float/boolean
            If >0 or !=False, apodize signal during peak analysis

        Returns
        -------
        peak_trace : numpy array [averages,4]
            Peak changes (amplitude, linewidth, frequency and phase) for each average in raw data
        peak_trace_rel2mean : numpy array [averages,4]
            Peak changes (amplitude, linewidth, frequency and phase) for each average in raw data relative to mean value
        peak_trace_rel2firstpt : numpy array [averages,4]
            Peak relative changes (amplitude, linewidth, frequency and phase) for each average in raw data relative to 1st point
        """
        log.debug("analyzing peak for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # apodize if needed
        s = self.copy()

        # first, find peak of interest in range, just to check
        s_avg = np.mean(s, axis=0)
        peak_ppm, _, _, _, _ = s_avg._analyze_peak_1d(peak_range, allowed_apodization)
        log.debug("found peak of interest at %0.2fppm!" % peak_ppm)

        # for each average in moving averaged data
        peak_trace = np.zeros([s.shape[0], 4])
        pbar = log.progressbar("analyzing", s.shape[0])
        for a in range(0, s.shape[0]):

            # call 1D peak analysis
            peak_ppm, peak_val, peak_lw, _, _ = s[a, :]._analyze_peak_1d(peak_range, allowed_apodization)

            # shift in ppm
            peak_trace[a, 2] = peak_ppm
            # amplitude
            peak_trace[a, 0] = np.real(peak_val)
            # linewidth in Hz
            peak_trace[a, 1] = peak_lw
            # phase in rad
            peak_trace[a, 3] = np.angle(peak_val)

            pbar.update(a)

        # normalize stuff

        # relative to mean
        peak_trace_rel2mean = np.zeros([s.shape[0], 4])
        peak_trace_rel2mean[:, 0] = peak_trace[:, 0] / peak_trace[:, 0].mean() * 100 - 100
        peak_trace_rel2mean[:, 1] = peak_trace[:, 1] - peak_trace[:, 1].mean()
        peak_trace_rel2mean[:, 2] = peak_trace[:, 2] - peak_trace[:, 2].mean()
        peak_trace_rel2mean[:, 3] = peak_trace[:, 3] - peak_trace[:, 3].mean()

        # relative to 1st pt
        peak_trace_rel2firstpt = np.zeros([s.shape[0], 4])
        peak_trace_rel2firstpt[:, 0] = peak_trace[:, 0] / peak_trace[0, 0] * 100 - 100
        peak_trace_rel2firstpt[:, 1] = peak_trace[:, 1] - peak_trace[0, 1]
        peak_trace_rel2firstpt[:, 2] = peak_trace[:, 2] - peak_trace[0, 2]
        peak_trace_rel2firstpt[:, 3] = peak_trace[:, 3] - peak_trace[0, 3]

        pbar.finish("done")
        return(peak_trace, peak_trace_rel2mean, peak_trace_rel2firstpt)

    def analyze_noise_nd(self, n_pts=100):
        """
        Measure noise level in time domain and store it in the "noise_level" attribute. This is usefull to keep track of the original noise level for later use, CRB normalization durnig quantification for example.

        * Works with multi-dimensional signals.
        * Returns a multi-dimensional signal.

        Parameters
        ----------
        n_pts : int
            Number of points at the end of the FID signal which we should use to estimate the noise STD

        Returns
        -------
        noise_lev : float
            Time-domain noise level
        """
        log.debug("estimating noise level for [%s]..." % self.display_label)

        s = self.copy()
        s_real = np.real(s)

        # average if needed
        while(s_real.ndim > 1):
            s_real = np.mean(s_real, axis=0)

        # init
        log.debug("estimating noise level in FID using last %d points..." % n_pts)
        # noise is the std of the last real points, but that is not so simple
        # we really want real noise, not zeros from zero_filling
        s_nonzero_mask = (s_real != 0.0)
        s_analyze = s_real[s_nonzero_mask]
        # now take the last n_pts points
        noise_lev = np.std(s_analyze[-n_pts:-1])
        log.info("noise level = %.2E" % noise_lev)

        # changing noise level attribute
        log.debug("updating noise level...")
        s._noise_level = noise_lev

        # if any ref data available, analyze noise there too
        if(s.data_ref is not None):
            s.data_ref = s.data_ref.analyze_noise_nd(n_pts)

        return(s)

    def analyze_physio_2d(self, peak_range=[4.5, 5], delta_time_range=1000.0, allowed_apodization=1.0, display=False):
        """
        Analyze the physiological signal and try to correlate it to a peak amplitude, linewidth, frequency shift and phase variations.

        * Works only with a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : array [2]
            Range in ppm used to analyze peak phase when no reference signal is specified
        delta_time_range : float
            Range in ms used to correlate / match the NMR and the physiological signal. Yes, since we are not really sure of the start timestamp we found in the TWIX header, we try to match perfectly the two signals.
        allowed_apodization : float/boolean
            If >0 or !=False, apodize signal during correction process. However, the final corrected signal will not be apodized.
        display : boolean
            Display correction process (True) or not (False)
        """
        log.debug("analyzing physiological signals for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        if(self._physio_file is None):
            # no physio signal here, exiting
            log.error("no error physiological recording file provided!")
            return()

        # read data
        [rt, rup, rd, rr] = io.read_physio_file(self._physio_file)
        resp_trace = [rt, rup, rd, rr]
        # physio signal
        resp_t = resp_trace[0]
        resp_s = resp_trace[3]

        # perform peak analysis
        peak_prop_abs, _, _ = self._analyze_peak_2d(peak_range, allowed_apodization)

        # init
        mri_t = np.linspace(self.sequence.timestamp, self.sequence.timestamp + self.sequence.tr * peak_prop_abs.shape[0], peak_prop_abs.shape[0])
        dt_array = np.arange(-delta_time_range / 2.0, delta_time_range / 2.0, 1.0)
        cc_2d = np.zeros([dt_array.shape[0], 4])

        # shift signal and calculate corr coeff
        pbar = log.progressbar("correlating signals", dt_array.shape[0])
        for idt, dt in enumerate(dt_array):
            # build time scale
            this_resp_t_interp = mri_t.copy() + dt
            this_resp_s_interp = np.interp(this_resp_t_interp, resp_t, resp_s)

            # now crop the signals to have the same length
            final_length = min(this_resp_s_interp.shape[0], peak_prop_abs.shape[0])
            this_mri_t = mri_t[0:final_length]
            this_params_trace = peak_prop_abs[0:final_length, :]
            this_resp_s_interp = this_resp_s_interp[0:final_length]

            # now remove points where resp trace is at 0 or 1
            mm = np.logical_and(this_resp_s_interp > 0, this_resp_s_interp < 1)
            this_resp_s_interp = this_resp_s_interp[mm]
            this_params_trace = this_params_trace[mm, :]

            # now, for each parameter
            for p in range(4):
                # estimate some R coeff
                cc = np.corrcoef(this_resp_s_interp, this_params_trace[:, p])
                cc_2d[idt, p] = cc[0, 1]

            pbar.update(idt)

        # find time shift that gives best correlation for each parameter
        best_dt_per_par = np.zeros(4)
        for p in range(4):
            i_maxcorr = np.argmax(np.abs(cc_2d[:, p]))
            best_dt = dt_array[i_maxcorr]
            best_dt_per_par[p] = best_dt

        # find time shift that gives best correlation for all 4 parameters
        cc_2d_all = np.sum(np.abs(cc_2d), axis=1)
        i_maxcorr = np.argmax(cc_2d_all)
        best_dt_all = dt_array[i_maxcorr]
        pbar.finish("done")

        # some info in the term
        st_ms = self.sequence.timestamp
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("data timestamp=\t" + str(st_ms) + "ms\t" + st_str)
        log.info("best start time for...")

        st_ms = self.sequence.timestamp + best_dt_per_par[0]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("amplitude=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.sequence.timestamp + best_dt_per_par[1]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("linewidth=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.sequence.timestamp + best_dt_per_par[2]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("frequency=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.sequence.timestamp + best_dt_per_par[3]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("phase=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.sequence.timestamp + best_dt_all
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("total=\t\t" + str(st_ms) + "ms\t" + st_str)

        imaxR = np.argmax(best_dt_per_par)
        best_dt = best_dt_per_par[imaxR]
        st_ms = self.sequence.timestamp + best_dt
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("max R for=\t\t" + str(st_ms) + "ms\t" + st_str)

        # time shift the signals with optimal shift
        # build time scale
        this_resp_t_interp = mri_t.copy() + best_dt
        this_resp_s_interp = np.interp(this_resp_t_interp, resp_t, resp_s)

        # now crop the signals to have the same length
        final_length = min(this_resp_s_interp.shape[0], peak_prop_abs.shape[0])
        this_mri_t = mri_t[0:final_length]
        this_params_trace = peak_prop_abs[0:final_length, :]
        this_resp_s_interp = this_resp_s_interp[0:final_length]

        # evaluate correlation coeff.
        # for each parameter
        cc_final = np.zeros(4)
        for p in range(4):
            # estimate some R coeff
            cc = np.corrcoef(this_resp_s_interp, this_params_trace[:, p])
            cc_final[p] = cc[0, 1]

        # now let's talk about FFT
        log.debug("FFT analysis...")
        nFFT = 2048
        # freq axis for resp trace
        resp_f_axis = np.fft.fftshift(np.fft.fftfreq(nFFT, d=(resp_t[1] - resp_t[0]) / 1000.0))  # Hz ou  / s
        resp_f_axis_bpm = resp_f_axis * 60.0  # / min
        # freq axis for params traces
        this_params_trace_f_axis = np.fft.fftshift(np.fft.fftfreq(nFFT, d=self.tr / 1000.0))  # Hz ou  / s
        this_params_trace_f_axis_bpm = this_params_trace_f_axis * 60.0  # / min

        # FFT of resp. trace
        resp_fft = np.abs(np.fft.fftshift(np.fft.fft((resp_s - np.mean(resp_s)) * signal.windows.hann(resp_s.shape[0]), nFFT, norm='ortho')))

        # FFT of params traces
        this_params_trace_fft = np.zeros([nFFT, 4])
        for p in range(4):
            this_params_trace_fft[:, p] = np.abs(np.fft.fftshift(np.fft.fft((this_params_trace[:, p] - np.mean(this_params_trace[:, p])) * signal.windows.hann(this_params_trace.shape[0]), nFFT, axis=0, norm='ortho'), axes=0))

        if(display):
            # display the cross-correlation plots
            fig_title = "Analyzing physiological signal [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='all')

            p = 0
            for ix in range(2):
                for iy in range(2):
                    axs[ix, iy].plot(dt_array, cc_2d[:, p], '-', linewidth=1)
                    axs[ix, iy].axvline(x=best_dt_per_par[p], color='r', linestyle='-')
                    axs[ix, iy].axvline(x=best_dt, color='r', linestyle='--')
                    axs[ix, iy].set_xlabel('time shift (ms)')
                    axs[ix, iy].grid('on')
                    p = p + 1

            axs[0, 0].set_ylabel('R amplitude vs. resp.')
            axs[0, 1].set_ylabel('R linewidth. vs. resp.')
            axs[1, 0].set_ylabel('R frequency vs. resp.')
            axs[1, 1].set_ylabel('R phase vs. resp.')

            fig.subplots_adjust()
            fig.show()

            # display time signals
            fig_title = "Physio. signal analysis (I) [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='all')

            p = 0
            for ix in range(2):
                for iy in range(2):
                    if(cc_final[p] > 0):
                        axs[ix, iy].plot(resp_t - best_dt, resp_s, 'k-', label='resp. original')
                        axs[ix, iy].plot(this_mri_t, this_resp_s_interp, 'b-', label='resp. resampled')
                    else:
                        axs[ix, iy].plot(resp_t - best_dt, 1.0 - resp_s, 'k-', label='resp. original')
                        axs[ix, iy].plot(this_mri_t, 1.0 - this_resp_s_interp, 'b-', label='resp. resampled')

                    ax2 = axs[ix, iy].twinx()
                    ax2.plot(this_mri_t, this_params_trace[:, p], 'rx-', label='MR peak property')
                    axs[ix, iy].set_xlabel('time (ms)')
                    axs[ix, iy].grid('on')
                    axs[ix, iy].legend(bbox_to_anchor=(1.05, 1), loc=2, mode='expand', borderaxespad=0)
                    ax2.legend(bbox_to_anchor=(1.05, 1), loc=3, mode='expand', borderaxespad=0)
                    p = p + 1

            axs[0, 0].set_ylabel('Rel. amplitude change (%)')
            axs[0, 1].set_ylabel('Abs. linewidth (Hz)')
            axs[1, 0].set_ylabel('Abs. frequency (Hz)')
            axs[1, 1].set_ylabel('Abs. phase shift (rd)')

            fig.subplots_adjust()
            fig.show()

            # display correlation plots
            fig_title = "Physio. signal analysis (II) [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='all')

            p = 0
            for ix in range(2):
                for iy in range(2):
                    axs[ix, iy].scatter(this_resp_s_interp, this_params_trace[:, p])
                    axs[ix, iy].set_xlabel('resp. (u.a)')
                    axs[ix, iy].grid('on')
                    axs[ix, iy].set_title("R=" + str(cc_final[p]))
                    p = p + 1

            axs[0, 0].set_ylabel('Rel. amplitude change (%)')
            axs[0, 1].set_ylabel('Abs. linewidth (Hz)')
            axs[1, 0].set_ylabel('Abs. frequency (Hz)')
            axs[1, 1].set_ylabel('Abs. phase shift (rd)')

            fig.subplots_adjust()
            fig.show()

            # display FFT plots
            fig_title = "Physio. signal analysis (III) [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='all')

            p = 0
            for ix in range(2):
                for iy in range(2):
                    axs[ix, iy].plot(resp_f_axis_bpm, resp_fft, 'k-', label='resp.')
                    ax2 = axs[ix, iy].twinx()
                    ax2.plot(this_params_trace_f_axis_bpm, this_params_trace_fft[:, p], 'r-', label='MR peak property')
                    axs[ix, iy].set_xlabel('frequency (BPM or 1 / min)')
                    axs[ix, iy].grid('on')
                    axs[ix, iy].set_xlim(0, 60.0)
                    axs[ix, iy].legend(bbox_to_anchor=(1.05, 1), loc=2, mode='expand', borderaxespad=0)
                    ax2.legend(bbox_to_anchor=(1.05, 1), loc=3, mode='expand', borderaxespad=0)
                    p = p + 1

            axs[0, 0].set_ylabel('Rel. amplitude change (FFT)')
            axs[0, 1].set_ylabel('Abs. linewidth (FFT)')
            axs[1, 0].set_ylabel('Abs. frequency (FFT)')
            axs[1, 1].set_ylabel('Abs. phase shift (FFT)')

            fig.subplots_adjust()
            fig.show()

        # done

    def analyze_snr_1d(self, peak_range, noise_range=[-2, -1], half_factor=False, magnitude_mode=False, display=False, display_range=[1, 6]):
        """
        Estimate the SNR of a peak in the spectrum ; chemical shift ranges for the peak and the noise regions are specified by the user. Can also look at time-domain SNR. Works only for a 1D MRSData2 objects.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to find a peak of interest
        noise_range : list [2]
            Range in ppm used to estimate noise
        half_factor : float
            If (True), will divide the SNR by 2 folowing an old definition of the SNR where the std of noise is multiplied by two. Btw, this outdated definition is used by LCModel.
        magnitude_mode : boolean
            analyze signal in magnitude mode (True) or the real part (False)
        display : boolean
            Display process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        snr : float
            Resulting SNR value
        s : float
            Resulting signal value
        n : float
            Resulting noise value
        """
        log.debug("analyzing SNR for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        # display
        if(display):
            fig_title = "SNR estimation [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 1, sharex='all', sharey='all')

        # find maximum peak in range
        sf = s.spectrum()
        # analyze peak WITHOUT ANY apodization (important)
        ppm_peak, peak_val, _, _, _ = s._analyze_peak_1d(peak_range, 0.0)
        if(magnitude_mode):
            log.debug("measuring the MAGNITUDE intensity at %0.2fppm!" % ppm_peak)
            snr_signal = np.abs(peak_val)
            sf_noise = np.abs(sf)
        else:
            log.debug("measuring the REAL intensity at %0.2fppm!" % ppm_peak)
            snr_signal = np.real(peak_val)
            sf_noise = np.real(sf)

        # estimate noise in user specified spectral region
        ppm = s.frequency_axis_ppm()
        log.debug("estimating noise from %0.2f to %0.2fppm region!" % (noise_range[0], noise_range[1]))
        ippm_noise_range = (noise_range[0] < ppm) & (ppm < noise_range[1])
        snr_noise = np.std(sf_noise[ippm_noise_range])
        if(half_factor):
            snr_noise = 2 * snr_noise

        if(display):
            axs[0].plot(ppm, np.real(sf), 'k-', linewidth=1)
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('real part')
            axs[0].grid('on')

            axs[1].plot(ppm, np.abs(sf), 'k-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('magnitude mode')
            axs[1].grid('on')

            if(magnitude_mode):
                ax = axs[1]
            else:
                ax = axs[0]

            # show peak of interest
            ax.plot(ppm_peak, snr_signal, 'ro')
            ax.axvline(x=ppm_peak, color='r', linestyle='--')
            # show noise region
            ax.plot(ppm[ippm_noise_range], sf_noise[ippm_noise_range], 'bo')

        # finish display
        if(display):
            fig.subplots_adjust()
            fig.show()

        # that's it
        snr = snr_signal / snr_noise
        log.info("results for [" + s.display_label + "] coming...")
        log.info("S = %.2E, N = %.2E, SNR = %0.2f!" % (snr_signal, snr_noise, snr))

        return(snr, snr_signal, snr_noise)

    def analyze_linewidth_1d(self, peak_range, magnitude_mode=False, display=False, display_range=[1, 6]):
        """
        Estimate the linewidth of a peak in the spectrum ; chemical shift ranges for the peak and the noise regions are specified by the user.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to find a peak of interest
        magnitude_mode : boolean
            analyze signal in magnitude mode (True) or the real part (False)
        display : boolean
            Display process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        lw : float
            Linewidth in Hz
        """
        log.debug("analyzing peak linewidth for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()

        # call 1D peak analysis WITHOUT ANY apodization (important)
        ppm_peak, _, lw, peak_seg_ppm, peak_seg_val = s._analyze_peak_1d(peak_range, False)
        if(magnitude_mode):
            log.debug("estimating the MAGNITUDE peak linewidth at %0.2fppm!" % ppm_peak)
        else:
            log.debug("estimating the REAL peak linewidth at %0.2fppm!" % ppm_peak)

        log.info("results for [" + s.display_label + "] coming...")
        log.info("LW = %0.2f Hz!" % lw)

        if(display):
            ppm = s.frequency_axis_ppm()
            sf = s.spectrum()

            fig_title = "FWHM estimation [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 1, sharex='all', sharey='all')

            axs[0].plot(ppm, np.real(sf), 'k-', linewidth=1)
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('real part')
            axs[0].grid('on')

            axs[1].plot(ppm, np.abs(sf), 'k-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('magnitude mode')
            axs[1].grid('on')

            if(magnitude_mode):
                axs[1].plot(peak_seg_ppm, np.abs(peak_seg_val), 'r-')
            else:
                axs[0].plot(peak_seg_ppm, np.real(peak_seg_val), 'r-')

            fig.subplots_adjust()
            fig.show()

        return(lw)

    def correct_intensity_scaling_nd(self, scaling_factor_rawdata=1e8, scaling_factor_dcm=1.0):
        """
        Amplify the FID signals. Sounds useless but can actually help during quantification! Yes, it is not a good idea to fit signals which have intensities around 1e-6 or lower because of various fit tolerances and also digital problems (epsilon).

        * Works with multi-dimensional signals.

        Parameters
        ----------
        scaling_factor_rawdata : float
            Amplification factor if data is raw
        scaling_factor_dcm : float
            Amplification factor if data is from a dcm file (already amplified)

        Returns
        -------
        s : MRSData2 numpy array [whatever dimensions]
            Resulting amplified MRSData2 object
        """
        log.debug("intensity scaling [%s]..." % self.display_label)

        if(self.is_rawdata):
            scaling_factor = scaling_factor_rawdata
        else:
            scaling_factor = scaling_factor_dcm

        # scale signal
        log.debug("multiplying time-domain signals by %E..." % scaling_factor)
        s_sc = self.copy() * scaling_factor

        # convert back to MRSData2
        s_sc = self.inherit(s_sc)

        # if any ref data available, we crop it too (silently)
        if(s_sc.data_ref is not None):
            s_sc.data_ref = s_sc.data_ref.correct_intensity_scaling_nd(scaling_factor)

        return(s_sc)

    def correct_fidmodulus_nd(self):
        """
        Calculate absolute mode of FID signals. I am not sure I am doing this correctly but it is a first attempt.

        * Works with multi-dimensional signals.

        Returns
        -------
        s : MRSData2 numpy array [whatever dimensions]
            Resulting FID modulus data stored in a MRSData2 object
        """
        # init
        log.debug("fid modulus [%s]..." % self.display_label)
        log.debug("calculation magnitude of signal...")
        s = self.copy()

        # return magnitude
        s = self.inherit(np.abs(s))
        return(s)

    def correct_zerofill_nd(self, nPoints_final=16384, display=False, display_range=[1, 6]):
        """
        Zero-fill MRS data signals along the time axis.

        Parameters
        ----------
        nPoints_final : int
            Final number of points
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        * Works with multi-dimensional signals.
        * Returns a multi-dimensional signal.

        Returns
        -------
        s_zf : MRSData2 numpy array [whatever,...,timepoints]
            Resulting zero-filled data stored in a MRSData2 object
        """
        # init
        log.debug("zero_filling [%s]..." % self.display_label)
        s = self.copy()

        # check
        nZeros = nPoints_final - s.np
        if(nZeros <= 0):
            log.warning("no zero-filling performed. The number of zeros to add was negative (%d)!" % nZeros)
            return(s)

        s_new_shape = list(s.shape)
        s_new_shape[-1] = nZeros
        log.debug("%d-pts signal + %d zeros = %d-pts zero-filled signal..." % (s.np, nZeros, nPoints_final))
        s_zf = self.inherit(np.concatenate((s, np.zeros(s_new_shape)), axis=s.ndim - 1))

        # sequence npts
        if(s_zf.sequence is not None):
            log.debug("updating sequence.npts...")
            s_zf.sequence.npts = nPoints_final
            s_zf.sequence._ready = False

        if(display):
            s_disp = s.copy()
            s_zf_disp = s_zf.copy()
            while(s_disp.ndim > 1):
                s_disp = np.mean(s_disp, axis=0)
                s_zf_disp = np.mean(s_zf_disp, axis=0)

            fig_title = "Zero-filling [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='row', sharey='row')

            # no time axis, we want to see the number of points
            axs[0, 0].plot(np.real(s_disp), 'k-', linewidth=1)
            axs[0, 0].set_xlabel('number of points')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(np.real(s_zf_disp), 'b-', linewidth=1)
            axs[0, 1].set_xlabel('number of points')
            axs[0, 1].set_ylabel('zero-filled')
            axs[0, 1].grid('on')

            axs[1, 0].plot(s_disp.frequency_axis_ppm(), s_disp.spectrum().real, 'k-', linewidth=1)
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('original')
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].grid('on')

            axs[1, 1].plot(s_zf_disp.frequency_axis_ppm(), s_zf_disp.spectrum().real, 'b-', linewidth=1)
            axs[1, 1].set_xlabel("chemical shift (ppm)")
            axs[1, 1].set_ylabel('zero-filled')
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].grid('on')

            fig.subplots_adjust()
            fig.show()

        # if any ref data available, we crop it too (silently)
        if(s_zf.data_ref is not None):
            s_zf.data_ref = s_zf.data_ref.correct_zerofill_nd(nPoints_final, False)

        return(s_zf)

    def correct_time_shift_nd(self, time_shift_us=-375, display=False, display_range=[1, 6]):
        """
        Shift time signals of Zero-fill MRS data signals along the time axis.

        Parameters
        ----------
        time_shift_us : float
            Time shift to apply to time signals. Negative means the beginning of the FIDs will be eaten up and circshifted at the end.
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        * Works with multi-dimensional signals.
        * Returns a multi-dimensional signal.

        Returns
        -------
        s_shifted : MRSData2 numpy array [averages,channels,timepoints]
            Resulting shifted data stored in a MRSData2 object
        """
        log.debug("time_shifting [%s]..." % self.display_label)

        # init
        s = self.copy()

        # prepare frequency vector
        f = s.frequency_axis()
        if(s.ndim > 1):
            f_tiles = list(s.shape)
            f_tiles[-1] = 1
            f = np.tile(f, f_tiles)

        # fft
        sf = np.fft.fftshift(np.fft.fft(s, axis=-1), axes=-1)
        # apply phase 1st order in frequency-frequency domain (which is equivalent to a time shift in time domain)
        sf_shifted = sf * np.exp(1j * 2.0 * np.pi * -time_shift_us / 1000000.0 * f)
        # fft back and prey
        s_shifted = np.fft.ifft(np.fft.ifftshift(sf_shifted, axes=-1), axis=-1)

        # convert back to MRSData2
        s_shifted = self.inherit(s_shifted)

        if(display):
            t = s.time_axis() * 1000000.0  # us
            ppm = s.frequency_axis_ppm()

            s_disp = s.copy()
            s_shifted_disp = s_shifted.copy()
            while(s_disp.ndim > 1):
                s_disp = np.mean(s_disp, axis=0)
                s_shifted_disp = np.mean(s_shifted_disp, axis=0)

            fig_title = "Time-shifting [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='row', sharey='row')

            axs[0, 0].plot(t, np.real(s_disp), 'k-', linewidth=1)
            axs[0, 0].set_xlabel('time (us)')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(t, np.real(s_shifted_disp), 'b-', linewidth=1)
            axs[0, 1].set_xlabel('time (us)')
            axs[0, 1].set_ylabel('time-shifted')
            axs[0, 1].grid('on')

            axs[1, 0].plot(ppm, s_disp.spectrum().real, 'k-', linewidth=1)
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('original')
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].grid('on')

            axs[1, 1].plot(ppm, s_shifted_disp.spectrum().real, 'b-', linewidth=1)
            axs[1, 1].set_xlabel("chemical shift (ppm)")
            axs[1, 1].set_ylabel('time-shifted')
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].grid('on')

            fig.subplots_adjust()
            fig.show()

        # if any ref data available, we crop it too (silently)
        if(s_shifted.data_ref is not None):
            s_shifted.data_ref = s_shifted.data_ref.correct_time_shift_nd(time_shift_us)

        return(s_shifted)

    def correct_phase_3d(self, use_ref_data=True, peak_range=[4.5, 5], average_per_channel_mode=False, first_point_fid_mode=False, phase_order=0, phase_offset=0.0, display=False, display_range=[1, 6]):
        """
        Well, that's a big one but basically it rephases the signal of interest.

        >In the case of multi-channel acquisition, the phase will be estimated for each channel, for each average using phase time evolution estimated on reference signal.
        >If the reference signal is not specified and weak water suppression was performed, 0th order phase correction will be done using the first point in the fid.
        >If strong water suppression was performed, 0th order phase correction will be done using the phase of a peak in the spectrum; the chemical shift range to find this peak is specified by the user.
        >>Note that those last two approaches can be performed for each average in the case of high SNR (rare) or by averaging all the scans for each channel in the case of lower SNR.

        * Works only with a 3D [averages,channels,timepoints] signal.
        * Returns a 3D [averages,channels,timepoints] signal.

        Parameters
        ----------
        use_ref_data : boolean
            Use reference data (usually non water suppressed) for phasing
        peak_range : list [2]
            Range in ppm used to peak-pick and estimate a phase
        average_per_channel_mode : boolean
            Average all the averages for each channel when doing peak analysis
        first_point_fid_mode : boolean
            Estimate phase from 1st point of FID
        phase_order : int
            Order of phasing: 0(th) or 1(st) order phasing
        phase_offset : float
            Phase added to signal (rad)
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_phased : MRSData2 numpy array [averages,channels,timepoints]
            Resulting phased data stored in a MRSData2 object
        """
        log.debug("phasing [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 3):
            log.error("this method only works for 3D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # check if any ref data
        if(self.data_ref is None and use_ref_data):
            log.warning("you want to phase data based on ref. data but no such data is available!")
            use_ref_data = False
        # dimensions check for reference data
        if(use_ref_data and self.data_ref.ndim != 3):
            log.error("this method only works for 3D signals! You are feeding it with %d-dimensional reference data. :s" % self.ndim)

        # init
        s = self.copy()
        s_phased = self.copy()

        # list of phasing methods
        list_phase_method = {}
        list_phase_method[0] = "0th & 1st order phase from ref. scan FID (method #0)"
        list_phase_method[1] = "0th order only phase from ref. scan FID (method #1)"
        list_phase_method[2] = "0th order phase from 1st pt in FID (method #2)"
        list_phase_method[3] = "0th order phase from 1st pt in averaged FID (method #3)"
        list_phase_method[4] = "0th order phase from peak in spectrum (method #4)"
        list_phase_method[5] = "0th order phase from peak in averaged spectrum (method #5)"
        t = s.time_axis()
        ppm = s.frequency_axis_ppm()

        # choose which method is adequate
        only_0th_order = (phase_order == 0)
        if(use_ref_data):
            # we have a ref scan, so method #0 or #1
            phase_method = 0 + only_0th_order
            t_ref = self.data_ref.time_axis()
        else:
            # we do not have a ref scan, so method #2, #3, #4, #5
            # that depends on water intensity and water suppression
            if(first_point_fid_mode):
                phase_method = 2
            else:
                phase_method = 4

            # and on SNR
            if(average_per_channel_mode):
                phase_method += 1

        if(display):
            # prepare subplot
            fig_title = "Phasing [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 3, sharex='col')

        # display chosen method
        log.debug("phasing using method: " + list_phase_method[phase_method])

        # for each channel
        for c in range(0, s.shape[1]):

            if(phase_method == 0 or phase_method == 1):
                # time-domain phase of reference signal for this channel
                sp_ref = np.angle(self.data_ref[0, c, :])

                if(display):
                    # display reference FID
                    axs[0, 0].cla()
                    axs[0, 0].plot(t_ref, np.real(self.data_ref[0, c, :]), linewidth=1, label='real part')
                    axs[0, 0].plot(t_ref, np.imag(self.data_ref[0, c, :]), linewidth=1, label='imag part')
                    axs[0, 0].set_xlabel('time (s)')
                    axs[0, 0].set_ylabel('intensity (u.a)')
                    axs[0, 0].grid('on')
                    axs[0, 0].legend()
                    axs[0, 0].set_title("Ref - channel #" + str(c + 1))

                    # display reference time-domain phase
                    axs[1, 0].cla()
                    axs[1, 0].plot(t_ref, sp_ref, linewidth=1, label='unwrapped phase')
                    axs[1, 0].set_xlabel('time (s)')
                    axs[1, 0].set_ylabel("phase (rd)")
                    axs[1, 0].grid('on')
                    axs[1, 0].legend()
                    axs[1, 0].set_title("Ref - channel #" + str(c + 1))

            elif(phase_method == 3):
                # phase of first point in averaged fid
                s_avg = s[:, c, 0].mean(axis=0)
                phase_fid_avg = np.angle(s_avg)

            elif(phase_method == 5):
                # phase of peak in averaged spectrum
                s_avg = s[:, c, :].mean(axis=0)

                # find maximum peak in range and its phase
                peak_ppm, peak_val, _, _, _ = s_avg._analyze_peak_1d(peak_range)
                phase_peak_avg = np.angle(peak_val)
                if(c == 0):
                    log.debug("measuring phase at %0.2fppm on 1st channel!" % peak_ppm)

            # late init progress bar
            if(c == 0):
                pbar = log.progressbar("phasing", s.shape[1] * s.shape[0])

            # now, for each average in meta signal acquired with this channel
            for a in range(0, s.shape[0]):

                # this spectrum
                this_s = s[a, c, :]
                this_sf = this_s.spectrum()

                if(phase_method == 0):
                    # correct phase using reference time-domain phase estimation
                    this_s_phased = this_s * np.exp(-1j * (sp_ref + phase_offset))
                elif(phase_method == 1):
                    # correct phase using first point of reference time-domain phase estimation
                    this_s_phased = this_s * np.exp(-1j * (sp_ref[0] + phase_offset))
                elif(phase_method == 2):
                    # phase of first point of this fid
                    phase_fid = np.angle(this_s[0])
                    # and apply it to correct the spectrum
                    this_s_phased = this_s * np.exp(-1j * (phase_fid + phase_offset))
                elif(phase_method == 3):
                    # apply to correct the spectrum
                    this_s_phased = this_s * np.exp(-1j * (phase_fid_avg + phase_offset))
                elif(phase_method == 4):
                    # find maximum peak in range and its phase
                    peak_ppm, peak_val, _, _, _ = this_s._analyze_peak_1d(peak_range)
                    phase_peak = np.angle(peak_val)

                    # apply phase to spectrum
                    this_sf_phased = this_sf * np.exp(-1j * (phase_peak + phase_offset))
                    # ifft back
                    this_s_phased = np.fft.ifft(np.fft.ifftshift(this_sf_phased))
                elif(phase_method == 5):
                    # apply phase to spectrum
                    this_sf_phased = this_sf * np.exp(-1j * (phase_peak_avg + phase_offset))
                    # ifft back
                    this_s_phased = np.fft.ifft(np.fft.ifftshift(this_sf_phased))

                # store
                s_phased[a, c, :] = this_s_phased

                if(display):

                    # convert back to MRSData2
                    this_s_phased = self.inherit(this_s_phased)

                    # display original meta FID
                    axs[0, 1].cla()
                    axs[0, 1].plot(t, np.real(this_s), linewidth=1, label='real part')
                    axs[0, 1].plot(t, np.imag(this_s), linewidth=1, label='imag part')
                    axs[0, 1].set_xlabel('time (s)')
                    axs[0, 1].set_ylabel('original')
                    axs[0, 1].grid('on')
                    axs[0, 1].legend()
                    axs[0, 1].set_title("Meta channel #" + str(c + 1) + " average #" + str(a + 1))

                    # display original meta spectrum
                    axs[0, 2].cla()
                    axs[0, 2].plot(ppm, np.real(this_sf), linewidth=1, label='real part')
                    if(phase_method == 4 or phase_method == 5):
                        axs[0, 2].plot(peak_ppm, np.real(peak_val), 'ro')
                        axs[0, 2].axvline(x=peak_ppm, color='r', linestyle='--')
                    axs[0, 2].set_xlim(display_range[1], display_range[0])
                    axs[0, 2].set_xlabel('time (s)')
                    axs[0, 2].set_ylabel('original')
                    axs[0, 2].grid('on')
                    axs[0, 2].legend()
                    axs[0, 2].set_title("Meta channel #" + str(c + 1) + " average #" + str(a + 1))

                    # display corrected meta FID
                    axs[1, 1].cla()
                    axs[1, 1].plot(t, np.real(this_s_phased), 'b-', linewidth=1, label='real part')
                    axs[1, 1].set_xlabel('time (s)')
                    axs[1, 1].set_ylabel('corrected')
                    axs[1, 1].grid('on')
                    axs[1, 1].legend()

                    # display corrected meta spectrum
                    axs[1, 2].cla()
                    axs[1, 2].plot(ppm, this_s_phased.spectrum().real, 'b-', linewidth=1, label='real part')
                    axs[1, 2].set_xlim(display_range[1], display_range[0])
                    axs[1, 2].set_xlabel('frequency (ppm)')
                    axs[1, 2].set_ylabel('corrected')
                    axs[1, 2].grid('on')
                    axs[1, 2].legend()

                    fig.subplots_adjust()
                    fig.show()
                    plt.pause(0.1)

                pbar.update(c * s.shape[0] + a + 1)

        pbar.finish("done")

        # convert back to MRSData2
        s_phased = self.inherit(s_phased)

        # if any ref data available, we phase it too (silently)
        if(s_phased.data_ref is not None):
            s_phased.data_ref = s_phased.data_ref.correct_phase_3d(True)

        return(s_phased)

    def correct_combine_channels_3d(self, use_ref_data=True, phasing=False, channels_onoff=[True]):
        """
        Recombine Rx channels using SVD stuff. If no reference signal is specified, the recombination weights will be calculated from the signal of interest (not optimal).

        * Works only with a 3D [averages,channels,timepoints] signal.
        * Returns a 2D [averages,timepoints] signal.

        Parameters
        ----------
        use_ref_data : boolean
            Use reference data (usually non water suppressed) for channel combining
        phasing : boolean
            Allow 0th order phasing during channel combine or not
        channels_onoff : boolean list [nChannels]
            Binary weights to apply to each channel for example to turn off some of them
        Returns
        -------
        s_combined : MRSData2 numpy array [averages,timepoints]
            Resulting channel combined data stored in a MRSData2 object
        """
        log.debug("channel combining [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 3):
            log.error("this method only works for 3D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # check if any ref data
        if(self.data_ref is None and use_ref_data):
            log.warning("you want to channel-combine data based on ref. data but no such data is available!")
            use_ref_data = False
        # dimensions check for reference data
        if(use_ref_data and self.data_ref.ndim != 3):
            log.error("this method only works for 3D signals! You are feeding it with %d-dimensional reference data. :s" % self.ndim)

        # init
        s = self.copy()
        if(s.shape[1] == 1):
            log.warning("this is a single-channel signal, no need to recombine this!")
            s_combined = np.mean(s, axis=1)
            log.warning("reshaped to " + str(s_combined.shape))
        else:
            if(phasing):
                if(use_ref_data):
                    log.debug("channel recombine WITH reference scan AND phasing (original suspect code)...")
                    weights = suspect.processing.channel_combination.svd_weighting(self.data_ref.mean(axis=0))
                else:
                    log.debug("channel recombine WITHOUT reference scan AND phasing (original suspect code)...")
                    s_dirty_mean = np.mean(s, axis=0)
                    weights = suspect.processing.channel_combination.svd_weighting(s_dirty_mean)
            else:
                if(use_ref_data):
                    log.debug("channel recombine WITH reference scan & NO phasing...")
                    p, _, v = np.linalg.svd(self.data_ref.mean(axis=0), full_matrices=False)
                    channel_weights = p[:, 0].conjugate()
                    weights = -channel_weights / np.sum(np.abs(channel_weights))
                else:
                    log.debug("channel recombine reference scan & NO phasing...")
                    s_dirty_mean = np.mean(s, axis=0)
                    p, _, v = np.linalg.svd(s_dirty_mean, full_matrices=False)
                    channel_weights = p[:, 0].conjugate()
                    weights = -channel_weights / np.sum(np.abs(channel_weights))

            # turn off some channels ?
            channels_onoff_np = np.array(channels_onoff)
            if((channels_onoff_np == False).any()):
                channels_onoff_float = channels_onoff_np.astype(float)
                log.debug("playing with channel weights...")
                log.debug(str(channels_onoff_float))
                weights = weights * channels_onoff_float

            s_combined = suspect.processing.channel_combination.combine_channels(s, weights)

        # convert back to MRSData2
        s_combined = self.inherit(s_combined)

        # if any ref data available, we combine it too (silently)
        if(s_combined.data_ref is not None):
            s_combined.data_ref = s_combined.data_ref.correct_combine_channels_3d(True, True, channels_onoff)

        return(s_combined)

    def concatenate_2d(self, data):
        """
        Concatenate current signal with another one along the averages axis.

        * Works only with a 2D [averages,timepoints] signal.
        * Returns a 2D [averages,timepoints] signal.

        Parameters
        ----------
        data : MRSData2 numpy array [averages,timepoints]
            MRS data to concatenate to the current data

        Returns
        -------
        s_concatenated : MRSData2 numpy array [averages,timepoints]
            Resulting concatenated signal
        """
        log.debug("concatenating [%s] to [%s]..." % (self.display_label, data.display_label))
        # dimensions check
        if(self.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)
        if(data.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % data.ndim)

        # init
        log.debug("concatenating dataset shapes " + str(self.shape) + " and " + str(data.shape) + " ...")

        s_concatenated = np.concatenate((self, data))
        # convert back to MRSData2
        s_concatenated = self.inherit(s_concatenated)
        log.debug("obtained a dataset shape " + str(s_concatenated.shape))

        # update some attributes
        s_concatenated._is_concatenated = True

        return(s_concatenated)

    def _build_moving_average_data_2d(self, nAvgWindow=5):
        """
        Build moving average data in the average dimension. Usefull for the analyze_peak and correct_realign_2d functions.

        * Works only with a 2D [averages,timepoints] signal.
        * Returns a 2D [averages,timepoints] signal.

        Parameters
        ----------
        nAvgWindow : int
            Size of the moving average window

        Returns
        -------
        s_ma : MRSData2 numpy array [averages,timepoints]
            Resulting moving average data stored in a MRSData2 object. The number of averages is the same as the original data BUT each of those average is actually an average of nAvgWindow spectra.
        """
        log.debug("calculating moving average for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        log.debug("moving averaging with window of %d samples!" % nAvgWindow)
        s = self.copy()

        # number of averages for moving average?
        if(np.mod(nAvgWindow, 2) == 0):
            nAvgWindow += 1

        # build moving averaged data
        s_ma = s.copy()
        moving_averages_half = int((nAvgWindow - 1) / 2)
        for a in range(0, s.shape[0]):
            ia = max(0, a - moving_averages_half)
            ib = min(s.shape[0], a + moving_averages_half + 1)
            s_ma[a, :] = np.mean(s[ia:ib, :], axis=0)

        return(s_ma)

    def correct_analyze_and_reject_2d(self, peak_analyze_range=[4.5, 5], peak_snr_range=[1.8, 2.2], peak_lw_range=[4.5, 5], moving_averages=1, reject_when_linewidth_fails=True, peak_properties_ranges={"amplitude (%)": None, "linewidth (Hz)": [5.0, 30.0], "chemical shift (ppm)": 0.5, "phase std. factor (%)": 60.0}, peak_properties_rel2mean=True, auto_method_list=None, auto_adjust_allowed_snr_change=1.0, allowed_apodization=0.0, display=False, display_range=[1, 6]):
        """
        Analyze peak in each average in terms intensity, linewidth, chemical shift and phase and reject data if one of these parameters goes out of the min / max bounds. Usefull to understand what the hell went wrong during your acquisition when you have the raw data and to try to improve things a little. You can choose to set the bounds manually or automatically based on a peak property (amplitude, linewidth, frequency, phase). And you can run several automatic adjusment methods, the one giving the highest SNR and/or the lowest peak linewidth will be selected. All this is very experimental and the code is long and not optimized, sorry ;).
    Special note about the optimization: when does it stop? First, the algorithm tries to optimize the data rejection to get a SNR higher than the (initial SNR * auto_adjust_allowed_snr_change). If the latter is not possible, then the algorithm tries to reduce the linewidth without reducing SNR compared to the initial SNR. If nothing works out, no data rejection is performed except maybe based on peak detection (see reject_when_linewidth_fails).

        * Works only with a 2D [averages,timepoints] signal.
        * Returns a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_analyze_range : list
            Range in ppm used to analyze peak properties (amplitude, linewidth, chemical shift, phase)
        peak_snr_range : list
            Range in ppm used to estimate SNR
        peak_lw_range : list
            Range in ppm used to estimate linewidth
        moving_averages : int
            Number of averages to perform when using moving average, need to be an odd number
        reject_when_linewidth_fails : boolean
            Reject data if linewidth estimatiopn fails (=0Hz). Usually when the peak looks so bad that the peak width cannot be measured...
        peak_properties_ranges : dict
            Dictionnary that contains 4 entries, 4 rejection criterias for
                "amplitude (%)": amplitude relative changes: keep data if within +/-val % range
                "linewidth (Hz)": linewidth changes: keep data is within values in Hz
                "chemical shift (ppm)": chemical shift changes: keep data is within +/-val ppm
                "phase std. factor (%)": phase changes: keep data if within +/- val/100 * std(phase) rd
        peak_properties_rel2mean : boolean
            Relative peak properties (amplitude, chemical shift and phase) should be caculated based on the mean value over the whole acquisition (True) or only the first acquired point (False)
        auto_method_list : list of data_rejection_method
            Automatic rejection bounds adjustment methods
        auto_adjust_allowed_snr_change : float
            Allowed change in SNR (%), a positive or negative relative to the initial SNR without data rejection
        allowed_apodization : float/boolean
            If >0 or !=False, apodize signal during correction process. However, the final corrected signal will not be apodized.
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_cor : MRSData2 numpy array [averages,timepoints]
           Data remaining after data rejection stored in a MRSData2 object.
        """
        log.debug("analyzing data and rejecting some for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        if(s.shape[0] == 1):
            log.warning("single-shot signal, nothing to analyze!")
            return(s)
        ppm = s.frequency_axis_ppm()

        # check if we did data rejection before
        if(self.data_rejection is None):
            iround_data_rej = 1
        else:
            iround_data_rej = len(self.data_rejection) + 1
            log.info("%dth round of data rejection!" % iround_data_rej)

        # estimate initial SNR and linewidth
        log.pause()
        s_avg = s.correct_average_2d()
        initial_snr, _, _ = s_avg.analyze_snr_1d(peak_snr_range)
        initial_lw = s_avg.analyze_linewidth_1d(peak_lw_range)
        log.resume()
        log.info("* Pre-data-rejection SNR = %.2f" % initial_snr)
        log.info("* Pre-data-rejection linewidth = %.2f Hz" % initial_lw)

        # build moving averaged data
        s_ma = self._build_moving_average_data_2d(moving_averages)

        # perform peak analysis (possibly with apodization to stabilize things)
        peak_prop_abs, peak_prop_rel2mean, peak_prop_rel2firstpt = s_ma._analyze_peak_2d(peak_analyze_range, allowed_apodization)

        # first set the data according to relative option: this is a user option
        if(peak_properties_rel2mean):
            peak_prop_rel = peak_prop_rel2mean
        else:
            peak_prop_rel = peak_prop_rel2firstpt

        # choose if absolute or relative will be analyzed: this is hard-coded
        peak_prop_analyze = peak_prop_abs * 0.0
        # amplitude: relative in %
        peak_prop_analyze[:, 0] = peak_prop_rel[:, 0]
        # linewidth: absolute in Hz
        peak_prop_analyze[:, 1] = peak_prop_abs[:, 1]
        # frequency: relative in ppm
        peak_prop_analyze[:, 2] = peak_prop_rel[:, 2]
        # phase: absolute in rad
        peak_prop_analyze[:, 3] = peak_prop_rel[:, 3]

        # choose if absolute or relative will be displayed
        peak_prop_disp = peak_prop_rel * 0.0
        # amplitude: relative in %
        peak_prop_disp[:, 0] = peak_prop_rel[:, 0]
        # linewidth: absolute in Hz
        peak_prop_disp[:, 1] = peak_prop_abs[:, 1]
        # frequency: absolute in ppm
        peak_prop_disp[:, 2] = peak_prop_abs[:, 2]
        # phase: absolute in rad
        peak_prop_disp[:, 3] = peak_prop_abs[:, 3]

        # our time scale
        t_ma = np.linspace(0, self.tr * s.shape[0], s_ma.shape[0]) / 1000.0  # s

        # stats
        log.info("peak analysis: means ± std. deviations")
        log.info("rel. peak amplitude = %.2f ± %.2f %%" % (peak_prop_disp[:, 0].mean(), peak_prop_disp[:, 0].std()))
        log.info("abs. linewidth = %.1f ± %.1f Hz (%.3f ± %.3f ppm)" % (peak_prop_disp[:, 1].mean(), peak_prop_disp[:, 1].std(), (peak_prop_disp[:, 1] / s_ma.f0).mean(), (peak_prop_disp[:, 1] / s_ma.f0).std()))
        log.info("abs. frequency = %.2f ± %.2f ppm (± %.1f Hz)" % (peak_prop_disp[:, 2].mean(), peak_prop_disp[:, 2].std(), (peak_prop_disp[:, 2] * s_ma.f0).std()))
        log.info("abs. phase = %.2f ± %.2f rad" % (peak_prop_disp[:, 3].mean(), peak_prop_disp[:, 3].std()))

        # check for Nones
        peak_properties_ranges_list = list(peak_properties_ranges.values())
        peak_properties_ranges_list = [np.inf if p is None else p for p in peak_properties_ranges_list]

        # special for linewidth: can be a max linewidth or a list
        if(type(peak_properties_ranges_list[1]) != list):
            peak_properties_ranges_list[1] = [1.0, peak_properties_ranges_list[1]]

        # reject when peak linewidth fails
        if(reject_when_linewidth_fails):
            peak_properties_ranges_list[1][0] = 1.0
        else:
            peak_properties_ranges_list[1][0] = -1.0

        # special for phase: rejection range is a factor of std
        phase_std = peak_prop_analyze[:, 3].std()
        phase_std_reject_range = peak_properties_ranges_list[3] / 100.0 * phase_std

        # prepare rejection min/max vectors
        peak_prop_min = [-peak_properties_ranges_list[0],
                         peak_properties_ranges_list[1][0],
                         -peak_properties_ranges_list[2],
                         -phase_std_reject_range]

        peak_prop_max = [+peak_properties_ranges_list[0],
                         peak_properties_ranges_list[1][1],
                         +peak_properties_ranges_list[2],
                         +phase_std_reject_range]

        # automatic rejection ?
        if(auto_method_list is not None):

            # init
            display_axes_ready = [False, False]
            properties_names = list(peak_properties_ranges.keys())
            auto_method_final_snr_list = np.array([0.0] * 4)
            auto_method_final_lw_list = np.array([np.inf] * 4)
            peak_prop_min_auto_res = peak_prop_min.copy()
            peak_prop_max_auto_res = peak_prop_max.copy()

            # for each auto method
            for this_auto_method in auto_method_list:
                # prepare min/max peak property range
                this_prop_min = np.abs(peak_prop_analyze[:, this_auto_method.value]).min()
                this_prop_max = np.abs(peak_prop_analyze[:, this_auto_method.value]).max()

                # get range resolution from constants
                if(this_auto_method == data_rejection_method.AUTO_AMPLITUDE):
                    this_prop_step = DATA_REJECTION_AMPLITUDE_STEP
                elif(this_auto_method == data_rejection_method.AUTO_LINEWIDTH):
                    this_prop_step = DATA_REJECTION_LINEWIDTH_STEP
                elif(this_auto_method == data_rejection_method.AUTO_FREQUENCY):
                    this_prop_step = DATA_REJECTION_FREQUENCY_STEP
                elif(this_auto_method == data_rejection_method.AUTO_PHASE):
                    this_prop_step = DATA_REJECTION_PHASE_STEP
                else:
                    log.error("upsyy! I am not aware of this automatic data rejection method: " + str(this_auto_method))

                # generate a range to test, using the resolution
                this_prop_range = np.arange(this_prop_min - this_prop_step, this_prop_max + this_prop_step, this_prop_step)

                # checking that there is actually a variation and a range to test
                if(this_prop_range.size == 0):
                    # no? let's skip this method
                    this_prop_range = np.array([this_prop_min])
                else:
                    # now let's be smart and reduce the number of values to test according to the peak property measurements
                    # regrid to nearest in previously computed range
                    this_prop_analyze_set = np.abs(peak_prop_analyze[:, this_auto_method.value])
                    this_prop_analyze_set = interpolate.interp1d(this_prop_range, this_prop_range, kind='nearest')(this_prop_analyze_set)
                    # remove one step resolution, remove duplicates and sort
                    this_prop_analyze_set = np.sort(np.array(list(set(this_prop_analyze_set - this_prop_step))))
                    # remove negative values
                    this_prop_analyze_set = this_prop_analyze_set[this_prop_analyze_set > 0]
                    # we should now have an optimized set of thresholds to test
                    this_prop_range = this_prop_analyze_set

                    # checking that there is actually a variation and a range to test (again, sorry)
                    if(this_prop_range.size == 0):
                        # no? let's skip this method
                        this_prop_range = np.array([this_prop_min])

                # iterate and test the resulting data
                pbar = log.progressbar("adjusting rejection threshold for [" + properties_names[this_auto_method.value] + "] in range [%.3f;%.3f] (n=%d)" % (this_prop_range.min(), this_prop_range.max(), this_prop_range.size), this_prop_range.shape[0])

                # add inf to be sure that we try without any thresholds
                this_prop_range = np.hstack((this_prop_range, +np.inf))

                test_snr_list = np.zeros(this_prop_range.shape)
                test_lw_list = np.zeros(this_prop_range.shape)
                test_nrej_list = np.zeros(this_prop_range.shape)

                # test each criteria bound
                for (i_prop_val, this_prop_val) in enumerate(this_prop_range):
                    # rebuild min/max rejection bounds including user values
                    peak_prop_min_auto = peak_prop_min.copy()
                    peak_prop_max_auto = peak_prop_max.copy()
                    peak_prop_max_auto[this_auto_method.value] = this_prop_val

                    # now see what we can reject
                    this_mask_reject_data = np.full([s_ma.shape[0], 4], False)
                    for a in range(0, s_ma.shape[0]):
                        for p in range(4):
                            if(peak_prop_analyze[a, p] < peak_prop_min_auto[p]):
                                this_mask_reject_data[a, p] = True
                            if(peak_prop_analyze[a, p] > peak_prop_max_auto[p]):
                                this_mask_reject_data[a, p] = True

                    # reject data now
                    this_mask_reject_data_sumup = (this_mask_reject_data.sum(axis=1) > 0)
                    this_s_cor = s[(this_mask_reject_data_sumup == False), :]

                    # analyze snr / lw and number of rejections
                    if(this_mask_reject_data_sumup.sum() < s_ma.shape[0]):
                        log.pause()
                        this_s_cor_avg = this_s_cor.correct_average_2d()
                        test_snr_list[i_prop_val], _, _ = this_s_cor_avg.analyze_snr_1d(peak_snr_range)
                        test_lw_list[i_prop_val] = this_s_cor_avg.analyze_linewidth_1d(peak_lw_range)
                        log.resume()
                        test_nrej_list[i_prop_val] = this_mask_reject_data_sumup.sum()

                    # progression
                    pbar.update(i_prop_val)

                pbar.finish("done")

                # relative SNR and minimum acceptable relative SNR change
                test_snr_threshold = initial_snr + initial_snr * auto_adjust_allowed_snr_change / 100.0
                test_snr_list_rel = test_snr_list / initial_snr * 100.0 - 100.0

                # relative LW change
                test_lw_list_rel = test_lw_list - initial_lw

                # first, try and find a higher SNR than the initial one (best case, we reject crappy data and improved final SNR)
                if(test_snr_list_rel.max() > auto_adjust_allowed_snr_change):
                    log.info("SNR change above threshold: %.2f%% > %.2f%% threshold! :)" % (test_snr_list_rel.max(), auto_adjust_allowed_snr_change))
                    ind_max_snr = np.argmax(test_snr_list)
                    optim_prop = this_prop_range[ind_max_snr]
                    optim_res_snr = test_snr_list[ind_max_snr]
                    optim_res_lw = test_lw_list[ind_max_snr]
                    log.info("optimal [" + properties_names[this_auto_method.value] + "] = %.1f" % optim_prop)
                else:
                    # no SNR above threshold found, so let's try to find at least a lower linewidth (intermediate case), for the same SNR or more
                    # check that we have a segment of the curve above the initial SNR
                    test_snr_list_mask = (test_snr_list_rel >= 0.0)
                    if(not test_snr_list_mask.any()):
                        # that was a bit ambitious, there was absolutely no SNR enhancement
                        log.info("sorry, this is only making your SNR worse...")
                        log.debug("the best SNR change we found was %.2f%% compared to initial :(" % test_snr_list_rel.max())
                        # set optimal LW to max (inf)
                        optim_prop = this_prop_max
                        optim_res_snr = test_snr_list[-1]
                        optim_res_lw = test_lw_list[-1]
                    else:
                        # we found relative SNR changes equal (no change) or above 0 (little SNR enchancement, still below expectation)
                        # let's choose the one with the lowest LW
                        min_lw_snr_masked = np.min(test_lw_list[test_snr_list_mask])
                        ind_min_lw_snr_masked = np.argmin(test_lw_list[test_snr_list_mask])
                        # if LW was actually reduced
                        if(min_lw_snr_masked < 0):
                            log.info("could not improve SNR above threshold but reduced peak linewidth! :)")
                            optim_prop = this_prop_range[test_snr_list_mask][ind_min_lw_snr_masked]
                            optim_res_snr = test_snr_list[test_snr_list_mask][ind_min_lw_snr_masked]
                            optim_res_lw = test_lw_list[test_snr_list_mask][ind_min_lw_snr_masked]
                            log.info("optimal [" + properties_names[this_auto_method.value] + "] = %.1f" % optim_prop)
                        else:
                            # that was a bit ambitious, there was absolutly no SNR enhancement and no LW enhancement too! :(
                            log.info("sorry, this is only making your SNR and LW worse...")
                            # set optimal LW to max (inf)
                            optim_prop = this_prop_max
                            optim_res_snr = test_snr_list[-1]
                            optim_res_lw = test_lw_list[-1]

                # display and save the final snr and lw
                log.info("* Post-data-rejection based on [" + properties_names[this_auto_method.value] + "] SNR = %.2f" % optim_res_snr)
                log.info("* Post-data-rejection based on [" + properties_names[this_auto_method.value] + "] linewidth = %.2f Hz" % optim_res_lw)
                auto_method_final_snr_list[this_auto_method.value] = optim_res_snr
                auto_method_final_lw_list[this_auto_method.value] = optim_res_lw

                # plot SNR / LW combinaisons and optimal choice
                if(display):
                    # plot the SNRs versus LWs
                    fig_title = "Data discarding [%s]: adjusting criteria" % self.display_label
                    if(iround_data_rej > 1):
                        fig_title += " (round #%d)" % iround_data_rej

                    fig = plt.figure(fig_title)
                    if(not display_axes_ready[0]):
                        # create the figure, let's create the axes
                        fig.clf()
                        fig.suptitle(fig_title)
                        fig.subplots(2, 2)
                        for a in fig.axes:
                            a.twinx()
                        display_axes_ready[0] = True

                    fig.axes[this_auto_method.value].plot(this_prop_range, test_snr_list, 'rx-', label='SNR')
                    fig.axes[this_auto_method.value].axvline(optim_prop, color='m', linestyle='--', label='Optimal')
                    fig.axes[this_auto_method.value].axhline(test_snr_threshold, color='g', linestyle='--', label='SNR threshold')
                    fig.axes[this_auto_method.value].set_xlabel(properties_names[this_auto_method.value][0].upper() + properties_names[this_auto_method.value][1:])
                    fig.axes[this_auto_method.value].set_ylabel('Estimated SNR (u.a)')
                    fig.axes[this_auto_method.value].grid('on')
                    fig.axes[this_auto_method.value].legend(loc='lower left')

                    fig.axes[this_auto_method.value + 4].plot(this_prop_range, test_lw_list, 'bx-', label='Linewidth')
                    fig.axes[this_auto_method.value + 4].set_ylabel('Estimated linewidth (Hz)')
                    fig.axes[this_auto_method.value + 4].legend(loc='lower right')

                    fig.subplots_adjust()
                    fig.show()

                # save the found optimal criteria to rejection critera vector
                if(this_auto_method == data_rejection_method.AUTO_LINEWIDTH):
                    peak_prop_max_auto_res[this_auto_method.value] = optim_prop
                else:
                    peak_prop_min_auto_res[this_auto_method.value] = -optim_prop
                    peak_prop_max_auto_res[this_auto_method.value] = optim_prop

            # tried all auto methods, now apply best one
            # the final snr and lw obtained are in auto_method_final_snr_list and auto_method_final_lw_list
            # the final bounds are in peak_prop_min_auto_res and peak_prop_max_auto_res

            # relative SNR change
            auto_method_final_snr_list_rel = auto_method_final_snr_list / initial_snr * 100.0 - 100.0
            auto_method_final_lw_list_rel = auto_method_final_lw_list - initial_lw

            # is this higher than the initial snr?
            if(auto_method_final_snr_list_rel.max() > auto_adjust_allowed_snr_change):
                # apply this method
                ind_max_snr_auto_method = np.argmax(auto_method_final_snr_list)
                optim_auto_method = data_rejection_method(ind_max_snr_auto_method)
                log.info("best adjustment done with " + str(optim_auto_method) + " regarding SNR! (round #%d)" % iround_data_rej)
            else:
                # no methods could give a SNR above threshold
                # so let's try to find at least a method that does not reduce SNR and gives a lower linewidth (intermediate case)
                auto_method_final_snr_list_mask = (auto_method_final_snr_list_rel >= 0.0)
                if(not auto_method_final_snr_list_mask.any()):
                    # that was a bit ambitious
                    log.info("sorry but only made your SNR worse...")
                    log.debug("the best SNR change we found was %.2f%% compared to initial :(" % auto_method_final_snr_list_rel.max())
                    log.info("automatic data rejection failed, no optimal method found, sorry! :(")
                    # no optimal method!
                    optim_auto_method = None
                else:
                    # we found relative SNR changes equal (no change) or above 0 (little SNR enchancement, still below expectation)
                    # let's choose the method that gives the the lowest LW
                    min_lw_auto_method = np.min(auto_method_final_lw_list_rel[auto_method_final_snr_list_mask])
                    ind_min_lw_auto_method = np.argmin(auto_method_final_lw_list_rel[auto_method_final_snr_list_mask])

                    # if LW was actually reduced
                    if(min_lw_auto_method < 0):
                        log.info("could not find a method that improves SNR above threshold, only reduces the peak linewidth! :)")
                        optim_auto_method = data_rejection_method(ind_min_lw_auto_method)
                        log.info("best adjustment done with " + str(optim_auto_method) + " regarding linewidth! (round #%d)" % iround_data_rej)
                    else:
                        # that was a bit ambitious, there was absolutly no SNR enhancement and no LW enhancement too! :(
                        log.info("automatic data rejection failed, no optimal method found, sorry! :(")
                        # no optimal method!
                        optim_auto_method = None

            if(optim_auto_method is not None):
                log.info("* Post-data-rejection SNR = %.2f" % auto_method_final_snr_list[optim_auto_method.value])
                log.info("* Post-data-rejection linewidth = %.2f Hz" % auto_method_final_lw_list[optim_auto_method.value])

                # apply automatically optimized bounds to rejection vectors
                peak_prop_min[optim_auto_method.value] = peak_prop_min_auto_res[optim_auto_method.value]
                peak_prop_max[optim_auto_method.value] = peak_prop_max_auto_res[optim_auto_method.value]

        else:
            optim_auto_method = None

        # for each average, check if peak parameters are in the min / max bounds
        mask_reject_data = np.full([s_ma.shape[0], 4], False)
        pbar = log.progressbar("rejecting data", s_ma.shape[0])
        for a in range(0, s_ma.shape[0]):
            for p in range(4):
                if(peak_prop_analyze[a, p] < peak_prop_min[p]):
                    mask_reject_data[a, p] = True
                if(peak_prop_analyze[a, p] > peak_prop_max[p]):
                    mask_reject_data[a, p] = True

            pbar.update(a)

        pbar.finish("done")

        # stats regarding data rejection, how many, for what reasons, overall percentage
        log.info("data rejection: summary (round #%d)" % iround_data_rej)
        log.info("number of averages rejected because of...")
        log.info("amplitude = %d" % mask_reject_data[:, 0].sum())
        log.info("linewidth = %d" % mask_reject_data[:, 1].sum())
        log.info("frequency = %d" % mask_reject_data[:, 2].sum())
        log.info("phase = %d" % mask_reject_data[:, 3].sum())

        # actually reject data now
        mask_reject_data_sumup = (mask_reject_data.sum(axis=1) > 0)
        s_cor = s[(mask_reject_data_sumup == False), :]
        # build rejected spectrum
        s_rej = s[(mask_reject_data_sumup == True), :]
        log.pause()
        s_rej_avg = s_rej.correct_average_2d()
        log.resume()

        log.info("TOTAL data rejection = %d / %d (%.0f%%)" % (mask_reject_data_sumup.sum(), s_ma.shape[0], (mask_reject_data_sumup.sum() / s_ma.shape[0] * 100)))

        # perform post-correction measurements
        peak_prop_abs, peak_prop_rel2mean, peak_prop_rel2firstpt = s_ma._analyze_peak_2d(peak_analyze_range, allowed_apodization)

        # first set the data according to relative option: this is a user option
        if(peak_properties_rel2mean):
            peak_prop_rel = peak_prop_rel2mean
        else:
            peak_prop_rel = peak_prop_rel2firstpt

        # choose if absolute or relative will be analyzed: this is hard-coded
        peak_prop_analyze_postcor = peak_prop_abs * 0.0
        # amplitude: relative in %
        peak_prop_analyze_postcor[:, 0] = peak_prop_rel[:, 0]
        # linewidth: absolute in Hz
        peak_prop_analyze_postcor[:, 1] = peak_prop_abs[:, 1]
        # frequency: relative in ppm
        peak_prop_analyze_postcor[:, 2] = peak_prop_rel[:, 2]
        # phase: absolute in rad
        peak_prop_analyze_postcor[:, 3] = peak_prop_rel[:, 3]

        # final display
        if(display):
            fig_title = "Data discarding [%s]: summary" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 3, sharex='all')

            k = 0
            for ix in range(2):
                for iy in range(2):
                    # original data
                    axs[ix, iy].plot(t_ma, peak_prop_analyze[:, k], 'k-x', linewidth=1)
                    # rejected data
                    t_ma_rej = t_ma[mask_reject_data[:, k]]
                    this_peak_prop_analyze_rej = peak_prop_analyze[mask_reject_data[:, k], k]
                    axs[ix, iy].plot(t_ma_rej, this_peak_prop_analyze_rej, 'ro', linewidth=1)
                    axs[ix, iy].axhline(y=peak_prop_min[k], color='r', linestyle='--')
                    axs[ix, iy].axhline(y=peak_prop_max[k], color='r', linestyle='--')
                    k = k + 1

            axs[0, 0].set_ylabel('Rel. amplitude change (%)')
            axs[0, 0].grid('on')
            axs[0, 0].set_title("Rel. amplitude = %.2f ± %.2f %%" % (peak_prop_disp[:, 0].mean(), peak_prop_disp[:, 0].std()))

            axs[0, 1].set_ylabel('Abs. linewidth (Hz)')
            axs[0, 1].grid('on')
            axs[0, 1].set_title("Abs. linewidth = %.1f ± %.1f Hz (%.3f ± %.3f ppm)" % (peak_prop_disp[:, 1].mean(), peak_prop_disp[:, 1].std(), (peak_prop_disp[:, 1] / s_ma.f0).mean(), (peak_prop_disp[:, 1] / s_ma.f0).std()))

            axs[1, 0].set_xlabel('Acq. time (s)')
            axs[1, 0].set_ylabel('Abs. frequency (ppm)')
            axs[1, 0].grid('on')
            axs[1, 0].set_title("Abs. frequency = %.2f ± %.2f ppm (± %.1f Hz)" % (peak_prop_disp[:, 2].mean(), peak_prop_disp[:, 2].std(), (peak_prop_disp[:, 2] * s_ma.f0).std()))

            axs[1, 1].set_xlabel('Acq. time (s)')
            axs[1, 1].set_ylabel('Abs. phase shift (rd)')
            axs[1, 1].grid('on')
            axs[1, 1].set_title("Abs. phase = %.2f ± %.2f rad" % (peak_prop_disp[:, 3].mean(), peak_prop_disp[:, 3].std()))

            # nice plot showing all raw data
            ax = plt.subplot(1, 3, 3)
            ppm = s_ma.frequency_axis_ppm()
            ystep = np.max(np.mean(s_ma.spectrum().real, axis=0))
            ystep = np.power(10, 1 + (np.floor(np.log10(ystep))))

            ampfactor = 4
            for k in range(s_ma.shape[0]):
                if(mask_reject_data_sumup[k]):
                    plt.plot(ppm, s_ma[k, :].spectrum().real * ampfactor + ystep * k, 'r-', linewidth=1)
                else:
                    plt.plot(ppm, s_ma[k, :].spectrum().real * ampfactor + ystep * k, 'g-', linewidth=1)

                # build lineshape segment
                _, _, _, peak_seg_ppm, peak_seg_val = s_ma[k, :]._analyze_peak_1d(peak_analyze_range)
                plt.plot(peak_seg_ppm, np.real(peak_seg_val) * ampfactor + ystep * k, 'k-', linewidth=1)

            plt.xlim(peak_analyze_range[1], peak_analyze_range[0])
            plt.xlabel('chemical shift (ppm)')
            plt.ylabel('individual spectra')

            # y ticks: spectrum index
            # TODO: maybe need to calculate this automatically
            n_yticks = 16
            step_yticks = 2 ** np.ceil(np.sqrt(s.shape[0] / n_yticks))
            yt_lbl_list = np.arange(0, s.shape[0], step_yticks).tolist()
            if(yt_lbl_list[-1] != s.shape[0]):
                yt_lbl_list = yt_lbl_list + [s.shape[0]]
            yt_lbl_list = [("%d" % yt) for yt in yt_lbl_list]

            yt_loc_list = np.arange(0, s.shape[0] * ystep, step_yticks * ystep).tolist()
            if(yt_loc_list[-1] != (s.shape[0] * ystep)):
                yt_loc_list = yt_loc_list + [s.shape[0] * ystep]

            ax.set_yticks(yt_loc_list)
            ax.set_yticklabels(labels=yt_lbl_list)

            plt.grid('on')
            fig.subplots_adjust()
            fig.show()

        # wait, are we removing all data ???
        if(mask_reject_data_sumup.sum() == s.shape[0]):
            log.error("all data is rejected! You need to readjust your rejection bounds...")

        # estimate final SNR and linewidth
        log.pause()
        s_cor_avg = s_cor.correct_average_2d()
        final_snr, _, _ = s_cor_avg.analyze_snr_1d(peak_snr_range)
        final_lw = s_cor_avg.analyze_linewidth_1d(peak_lw_range)
        log.resume()
        log.info("* Final post-data-rejection SNR = %.2f" % final_snr)
        log.info("* Final post-data-rejection linewidth = %.2f Hz" % final_lw)

        # fill up dict about this data rejection
        data_rej_dict = {}
        data_rej_dict["Pre-rejection"] = {}
        data_rej_dict["Pre-rejection"]["snr"] = initial_snr
        data_rej_dict["Pre-rejection"]["lw"] = initial_lw
        data_rej_dict["Pre-rejection"]["na"] = s.shape[0]
        data_rej_dict["Pre-rejection"]["measurements"] = peak_prop_analyze
        data_rej_dict["Post-rejection"] = {}
        data_rej_dict["Post-rejection"]["snr"] = final_snr
        data_rej_dict["Post-rejection"]["lw"] = final_lw
        data_rej_dict["Post-rejection"]["na"] = s_cor.shape[0]
        data_rej_dict["Post-rejection"]["measurements"] = peak_prop_analyze_postcor
        # final rejection bounds
        final_peak_properties_ranges = peak_properties_ranges.copy()
        final_peak_properties_ranges["amplitude (%)"] = np.abs(peak_prop_max[0])
        final_peak_properties_ranges["linewidth (Hz)"] = [peak_prop_min[1], peak_prop_max[1]]
        final_peak_properties_ranges["chemical shift (ppm)"] = np.abs(peak_prop_max[2])
        final_peak_properties_ranges["phase std. factor (%)"] = np.abs(peak_prop_max[3]) / phase_std * 100.0  # special for phase
        data_rej_dict["Rejection bounds"] = final_peak_properties_ranges
        # auto methods
        data_rej_dict["Automatic data rejection methods"] = {}
        data_rej_dict["Automatic data rejection methods"]["Methods tried"] = auto_method_list
        data_rej_dict["Automatic data rejection methods"]["Best method"] = optim_auto_method
        data_rej_dict["Automatic data rejection methods"]["SNR change threshold (%)"] = auto_adjust_allowed_snr_change

        # check if empty or not (if first data rejection or not)
        if(s_cor._data_rejection is None):
            s_cor._data_rejection = [data_rej_dict]
        else:
            s_cor._data_rejection.append(data_rej_dict)

        return(s_cor)

    def correct_realign_2d(self, peak_range=[4.5, 5], moving_averages=1, inter_corr_mode=False, freq_shift_max=25, allowed_apodization=5.0, display=False, display_range=[1, 6]):
        """
        Realign each signal of interest in frequency by taking as a reference the first spectra in absolute mode using pick-picking or inter-correlation (experimental).

        * Works only with a 2D [averages,timepoints] signal.
        * Returns a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to analyze peak phase
        moving_averages : int
            Number of averages to perform when using moving average, need to be an odd number
        inter_corr_mode : boolean
            Use inter-correlation technique to adjust frequency shifts. Could be more robust when SNR is low.
        freq_shift_max : float
            Max allowed frequency shift during realignment (Hz).
        allowed_apodization : float/boolean
            If >0 or !=False, apodize signal during correction process. However, the final corrected signal will not be apodized.
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_realigned : MRSData2 numpy array [averages,timepoints]
            Resulting frequency realigned data stored in a MRSData2 object
        """
        log.debug("frequency realigning [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()
        s_realigned = s.copy()

        if(s.shape[0] == 1):
            log.warning("single-shot signal, cannot realign this!")
        else:
            # build moving averaged data
            s_ma = self._build_moving_average_data_2d(moving_averages)

            # init
            s_avg = np.mean(s, axis=0)
            if(inter_corr_mode):
                # let's fix a +/-0.5ppm range
                f_shifts_min = - np.abs(peak_range[1] - peak_range[0]) * s_ma.f0
                f_shifts_max = + np.abs(peak_range[1] - peak_range[0]) * s_ma.f0
                # let's fix a the resolution here for inter-correlation tests
                f_shifts_step = RECO_CORRECT_REALIGN_INTER_CORR_MODE_DF * s_ma.f0
                f_shifts_list = np.arange(f_shifts_min, f_shifts_max, f_shifts_step)
            else:
                # find peak in average spectrum absolute mode
                ppm_peak_avg, peak_val, _, _, _ = s_avg._analyze_peak_1d(peak_range, allowed_apodization)
                log.debug("measuring peak properties at %0.2fppm!" % ppm_peak_avg)

            # for each average in moving averaged data
            s_realigned_ma = s_ma.copy()
            df_trace = np.zeros(s_ma.shape[0])
            pbar = log.progressbar("realigning", s_ma.shape[0])
            for a in range(0, s_ma.shape[0]):

                if(inter_corr_mode):
                    # compare this individual spectrum with the first,  using inter-correlation

                    # zero-fill and apodize moving average signal if needed
                    # btw, I only do this in the case of the intercorrelation mode because it is done internally for the peak-picking mode by the the method _analyze_peak_1d
                    log.pause()
                    s_ma_ic = s_ma.correct_zerofill_nd().correct_apodization_nd(allowed_apodization)
                    log.resume()

                    # first spectrum as reference
                    s_ma_ic_ref = np.abs(s_ma_ic[0, :].spectrum())
                    # use the peak_range as a range for inter-corr tests
                    cc_2d = f_shifts_list * 0.0
                    for ifs, fs in enumerate(f_shifts_list):
                        s_ma_ic_shifted = np.abs(s_ma_ic[a, :].adjust_frequency(fs).spectrum())
                        cc = np.corrcoef(s_ma_ic_ref, s_ma_ic_shifted)
                        cc_2d[ifs] = np.abs(cc[0, 1])

                    # find max correlation
                    optimal_fs_ind = np.argmax(cc_2d)
                    optimal_fs = f_shifts_list[optimal_fs_ind]
                    df_trace[a] = optimal_fs
                else:
                    # measure shift on moving average data
                    ppm_peak, _, _, _, _ = s_ma[a, :]._analyze_peak_1d(peak_range, allowed_apodization)

                    # estimate frequency shift in Hz compared to average spectrum
                    dppm = -(ppm_peak_avg - ppm_peak)
                    df_trace[a] = dppm * s_ma.f0

                # check max shift
                if(np.abs(df_trace[a]) > freq_shift_max):
                    # that is too much, do not realign this spectrum
                    df_trace[a] = 0.0

                # correct moving averaged data (for display only, less heavy)
                s_realigned_ma[a, :] = s_ma[a, :].adjust_frequency(df_trace[a])
                # correct original data
                s_realigned[a, :] = s[a, :].adjust_frequency(df_trace[a])

                pbar.update(a)

            pbar.finish("done")

            # final display
            if(display):
                fig_title = "Realigning individual spectra [%s]" % self.display_label
                fig = plt.figure(fig_title)
                fig.clf()
                fig.suptitle(fig_title)
                axs = fig.subplots(2, 3, sharex='all', sharey='all')

                # display original averaged spectrum
                axs[0, 0].plot(ppm, np.abs(s_avg.spectrum()), 'k-', linewidth=1)
                axs[0, 0].set_xlim(display_range[1], display_range[0])
                axs[0, 0].set_ylabel('averaged')
                axs[0, 0].grid('on')
                if(not inter_corr_mode):
                    # add peak position
                    axs[0, 0].plot(ppm_peak_avg, np.abs(peak_val), 'ro')
                    axs[0, 0].axvline(x=ppm_peak_avg, color='r', linestyle='--')

                # display original data
                axs[0, 1].plot(ppm, np.abs(s_ma.spectrum().transpose()), 'k-', linewidth=1)
                axs[0, 1].set_xlim(display_range[1], display_range[0])
                axs[0, 1].set_ylabel('original')
                axs[0, 1].grid('on')

                # display corrected spectra
                axs[1, 1].plot(s_realigned_ma.frequency_axis_ppm(), np.abs(s_realigned_ma.spectrum().transpose()), 'b-', linewidth=1)
                axs[1, 1].set_xlim(display_range[1], display_range[0])
                axs[1, 1].set_xlabel('chemical shift (ppm)')
                axs[1, 1].set_ylabel('corrected')
                axs[1, 1].grid('on')

                # display corrected averaged spectrum
                axs[1, 0].plot(s_realigned_ma.frequency_axis_ppm(), np.abs(np.mean(s_realigned_ma, axis=0).spectrum().transpose()), 'b-', linewidth=1)
                axs[1, 0].set_xlim(display_range[1], display_range[0])
                axs[1, 0].set_xlabel('chemical shift (ppm)')
                axs[1, 0].set_ylabel('averaged & corrected')
                axs[1, 0].grid('on')

                plt.subplot(1, 3, 3)
                plt.plot(df_trace, np.arange(s_ma.shape[0]), 'k-x', linewidth=1)
                plt.xlabel('estimated frequency shift (Hz)')
                plt.ylabel('average index')
                plt.grid('on')
                fig.subplots_adjust()
                fig.show()

        # convert back to MRSData2
        s_realigned = self.inherit(s_realigned)

        return(s_realigned)

    def correct_average_2d(self, na=None, display=False, display_range=[1, 6]):
        """
        Average all averages data into one 1D MRS signal.

        * Works only with a 2D [averages,timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        na : int
            Number of signals to average
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_mean : MRSData2 numpy array [timepoints]
            Resulting frequency realigned data stored in a MRSData2 object
        """
        log.debug("averaging [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()

        if(s.shape[0] == 1):
            log.warning("single-shot signal, nothing to average!")
            s_mean = np.mean(s, axis=0)
            log.warning("reshaped to a " + str(s_mean.shape) + " vector")
        else:
            log.debug("averaging data...")

            if(na is not None):
                log.debug("only " + str(na) + "...")
                if(na == 1):
                    s_mean = s[0, :]
                else:
                    s_mean = np.mean(s[0:(na - 1), :], axis=0)
            else:
                s_mean = np.mean(s, axis=0)

            if(display):
                ppm = s.frequency_axis_ppm()
                ppm_mean = s_mean.frequency_axis_ppm()

                fig_title = "Averaging [%s]" % self.display_label
                fig = plt.figure(fig_title)
                fig.clf()
                fig.suptitle(fig_title)
                axs = fig.subplots(2, 1, sharex='all', sharey='all')

                axs[0].plot(ppm, s.spectrum().real.transpose(), 'k-', linewidth=1)
                axs[0].set_xlim(display_range[1], display_range[0])
                axs[0].set_xlabel('chemical shift (ppm)')
                axs[0].set_ylabel('all spectra')
                axs[0].grid('on')

                axs[1].plot(ppm_mean, s_mean.spectrum().real.transpose(), 'b-', linewidth=1)
                axs[1].set_xlim(display_range[1], display_range[0])
                axs[1].set_xlabel('chemical shift (ppm)')
                axs[1].set_ylabel('averaged spectrum')
                axs[1].grid('on')

                fig.subplots_adjust()
                fig.show()

        # convert back to MRSData2
        s_mean = self.inherit(s_mean)

        # if any ref data available, we average it too (silently)
        if(s_mean.data_ref is not None):
            s_mean.data_ref = s_mean.data_ref.correct_average_2d(None, False)

        return(s_mean)

    def correct_phase_1d(self, suspect_method=suspect_phasing_method.MATCH_MAGNITUDE_REAL, ppm_range=[0, 6], allowed_apodization=1.0, display=False, display_range=[1, 6]):
        """
        Phase signal using suspect's (hidden) phasing functions. You can choose between 3 different types of phasing methods. See suspect/processing/phase.py file.

        * Works only with a 1D [timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        suspect_method : suspect_phasing_method
            Suspect phasing method to use here
        ppm_range : list [2]
            Range in ppm when analyzing spectra for phasing
        allowed_apodization : float/boolean
            If >0 or !=False, apodize signal during phase analysis process. However, the final corrected signal will not be apodized.
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_phased : MRSData2 numpy array [timepoints]
            Resulting phased data stored in a MRSData2 object
        """
        log.debug("phasing using suspect functions [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        s_analyze = self.correct_apodization_nd(allowed_apodization)

        # estimate phases
        if(suspect_method == suspect_phasing_method.MATCH_MAGNITUDE_REAL):
            phi0, phi1 = suspect.processing.phase.mag_real(s_analyze, range_ppm=ppm_range)
        elif(suspect_method == suspect_phasing_method.MIN_IMAG_INTEGRAL):
            phi0, phi1 = suspect.processing.phase.ernst(s_analyze)
        elif(suspect_method == suspect_phasing_method.ACME):
            phi0, phi1 = suspect.processing.phase.acme(s_analyze, range_ppm=ppm_range)
        else:
            log.error("hey, I do not know this suspect phasing method!?")

        # apply phase corrections
        s_phased = s.adjust_phase(phi0, phi1)

        # convert back to MRSData2
        s_phased = self.inherit(s_phased)

        if(display):
            fig_title = "Phasing [%s] using suspect's routines" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 1, sharex='all', sharey='all')

            axs[0].plot(s.frequency_axis_ppm(), s.spectrum().real, 'k-', linewidth=1)
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('original')
            axs[0].grid('on')
            # add low/high cuts
            axs[0].axvline(x=ppm_range[0], color='r', linestyle='--')
            axs[0].axvline(x=ppm_range[1], color='r', linestyle='--')

            axs[1].plot(s_phased.frequency_axis_ppm(), s_phased.spectrum().real, 'b-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('phased')
            axs[1].grid('on')
            # add low/high cuts
            axs[1].axvline(x=ppm_range[0], color='r', linestyle='--')
            axs[1].axvline(x=ppm_range[1], color='r', linestyle='--')

            fig.subplots_adjust()
            fig.show()

        # if any ref data available, we phase it too (silently)
        if(s_phased.data_ref is not None):
            s_phased.data_ref = s_phased.data_ref.correct_phase_1d(suspect_method, ppm_range)

        return(s_phased)

    def correct_first_order_phase_1d(self, coeff_rad_ppm=0.15, display=False, display_range=[1, 6]):
        """
        Correct for first-order phase along the chemical shift axis using a coefficient set manually.

        * Works only with a 1D [timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        coeff_rad_ppm : float
            First-order phase coefficient in rad/ppm
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_phased : MRSData2 numpy array [timepoints]
            Resulting phased data stored in a MRSData2 object
        """
        log.debug("phasing using first-order phasing [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()

        # get spectrum and apply linear phase along ppm axis
        sf_phased = s.spectrum() * np.exp(1j * np.pi * (s.frequency_axis_ppm() - s.ppm0) * coeff_rad_ppm)

        # fft back and convert back to MRSData2
        s_phased = np.fft.ifft(np.fft.ifftshift(sf_phased))

        # convert back to MRSData2
        s_phased = self.inherit(s_phased)

        if(display):
            fig_title = "First-order Phasing [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='all', sharey='col')

            axs[0, 0].plot(s.frequency_axis_ppm(), s.spectrum().real, 'k-', linewidth=1)
            axs[0, 0].set_xlim(display_range[1], display_range[0])
            axs[0, 0].set_xlabel('chemical shift (ppm)')
            axs[0, 0].set_ylabel('original real spectrum')
            axs[0, 0].grid('on')

            axs[0, 1].plot(s.frequency_axis_ppm(), np.unwrap(np.angle(s.spectrum())), 'k-', linewidth=1)
            axs[0, 1].set_xlim(display_range[1], display_range[0])
            axs[0, 1].set_xlabel('chemical shift (ppm)')
            axs[0, 1].set_ylabel('original unwrapped phase (rad)')
            axs[0, 1].grid('on')

            axs[1, 0].plot(s_phased.frequency_axis_ppm(), s_phased.spectrum().real, 'b-', linewidth=1)
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('phased real spectrum')
            axs[1, 0].grid('on')

            axs[1, 1].plot(s_phased.frequency_axis_ppm(), np.unwrap(np.angle(s_phased.spectrum())), 'b-', linewidth=1)
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].set_xlabel('chemical shift (ppm)')
            axs[1, 1].set_ylabel('phased real spectrum')
            axs[1, 1].grid('on')

            fig.subplots_adjust()
            fig.show()

        # if any ref data available, we phase it too (silently)
        if(s_phased.data_ref is not None):
            s_phased.data_ref = s_phased.data_ref.correct_first_order_phase_1d(coeff_rad_ppm)

        return(s_phased)

    def correct_apodization_nd(self, apo_factor=1.0, display=False, display_range=[1, 6]):
        """
        Apodize signal using an exponential window adjusted by a linewidth parameter in Hz.

        * Works with multi-dimensional signals.

        Parameters
        ----------
        apo_factor : float
            Apodization factor in Hz
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_apo : MRSData2 numpy array [whatever,...,timepoints]
            Resulting apodized data stored in a MRSData2 object
        """
        log.debug("apodizing [%s]..." % self.display_label)

        # apodize each individual signal
        s = self.copy()

        # check apodization factor
        if(apo_factor <= 0):
            log.warning("apodization factor is zero or negative, skipping!")
            return(s)

        t = s.time_axis()
        w_apo = np.exp(-apo_factor * t)
        if(s.ndim == 1):
            w_apo_nd = w_apo
        else:  # >1
            w_apo_nd = np.tile(w_apo, list(s.shape[:-1]) + [1])
        s_apo = s * w_apo_nd

        if(display):
            # reshaping
            if(s.ndim == 3):
                s_disp = s.reshape([s.shape[0] * s.shape[1], s.shape[2]])
                s_apo_disp = s_apo.reshape([s_apo.shape[0] * s_apo.shape[1], s_apo.shape[2]])
            else:
                s_disp = s.copy()
                s_apo_disp = s_apo.copy()

            ppm = s_disp.frequency_axis_ppm()
            ppm_apo = s_apo_disp.frequency_axis_ppm()

            fig_title = "Apodizing [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='row', sharey='row')

            axs[0, 0].plot(t, np.abs(s_disp).transpose(), 'k-', linewidth=1, label='fid')
            axs[0, 0].plot(t, w_apo * np.abs(s_disp.max()), 'r-', linewidth=1, label='apodization window')
            axs[0, 0].set_xlabel('time (s)')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(t, np.abs(s_apo_disp).transpose(), 'b-', linewidth=1)
            axs[0, 1].set_xlabel('time (s)')
            axs[0, 1].set_ylabel('apodized')
            axs[0, 1].grid('on')

            axs[1, 0].plot(ppm, s_disp.spectrum().real.transpose(), 'k-', linewidth=1)
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('original spectrum')
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].grid('on')

            axs[1, 1].plot(ppm_apo, s_apo_disp.spectrum().real.transpose(), 'b-', linewidth=1)
            axs[1, 1].set_xlabel("chemical shift (ppm)")
            axs[1, 1].set_ylabel('apodized spectrum')
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].grid('on')

            fig.subplots_adjust()
            fig.show()

        # convert back to MRSData2
        s_apo = self.inherit(s_apo)

        # if any ref data available, we apodize it too (silently)
        if(s_apo.data_ref is not None):
            s_apo.data_ref = s_apo.data_ref.correct_apodization_nd(apo_factor, False)

        return(s_apo)

    def correct_crop_1d(self, nPoints_final=6144, display=False, display_range=[1, 6]):
        """
        Crop signal in time-domain to remove last points.

        * Works only with a 1D [timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        nPoints_final : int
            Final number of points (after crop)
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_crop : MRSData2 numpy array [timepoints]
            Resulting cropped data stored in a MRSData2 object
        """
        log.debug("cropping [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()

        # crop
        if(nPoints_final < s.shape[0]):
            log.debug("cropping data from %d to %d points..." % (s.shape[0], nPoints_final))
            s_crop = s[0:nPoints_final]
        else:
            s_crop = self.copy()
            log.debug("no cropping needed, getting bored...")

        if(display):
            t = s.time_axis()
            t_crop = s_crop.time_axis()
            ppm = s.frequency_axis_ppm()
            ppm_crop = s_crop.frequency_axis_ppm()

            fig_title = "Cropping [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2, sharex='row', sharey='row')

            axs[0, 0].plot(t, np.abs(s), 'k-', linewidth=1, label='fid')
            axs[0, 0].set_xlabel('time (s)')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(t_crop, np.abs(s_crop), 'b-', linewidth=1)
            axs[0, 1].set_xlabel('time (s)')
            axs[0, 1].set_ylabel('cropped')
            axs[0, 1].grid('on')

            axs[1, 0].plot(ppm, s.spectrum().real, 'k-', linewidth=1)
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('original spectrum')
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].grid('on')

            axs[1, 1].plot(ppm_crop, s_crop.spectrum().real, 'b-', linewidth=1)
            axs[1, 1].set_xlabel("chemical shift (ppm)")
            axs[1, 1].set_ylabel('cropped spectrum')
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].grid('on')

            fig.subplots_adjust()
            fig.show()

        # convert back to MRSData2
        s_crop = self.inherit(s_crop)

        # now we have a MRSData2 obj, modify sequence attribute
        if(s_crop.sequence is not None):
            log.debug("updating sequence.npts...")
            s_crop.sequence.npts = nPoints_final
            s_crop.sequence._ready = False

        # if any ref data available, we crop it too (silently)
        if(s_crop.data_ref is not None):
            s_crop.data_ref = s_crop.data_ref.correct_crop_1d(nPoints_final, False)

        return(s_crop)

    def correct_peak_removal_1d(self, hsvd_nComponents=5, hsvd_range=[4.6, 4.8], display=False, display_range=[1, 6]):
        """
        Remove any peak(s) within a ppm range using HSVD. Usually used to remove residual water peak.

        * Works only with a 1D [timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        hsvd_nComponents : int
            Number of components for HSVD
        hsvd_range : list [2]
            Range in ppm of HSVD components
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_peak_removed : MRSData2 numpy array [timepoints]
            Resulting water HSVD suppressed data stored in a MRSData2 object
        """
        log.debug("removing peak for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()
        pbar = log.progressbar("removing peak(s) with HSVD", 5)

        # estimate HSVD components
        components = suspect.processing.water_suppression.hsvd(s, hsvd_nComponents)
        pbar.update(1)

        # filter them by keeping the ones contributing to the residual water peak and its sidebands
        water_components = [component for component in components if ((4.7 - component["frequency"] / self.f0) > hsvd_range[0] and (4.7 - component["frequency"] / self.f0) < hsvd_range[1])]
        pbar.update(2)

        # reconstruct the estimated water peak
        hsvd_fid = suspect.processing.water_suppression.construct_fid(water_components, s.time_axis())
        pbar.update(3)

        # rebuild object
        hsvd_fid = s.inherit(hsvd_fid)
        pbar.update(4)

        # and substract it from the fid
        s_peak_removed = s - hsvd_fid
        pbar.update(5)

        # display this over the data
        if(display):
            fig_title = "Removing some peak(s) [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 1, sharex='all', sharey='all')

            # original spectrum
            axs[0].plot(ppm, s.spectrum().real, 'k-', linewidth=1, label='original data (real part)')
            axs[0].plot(ppm, hsvd_fid.spectrum().real, 'r-', linewidth=1, label='estimated HSVD peak')
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('original spectrum')
            axs[0].grid('on')
            axs[0].legend()

            # water removed spectrum
            axs[1].plot(ppm, s_peak_removed.spectrum().real, 'b-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('peak removed spectrum')
            axs[1].grid('on')
            axs[1].legend()

            fig.subplots_adjust()
            fig.show()

        pbar.finish("done")

        # convert back to MRSData2
        s_peak_removed = self.inherit(s_peak_removed)

        return(s_peak_removed)

    def correct_freqshift_1d(self, peak_range=[4.5, 5], peak_real_ppm=4.7, allowed_apodization=1.0, display=False, display_range=[1, 6]):
        """
        Shift the spectrum in frequency in order to get the right peaks at the right chemical shifts.

        * Works only with a 1D [timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to find a peak of interest
        peak_real_ppm : float
            Chemical shift to set to the peak found
        allowed_apodization : float/boolean
            If >0 or !=False, apodize signal during peak analysis process. However, the final corrected signal will not be apodized.
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_shifted : MRSData2 numpy array [timepoints]
            Resulting frequency calibrated data stored in a MRSData2 object
        """
        log.debug("calibrating [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()

        # find maximum peak in range and its chemical shift
        ppm_peak, peak_val, _, _, _ = s._analyze_peak_1d(peak_range, allowed_apodization)
        log.debug("peak detected at %0.2fppm -> %0.2fppm!" % (ppm_peak, peak_real_ppm))

        # estimate frequency shift in Hz
        log.debug("frequency shifting data...")
        dppm = (peak_real_ppm - ppm_peak)
        df = dppm * s.f0
        s_shifted = s.adjust_frequency(-df)

        if(display):
            fig_title = "Calibrating/frequency shifting [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 1, sharex='all', sharey='all')

            axs[0].plot(ppm, s.spectrum().real, 'k-', linewidth=1)
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('original')
            axs[0].grid('on')
            # add peak position
            axs[0].plot(ppm_peak, np.real(peak_val), 'ro')
            axs[0].axvline(x=ppm_peak, color='r', linestyle='--')

            axs[1].plot(ppm, s_shifted.spectrum().real, 'b-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('shifted')
            axs[1].grid('on')
            # add new peak position
            axs[1].plot(peak_real_ppm, np.real(peak_val), 'ro')
            axs[1].axvline(x=peak_real_ppm, color='r', linestyle='--')

            fig.subplots_adjust()
            fig.show()

        # convert back to MRSData2
        s_shifted = self.inherit(s_shifted)

        return(s_shifted)

    def correct_bandpass_filtering_1d(self, range_ppm=[0, 6], window_func=np.hanning, display=False, display_range=[1, 6]):
        """
        Filter the signal using FFT windowing.

        * Works only with a 1D [timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        range_ppm : list [2]
            Range in ppm used for band-pass filtering
        window_func : numpy windowing function
            Apodization window function
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_filtered : MRSData2 numpy array [timepoints]
            Resulting frequency filtered signal
        """
        log.debug("band-pass filtering [%s]: keeping the %.2f-%.2fppm region..." % (self.display_label, range_ppm[0], range_ppm[1]))
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()

        # build apodization window
        ind_low = np.argmin(np.abs(ppm - range_ppm[0]))
        ind_high = np.argmin(np.abs(ppm - range_ppm[1]))
        n_window = ind_low - ind_high
        window_segment = window_func(n_window)
        window_full = ppm * 0.0
        window_full[ind_high:ind_low] = window_segment

        # apply window
        sf = s.spectrum()
        sf_filtered = sf * window_full
        s_filtered = np.fft.ifft(np.fft.ifftshift(sf_filtered))

        # convert back to MRSData2
        s_filtered = self.inherit(s_filtered)

        if(display):
            fig_title = "Band-pass filtering [%s]" % self.display_label
            fig = plt.figure(fig_title)
            fig.clf()
            fig.suptitle(fig_title)
            axs = fig.subplots(2, 2)

            axs[0, 0].plot(s.time_axis(), np.real(s), 'k-', linewidth=1)
            axs[0, 0].set_xlabel('time (s)')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(ppm, s.spectrum().real, 'k-', linewidth=1)
            axs[0, 1].plot(ppm, window_full * np.max(s.spectrum().real), 'r-', linewidth=1)
            axs[0, 1].set_xlim(display_range[1], display_range[0])
            axs[0, 1].set_xlabel('chemical shift (ppm)')
            axs[0, 1].set_ylabel('original')
            axs[0, 1].grid('on')
            # add low/high cuts
            axs[0, 1].axvline(x=range_ppm[0], color='r', linestyle='--')
            axs[0, 1].axvline(x=range_ppm[1], color='r', linestyle='--')

            axs[1, 0].plot(s.time_axis(), np.real(s_filtered), 'k-', linewidth=1)
            axs[1, 0].set_xlabel('time (s)')
            axs[1, 0].set_ylabel('filtered')
            axs[1, 0].grid('on')

            axs[1, 1].plot(ppm, s_filtered.spectrum().real, 'b-', linewidth=1)
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].set_xlabel('chemical shift (ppm)')
            axs[1, 1].set_ylabel('filtered')
            axs[1, 1].grid('on')
            # add low/high cuts
            axs[1, 1].axvline(x=range_ppm[0], color='r', linestyle='--')
            axs[1, 1].axvline(x=range_ppm[1], color='r', linestyle='--')

            fig.subplots_adjust()
            fig.show()

        return(s_filtered)

    def display_voi_anatomy_nd(self):
        """
        Display the VOI on a anatomical image of your choice. Experimental for now, just following tutorial here: https://suspect.readthedocs.io/en/latest/notebooks/tut06_mpl.html

        * Works with multi-dimensional signals.
        """
        log.debug("displaying VOI [%s]..." % self.display_label)

        # read dicom images
        t2w = suspect.image.load_dicom_volume(self.anatomy_folderpath + "/original-primary-m-norm-nd_e01_0001.dcm")

        # find best slice
        pcg_centre = self.to_scanner(0, 0, 0)
        pcg_centre_index = t2w.from_scanner(*pcg_centre).round().astype(int)

        # VOI drawing coordinates
        # TODO: this works only for sagittal images (spinal cord)
        vx = self.sequence.voxel_size  # mm (will need cm)
        corner_coords_pcg = [[0,    -vx[1] / 20.0,   -vx[2] / 20.0],
                             [0,    -vx[1] / 20.0,   vx[2] / 20.0],
                             [0,    vx[1] / 20.0,    vx[2] / 20.0],
                             [0,    vx[1] / 20.0,    -vx[2] / 20.0],
                             [0,    -vx[1] / 20.0,   -vx[2] / 20.0]]
        corner_coords = np.array([t2w.from_scanner(*self.to_scanner(*coord)) for coord in corner_coords_pcg])


        fig_title = "Anatomical image display [%s]" % self.display_label
        fig = plt.figure(fig_title)
        fig.clf()
        fig.suptitle(fig_title)
        plt.imshow(t2w[pcg_centre_index[2]], cmap=plt.cm.gray)
        plt.plot(corner_coords[:, 0], corner_coords[:, 1], 'red')
        plt.xticks([])
        plt.yticks([])
        plt.show()

    def display_spectrum_1d(self, ifig="Displaying final spectra", display_range=[1, 6], allowed_apodization=5.0, magnitude_mode=False):
        """
        Display spectrum in figure 'ifig', overlaying if needed.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        ifig: int or str
            The figure index that shoud host the plot
        s : MRSData2 numpy array [timepoints]
            MRS data to display
        display_range : list [2]
            Range in ppm used for display
        allowed_apodization : float
            Apodization factor used for display (Hz)
        magnitude_mode : boolean
            Displays in magnitude mode (True) or the real part (False)

        Returns
        -------
        fig : matplotlib.figure
            Resulting matplotlib figure
        """
        log.debug("displaying [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.correct_apodization_nd(allowed_apodization)
        log.debug("displaying stuff!")

        fig = plt.figure(ifig)
        fig.suptitle("Displaying final spectra")
        if(magnitude_mode):
            plt.plot(s.frequency_axis_ppm(), np.abs(s.spectrum()) + self.display_offset, linewidth=1, label=self.display_label)
        else:
            plt.plot(s.frequency_axis_ppm(), s.spectrum().real + self.display_offset, linewidth=1, label=self.display_label)

        # add ytick if offset
        if(self.display_offset != 0):
            yt = plt.yticks()
            yt2 = np.hstack((yt[0], self.display_offset))
            yt3 = np.sort(yt2)
            plt.yticks(yt3)

        if any(display_range):
            plt.xlim(display_range[1], display_range[0])

        plt.xlabel('chemical shift (ppm)')
        plt.ylabel('spectrum')
        plt.grid('on')
        plt.legend()

        plt.subplots_adjust()
        plt.show()

        return(plt.figure(ifig))

    def save_ismrmd(self, h5_filepath):
        """
        Save the MRSData2 object to a ISMRMRD format file.

        Parameters
        ----------
        h5_filepath: string
            Full absolute file path pointing to h5 file
        """
        io.write_ismrmd(self, h5_filepath)

    def save_nifti_mrs(self, nifti_mrs_filepath):
        """
        Save this MRS signal to a NIFTI MRS file.

        Parameters
        ----------
        nifti_mrs_filepath: string
            Full absolute file path pointing to the nifti file
        """
        io.write_nifti_mrs(self, nifti_mrs_filepath)

    def save_mat(self, mat_filepath):
        """
        Save the numpy array content to a MATLAB mat file.

        Parameters
        ----------
        mat_filepath: string
            Full absolute file path pointing to the mat file
        """
        io.write_mat(self, mat_filepath)

    def save_pkl(self, pkl_filepath):
        """
        Save the whole object to a pickle file.

        Parameters
        ----------
        pkl_filepath: string
            Full absolute file path pointing to pkl file
        """
        io.write_pkl(self, pkl_filepath)

    def __reduce__(self):
        """Reduce internal pickling method used when dumping. Modified so that MRSData2 attributes are not forgotten. See for more info: https://docs.python.org/3/library/pickle.html ."""
        # get numpy reduce tuple
        rd = super().__reduce__()
        # add MRSData2 attributes
        rd2 = rd[2] + (self.__dict__,)
        # return the new reduce tuple version
        return(rd[0], rd[1], rd2)

    def __setstate__(self, d):
        """Set new state to object. Internal pickling method used when loading. Modified so that MRSData2 attributes are not forgotten. See for more info: https://docs.python.org/3/library/pickle.html ."""
        # load MRSData2 attributes
        self.__dict__ = d[-1]
        # load all the rest relative to numpy
        super().__setstate__(d[0:-1])
        return(self)

    def to_dataframe(self, include_obj=True, prefix_str="data_"):
        """
        Convert the object's attributes to dataframe. Can include the object itself.

        Parameters
        ----------
        include_obj : boolean
            Include self to the dataframe row
        prefix_str : string
            Prefix string to add to column names

        Returns
        -------
        df : Dataframe
            Containing the attributes as columns (a single row)
        """
        log.debug("converting to dataframe...")

        # get all attributes but remove the private ones
        df = pd.DataFrame.from_dict([vars(self)])
        df = df.filter(regex=("^(?!_).*"))

        # add some specific private attributes I need
        df["display_label"] = self.display_label
        df["noise_level"] = self.noise_level
        df["rejection"] = [self.data_rejection]  # can be a list
        df["file_hash"] = self.data_file_hash
        df["is_rawdata"] = self.is_rawdata
        df = pd.concat([df, pd.DataFrame.from_dict([self.patient])], axis=1)
        df = pd.concat([df, self.sequence.to_dataframe(True)], axis=1)

        if(include_obj):
            df["obj"] = [self]

        df = df.add_prefix(prefix_str)

        return(df)


class pipeline:
    """The pipeline class is used to store all the reconstruction parameters needed for a specific bunch of acquired signals. Once the parameters are set, the pipeline can be run using one of the methods."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, template_name=None):
        """
        Initialize the reconstruction pipeline, using a template if needed.

        Parameters
        ----------
        template_name: string
            Reco pipeline template file
        """
        # --- initializing dataset dict ---
        self.dataset = [{}] * MAX_NUM_DATASETS
        for i in range(MAX_NUM_DATASETS):
            self.dataset[i] = {"legend": None,
                               "raw": {"files": [None, None], "data": None, "analysis_results": None, "ref_data_analysis_results": None},
                               "dcm": {"files": [None, None], "data": None, "analysis_results": None, "ref_data_analysis_results": None},
                               "physio_file": None,
                               "imaging_file": None,
                               "resp_bpm": None,
                               "heart_bpm": None,
                               "comment": None}

        # --- global settings ---
        self.settings = {   # option to process only a set of datasets: list of indexes
                            "datasets_indexes": None,
                            # folder to search for supplementary datasets
                            "folder_additional_datasets": None,
                            # ppm scale reference
                            "ppm0": 4.7,
                            # ppm range to search for peak used for phasing, etc.
                            "POI_range_ppm": [4.5, 5.2],
                            # ppm range to search for ppm scale calibration
                            "POI_shift_range_ppm": [4.5, 5.2],
                            # real ppm value the above peak
                            "POI_shift_true_ppm": 4.7,
                            # ppm range to search for peak for SNR estimation
                            "POI_SNR_range_ppm": [1.8, 2.2],
                            # ppm range to search for peak for FWHM estimation
                            "POI_LW_range_ppm": [4.5, 5.2],
                            # ppm range to for SNR/LW estimation in ref. data
                            "POI_ref_range_ppm": [4.5, 5.2],
                            # apodization factor used during signal analysis, never actually applied for signal correction
                            "allowed_apodization": 1.0,
                            # no phasing using reference data for already reconstructed data
                            "no_phasing_using_ref_for_dcm_data": True,
                            # path to pkl file to store processed data
                            "storage_file": None,
                            # force display off if needed
                            "display": None,
                            # ppm range used for display
                            "display_range_ppm": [1, 6],
                            # y offset used for display
                            "display_offset": 0.0,
                            # automatically check which dataset is noWS or WS and swap if needed
                            "auto_detect_ref_scan": True,
                            # raise error if raw data looks worse than dcm
                            "raise_error_on_badreco": True}

        # --- available jobs and their parameters ---
        self.job = {}
        # --- job: spectrum final display ---
        self.job["displaying"] = {"job_func": MRSData2.display_spectrum_1d, "job_name": "displaying",
                                  # figure index
                                  "fig_index": "Displaying final spectra",
                                  # ppm range used for display
                                  "display_range_ppm": pipeline._get_setting,
                                  # apodization factor used for display (Hz)
                                  "allowed_apodization": pipeline._get_setting,
                                  # display spectrum in magnitude mode?
                                  "magnitude_mode": False
                                  }

        # --- job: VOI display on anatomical images ---
        self.job["displaying_anatomy"] = {"job_func": MRSData2.display_voi_anatomy_nd, "job_name": "overlaying VOI on anatomical image"
                                  }

        # --- job: automatic rephasing ---
        self.job["phasing"] = {"job_func": MRSData2.correct_phase_3d, "job_name": "phasing",
                               # use reference data is available?
                               "using_ref_data": True,
                               # ppm range to look fo peak used to estimate phase
                               "POI_range_ppm": pipeline._get_setting,
                               # average all averages per channel
                               "average_per_channel_mode": False,
                               # measure phase from 1st time point
                               "first_point_fid_mode": False,
                               # order of phasing in time: 0th or 1st order
                               "order": 0,
                               # add an additional 0th order phase (rd)
                               "offset": 0.0,
                               # display all this process to check what the hell is going on
                               "display": False,
                               "display_range_ppm": pipeline._get_setting
                               }

        # --- job: amplification ---
        self.job["scaling"] = {"job_func": MRSData2.correct_intensity_scaling_nd, "job_name": "scaling intensity",
                               "scaling_factor_rawdata": 1e8,
                               "scaling_factor_dcm": 1.0
                               }

        # --- job: FID modulus ---
        self.job["FID modulus"] = {"job_func": MRSData2.correct_fidmodulus_nd, "job_name": "FID modulus"
                                   }

        # --- job: time_shifting ---
        self.job["time_shifting"] = {"job_func": MRSData2.correct_time_shift_nd, "job_name": "time_shifting",
                                     # time shift in us
                                     "time_shift_us": 375,
                                     # display all this process to check what the hell is going on
                                     "display": True,
                                     "display_range_ppm": pipeline._get_setting
                                     }

        # --- job: channel combination ---
        self.job["channel_combining"] = {"job_func": MRSData2.correct_combine_channels_3d, "job_name": "channel-combining",
                                         # use non water-suppressed data to recombine and rephase channels
                                         "using_ref_data": True,
                                         # should we rephase (0th order) data while combining?
                                         "phasing": False,
                                         # boolean mask to switch on/off some Rx channels
                                         "weights": [True]
                                         }

        # --- job: concatenate ---
        self.job["concatenate"] = {"job_func": MRSData2.concatenate_2d, "job_name": "concatenate"
                                   }

        # --- job: zero_filling ---
        self.job["zero_filling"] = {"job_func": MRSData2.correct_zerofill_nd, "job_name": "zero-filling",
                                    # number of signal points after zf
                                    "npts": 8192 * 2,
                                    # display all this process to check what the hell is going on
                                    "display": True,
                                    "display_range_ppm": pipeline._get_setting
                                    }

        # --- job: analyze physio signal ---
        self.job["physio_analysis"] = {"job_func": MRSData2.analyze_physio_2d, "job_name": "analyzing physio. signals",
                                       # ppm range to look for a peak to analyze
                                       "POI_range_ppm": pipeline._get_setting,
                                       # time range in (ms) to look around timestamp for correlation physio/MRS
                                       "delta_time_ms": 1000.0,
                                       # apodization factor used during signal analysis stage
                                       "allowed_apodization": pipeline._get_setting,
                                       # display all this process to check what the hell is going on
                                       "display": True
                                       }

        # --- job: automatic data rejection based on criterias ---
        self.job["data_rejecting"] = {"job_func": MRSData2.correct_analyze_and_reject_2d, "job_name": "data rejecting",
                                      # ppm range to look for a peak to analyze
                                      "POI_range_ppm": pipeline._get_setting,
                                      # ppm range to estimate SNR
                                      "POI_SNR_range_ppm": pipeline._get_setting,
                                      # ppm range to estimate LW
                                      "POI_LW_range_ppm": pipeline._get_setting,
                                      # size of moving average window
                                      "moving_averages": 1,
                                      # if True, rejects if linewidth could not be estimated
                                      "reject_when_linewidth_fails": True,
                                      # rejection criterias for
                                      # amplitude relative changes: keep data if within +/-val % range
                                      # linewidth changes: keep data is below val Hz
                                      # chemical shift changes: keep data is within +/-val ppm
                                      # phase changes: keep data if within +/-val/100 * std(phase) rd
                                      "ranges": {"amplitude (%)": None,
                                                 "linewidth (Hz)": None,
                                                 "chemical shift (ppm)": None,
                                                 "phase std. factor (%)": None},
                                      # for amplitude, chemical shift and phase, the rejection of data is based on ranges of relative changes of those metrics. Relative to what? The man value over the whole acquisition (True) or the first acquired point (False)
                                      "rel2mean": True,
                                      # method for automatic adjustement
                                      "auto_method_list": [data_rejection_method.AUTO_AMPLITUDE,
                                                           data_rejection_method.AUTO_LINEWIDTH,
                                                           data_rejection_method.AUTO_FREQUENCY,
                                                           data_rejection_method.AUTO_PHASE],
                                      # minimum allowed SNR change (%) when adjusting the linewidth criteria, this can be positive (we want to increase SNR +10% by rejecting crappy dat) or negative (we are ok in decreasing the SNR -10% in order to get better resolved spectra)
                                      "auto_allowed_snr_change": 1.0,
                                      # apodization factor used during signal analysis stage
                                      "allowed_apodization": pipeline._get_setting,
                                      # intercorrelation mode for realignment?
                                      "display": True,
                                      "display_range_ppm": pipeline._get_setting
                                      }

        # --- job: automatic data frequency realignment ---
        self.job["realigning"] = {"job_func": MRSData2.correct_realign_2d, "job_name": "frequency realigning",
                                  # ppm range to look for a peak to analyze
                                  "POI_range_ppm": pipeline._get_setting,
                                  # size of moving average window
                                  "moving_averages": 1,
                                  # use correlation mode
                                  "inter_corr_mode": False,
                                  # maximum frequency shift allowed in Hz
                                  "freq_shift_max": 25,
                                  # apodization factor used during signal analysis stage
                                  "allowed_apodization": pipeline._get_setting,
                                  # display all this process to check what the hell is going on
                                  "display": True,
                                  "display_range_ppm": pipeline._get_setting
                                  }

        # --- job: spectral filtering ---
        self.job["filtering"] = {"job_func": MRSData2.correct_bandpass_filtering_1d, "job_name": "FFT filtering",
                                 # ppm range to keep
                                 "range_ppm": [1, 6],
                                 # type of apodization window (take it from numpy/scipy)
                                 "window_func": signal.tukey,
                                 # display all this process to check what the hell is going on
                                 "display": True,
                                 "display_range_ppm": pipeline._get_setting
                                 }

        # --- job: data averaging ---
        self.job["averaging"] = {"job_func": MRSData2.correct_average_2d, "job_name": "averaging",
                                 # number of averages to mean (None = all)
                                 "na": None,
                                 # display all this process to check what the hell is going on
                                 "display": True,
                                 "display_range_ppm": pipeline._get_setting
                                 }

        # --- job: phasing using suspect ---
        self.job["phasing_suspect"] = {  "job_func": MRSData2.correct_phase_1d, "job_name": "phasing (suspect)",
                                         # phasing method
                                         "suspect_method": suspect_phasing_method.MATCH_MAGNITUDE_REAL,
                                         # ppm range to analyze phase
                                         "range_ppm": [1, 6],
                                         # apodization factor used during signal analysis stage
                                         "allowed_apodization": pipeline._get_setting,
                                         # display all this process to check what the hell is going on
                                         "display": True,
                                         "display_range_ppm": pipeline._get_setting
                                         }

        # --- job: first-order phasing ---
        self.job["phasing_first_order"] = {  "job_func": MRSData2.correct_first_order_phase_1d, "job_name": "first-order phasing",
                                             # phasing coeficient
                                             "coeff_rad_ppm": 0.15,
                                             # display all this process to check what the hell is going on
                                             "display": True,
                                             "display_range_ppm": pipeline._get_setting
                                             }

        # --- job: noise level analysis ---
        self.job["noise_estimation"] = {"job_func": MRSData2.analyze_noise_nd, "job_name": "estimating noise level",
                                        # estimate noise std time-domain on the last 100 pts of the FID
                                        "npts": 100
                                        }

        # --- job: data apodization ---
        self.job["apodizing"] = {"job_func": MRSData2.correct_apodization_nd, "job_name": "apodizing",
                                 # exponential damping factor for apodization (Hz)
                                 "damping_hz": 5,
                                 # display all this process to check what the hell is going on
                                 "display": True,
                                 "display_range_ppm": pipeline._get_setting
                                 }

        # --- job: data cropping ---
        self.job["cropping"] = {"job_func": MRSData2.correct_crop_1d, "job_name": "cropping",
                                # final number of signal points after crop
                                "final_npts": 6144,
                                # display all this process to check what the hell is going on
                                "display": True,
                                "display_range_ppm": pipeline._get_setting
                                }

        # --- job: water post-acquisition removal ---
        self.job["water_removal"] = {"job_func": MRSData2.correct_peak_removal_1d, "job_name": "removing water peak",
                                     # number of components when running HSVD
                                     "hsvd_components": 5,
                                     # ppm range where all components will be remove
                                     "POI_range_ppm": pipeline._get_setting,
                                     # display all this process to check what the hell is going on
                                     "display": True,
                                     "display_range_ppm": pipeline._get_setting
                                     }

        # --- job: spectrum chemical shift calibration ---
        self.job["calibrating"] = {"job_func": MRSData2.correct_freqshift_1d, "job_name": "frequency shifting",
                                   # ppm range to look for the peak of interest (NAA by default)
                                   "POI_shift_range_ppm": pipeline._get_setting,
                                   # real ppm value for this peak
                                   "POI_shift_true_ppm": pipeline._get_setting,
                                   # apodization factor used during signal analysis stage
                                   "allowed_apodization": pipeline._get_setting,
                                   # display all this process to check what the hell is going on
                                   "display": True,
                                   "display_range_ppm": pipeline._get_setting
                                   }

        # --- job: SNR analysis ---
        self.job["analyzing_snr"] = {"job_func": MRSData2.analyze_snr_1d, "job_name": "analyzing SNR",
                                     # ppm range to look for a peak to analyze
                                     "POI_SNR_range_ppm": pipeline._get_setting,
                                     # ppm range to look for pure noise
                                     "n_range_ppm": [-2, -1],
                                     # divide SNR by 2, like LCModel do
                                     "half_factor": False,
                                     # should we look at the magnitude or real spectrum?
                                     "magnitude_mode": False,
                                     # display all this process to check what the hell is going on
                                     "display": True,
                                     "display_range_ppm": pipeline._get_setting
                                     }

        # --- job: LW analysis ---
        self.job["analyzing_lw"] = {"job_func": MRSData2.analyze_linewidth_1d, "job_name": "analyzing peak-linewidth",
                                    # ppm range to look for a peak to analyze
                                    "POI_LW_range_ppm": pipeline._get_setting,
                                    # should we look at the magnitude or real spectrum?
                                    "magnitude_mode": False,
                                    # display all this process to check what the hell is going on
                                    "display": True,
                                    "display_range_ppm": pipeline._get_setting
                                    }

        # --- job: ref data SNR analysis ---
        self.job["ref_data_analyzing_snr"] = {
                                    "job_func": MRSData2.analyze_snr_1d, "job_name": "analyzing ref. data SNR",
                                     # ppm range to look for a peak to analyze
                                     "POI_ref_range_ppm": pipeline._get_setting,
                                     # ppm range to look for pure noise
                                     "n_range_ppm": [-2, -1],
                                     # divide SNR by 2, like LCModel do
                                     "half_factor": False,
                                     # should we look at the magnitude or real spectrum?
                                     "magnitude_mode": False,
                                     # display all this process to check what the hell is going on
                                     "display": False,
                                     "display_range_ppm": pipeline._get_setting
                                     }

        # --- job: ref data LW analysis ---
        self.job["ref_data_analyzing_lw"] = {
                                    "job_func": MRSData2.analyze_linewidth_1d, "job_name": "analyzing ref. data peak-linewidth",
                                    # ppm range to look for a peak to analyze
                                    "POI_ref_range_ppm": pipeline._get_setting,
                                    # should we look at the magnitude or real spectrum?
                                    "magnitude_mode": False,
                                    # display all this process to check what the hell is going on
                                    "display": False,
                                    "display_range_ppm": pipeline._get_setting
                                    }

        # --- job list ---
        # list of data processing to apply to the data
        # beware, you need to know what you are doing here
        # also, be careful with the dimensionality of data 3D, 2D, 1D along the data processing
        # order is important!
        self.job_list = [self.job["phasing"],
                         self.job["scaling"],
                         self.job["FID modulus"],
                         self.job["channel_combining"],
                         self.job["concatenate"],
                         self.job["zero_filling"],
                         self.job["physio_analysis"],
                         self.job["data_rejecting"],
                         self.job["realigning"],
                         self.job["averaging"],
                         self.job["noise_estimation"],
                         self.job["apodizing"],
                         self.job["cropping"],
                         self.job["water_removal"],
                         self.job["calibrating"],
                         self.job["displaying"]]

        # --- analyze job list ---
        # SNR/LW analysis job list
        self.analyze_job_list = [self.job["channel_combining"],
                                 self.job["averaging"],
                                 self.job["calibrating"]]

        # --- SNR/LW analysis ---
        self.analyze_enable = True

        # --- template loading if needed ---
        self.template_name = template_name
        if(template_name is not None):
            # overwrite everything above with the template
            self.load_template(template_name)

        # freeze the object and prevent the creation of new attributes
        self.__isfrozen = True

    def _get_setting(self, setting_key):
        """
        Return value of setting from settings dict. Sounds like a stupid method but I use it to replace pointers. To make it a little easier for the user, a settings dict contains all the main parameters that I repeatably used during the reco, especially the ppm range to look for a peak of interest (POI). For each job defined above, the default parameter value are linked to this settings dict using this method.

        Parameters
        ----------
        setting_key : dict key from self.settings
            Name of the setting

        Returns
        -------
        self.settings[setting_key] : ?
            Value of the setting
        """
        return(self.settings[setting_key])

    def _run_job(self, job, data, default_args=False):
        """
        Estimate SNR and/or peak linewidth for this dataset. Values are stored. A mini default pipeline is applied before SNR/LW measurements and can be set with self.analyze_job_list.

        Parameters
        ----------
        job : dict entry from self.job
            The job to run on the data
        data : MRSData2 object [whatever,...,timepoints]
            Data to process
        default_args : boolean
            Should we ignore the pipeline job parameters and run with default arguments (True)?

        Returns
        -------
        job_result : ?
            Stuff returned by the job
        """
        # get job name
        log.info_line________________________()
        job_name = job["job_name"]
        log.info("%s on [%s]..." % (job_name, data.display_label))
        # get function
        job_func = job["job_func"]
        if(default_args):
            job_args = [data]
        else:
            # get arguments
            job_args = job.copy()
            del job_args["job_name"]
            del job_args["job_func"]
            job_args = [data] + list(job_args.values())

        # call job on data
        job_result = job_func(*job_args)

        # return
        return(job_result)

    def _analyze(self, data, current_job, already_done_jobs):
        """
        Estimate SNR and/or peak linewidth for this dataset. Values are stored. A mini default pipeline is applied before SNR/LW measurements and can be set with self.analyze_job_list.

        Parameters
        ----------
        data : MRSData2 object [whatever,...,timepoints]
            Dataset
        current_job : MRSData2 method function
            The job that was just applied to data before calling this function
        already_done_jobs : list (stack)
            List of already applied processing functions to this dataset

        Returns
        -------
        data_snr : float
            SNR estimated on data
        data_lw : float
            Peak linewidth estimated on data
        ref_data_snr : float
            SNR estimated on ref. data if any
        ref_data_lw : float
            Peak linewidth estimated on ref. data if any
        """
        log.debug("estimating SNR and LW for [%s]..." % data.display_label)

        # init job list
        this_analyze_job_list = [j for j in self.analyze_job_list + already_done_jobs if (j in self.analyze_job_list) and (j not in already_done_jobs)]

        # run mini-pipeline with default arguments (= no display)
        # and with no log output
        log.pause()
        for j in this_analyze_job_list:
            # run job on this dataset with default arguments
            data = self._run_job(j, data, True)

        # measure snr
        this_job = self.job["analyzing_snr"]
        data_snr, _, _ = self._run_job(this_job, data)
        # measure lw
        this_job = self.job["analyzing_lw"]
        data_lw = self._run_job(this_job, data)

        if(data.data_ref is not None):
            # measure ref data snr
            this_job = self.job["ref_data_analyzing_snr"]
            ref_data_snr, _, _ = self._run_job(this_job, data.data_ref)
            # measure ref data lw
            this_job = self.job["ref_data_analyzing_lw"]
            ref_data_lw = self._run_job(this_job, data.data_ref)
        else:
            ref_data_snr = np.nan
            ref_data_lw = np.nan

        # output
        log.resume()
        job_label = "post-" + current_job["job_name"]
        log.debug(job_label + " SNR of [%s] = %.2f" % (data.display_label, data_snr))
        log.debug(job_label + " LW of [%s] = %.2f" % (data.display_label, data_lw))

        return(data_snr, data_lw, ref_data_snr, ref_data_lw)

    def _detect_data_ref(self, s, s_ref):
        """
        Detect highest SNR between these two datasets and return them in the right order.

        Parameters
        ----------
        s : MRSData2 object
            Data supposed to be WS
        s_ref : MRSData2 object
            Data supposed to be noWS and used as a reference signal for phasing, etc.

        Returns
        -------
        s2 : MRSData2 object
            Data supposed to be WS
        s_ref2 : MRSData2 object
            Data supposed to be noWS and used as a reference signal for phasing, etc.
        """
        # compare FID max in magnitude (dirty but robust)
        s_ref_tmp = np.mean(np.abs(s_ref), axis=(0, 1))
        s_tmp = np.mean(np.abs(s), axis=(0, 1))

        if(np.max(s_tmp) > np.max(s_ref_tmp)):
            # need to swap
            log.debug("swapping signal of interest (WS) and reference signal (noWS)!")
            return(s_ref, s)
        else:
            return(s, s_ref)

    def _check_datasets(self, shrink_list=True):
        """
        Check if the datasets list is well initialized by user.

        Parameters
        ----------
        shrink_list : boolean
            After checking, keep only non-empty data entries in the list.
        """
        log.debug("checking datasets...")

        # for each data set, check dict fields
        ind_dataset_ok = []
        for i, d in enumerate(self.dataset):
            legend_ok = (d["legend"] is not None)
            rawdata_ok = (d["raw"]["data"] is not None)
            rawfiles_ok = (d["raw"]["files"][0] is not None)
            dicomfiles_ok = (d["dcm"]["files"][0] is not None)
            physiofile_ok = (d["physio_file"] is not None)
            imagingfile_ok = (d["imaging_file"] is not None)

            # checking that this is an entry
            if(legend_ok or rawfiles_ok or dicomfiles_ok or physiofile_ok or imagingfile_ok):
                # seems like a entry

                # checking legend
                if(not legend_ok):
                    log.error("dataset[%d] is missing a legend! :(" % i)
                # checking datafiles
                if(not rawdata_ok and not rawfiles_ok and not dicomfiles_ok):
                    log.error("dataset[%d] is missing data, no raw or dicom files were set! :(" % i)
                # checking number of datafiles
                if(rawfiles_ok and len(d["raw"]["files"]) > 2):
                    log.error("dataset[%d] has %d raw-files! There should be 1 or 2. :(" % (i, len(d["raw"]["files"])))
                if(dicomfiles_ok and len(d["dcm"]["files"]) > 2):
                    log.error("dataset[%d] has %d dicom-files! There should be 1 or 2. :(" % (i, len(d["dcm"]["files"])))

                # clean and print strings
                log.debug("dataset[%d]:" % i)
                # legend
                if(legend_ok):
                    self.dataset[i]["legend"] = d["legend"].strip()
                log.debug("  legend = " + self.dataset[i]["legend"])
                # raw data
                if(self.dataset[i]["raw"]["files"][0] is not None):
                    self.dataset[i]["raw"]["files"][0] = self.dataset[i]["raw"]["files"][0].strip()
                    log.debug("  raw #0 = " + self.dataset[i]["raw"]["files"][0])
                if(len(self.dataset[i]["raw"]["files"]) == 1):
                    self.dataset[i]["raw"]["files"].append(None)
                if(self.dataset[i]["raw"]["files"][1] is not None):
                    self.dataset[i]["raw"]["files"][1] = self.dataset[i]["raw"]["files"][1].strip()
                    log.debug("  raw #1 [REF] = " + self.dataset[i]["raw"]["files"][1])
                else:
                    log.debug("  raw #1 [REF] = None")
                # dicom
                if(self.dataset[i]["dcm"]["files"][0] is not None):
                    self.dataset[i]["dcm"]["files"][0] = self.dataset[i]["dcm"]["files"][0].strip()
                    log.debug("  dicom #0 = " + self.dataset[i]["dcm"]["files"][0])
                if(len(self.dataset[i]["dcm"]["files"]) == 1):
                    self.dataset[i]["dcm"]["files"].append(None)
                if(self.dataset[i]["dcm"]["files"][1] is not None):
                    self.dataset[i]["dcm"]["files"][1] = self.dataset[i]["dcm"]["files"][1].strip()
                    log.debug("  dicom #1 [REF] = " + self.dataset[i]["dcm"]["files"][1])
                else:
                    log.debug("  dicom #1 [REF] = None")
                # physio
                if(physiofile_ok):
                    self.dataset[i]["physio_file"] = d["physio_file"].strip()
                    log.debug("  physio = " + self.dataset[i]["physio_file"])
                # imaging
                if(imagingfile_ok):
                    self.dataset[i]["imaging_file"] = d["imaging_file"].strip()
                    log.debug("  imaging = " + self.dataset[i]["imaging_file"])

                # index of non-empty datasets
                ind_dataset_ok.append(i)

        # keep only the non-empty datasets
        if(shrink_list):
            self.dataset = [self.dataset[i] for i in ind_dataset_ok]

    def _complete_missing_datasets(self, accepted_dcm_file_ext_list=[".dcm", ".ima"], accepted_raw_file_ext_list=[".dat"]):
        """
        Attempt to complete datasets by finding corresponding reconstructed (DCM/IMA) files or raw data (DAT) files. This method really works with SIEMENS data.

        Parameters
        ----------
        accepted_dcm_file_ext_list : list
            List of file extension to consider for reconstructed data (dicoms generally for now)
        accepted_raw_file_ext_list : list
            List of file extension to consider for raw data (TWIX .dat files for now)
        """
        log.debug("find missing datasets...")

        # now browse folder and subfolders
        if(self.settings["folder_additional_datasets"] is None):
            log.warning("no [folder_additional_datasets] setting specified: will not look for additional datasets!")
            return

        # first let's check if we have missing files
        n_missing_dcm_files = 0
        n_missing_raw_files = 0
        for i, d in enumerate(self.dataset):
            for dc in d["dcm"]["files"]:
                if(dc is None):
                    n_missing_dcm_files += 1
            for dr in d["raw"]["files"]:
                if(dr is None):
                    n_missing_raw_files += 1

        if((n_missing_dcm_files == 0) and (n_missing_raw_files == 0)):
            log.error("no files in dataset list!")
        elif(n_missing_dcm_files != n_missing_raw_files):
            if(n_missing_dcm_files > n_missing_raw_files):
                log.debug("%d reconstructed data (DCM/IMA) files missing!" % (n_missing_dcm_files - n_missing_raw_files))
            else:
                log.debug("%d reconstructed data (DAT/TWIX) files missing!" % (n_missing_raw_files - n_missing_dcm_files))
        else:
            log.debug("no files missing :)")
            return

        # now load previous dataframe if any
        csv_filename = os.path.join(self.settings["folder_additional_datasets"], ".pastis.csv")
        if(os.path.isfile(csv_filename)):
            df = pd.read_csv(csv_filename)
            log.debug("found cache file [%s] containing %d file IDs :)" % (csv_filename, len(df)))
        else:
            log.info("no cache file found: please be patient will browsing folders (next time will be faster)...")
            df = pd.DataFrame({}, columns=['filename', 'dataType', 'lProtID'])

        # convert all relative paths to absolute
        df["filename"] = [os.path.normpath(os.path.join(self.settings["folder_additional_datasets"], f)) for f in df["filename"]]

        # initialize lists
        file_name_list = []
        file_type_list = []
        file_lprotid_list = []
        n_files_read = 0
        for root, dirs, files in os.walk(self.settings["folder_additional_datasets"], topdown=False):

            # read file only if .dat or if .dcm/.ima and if only .dcm/.ima in folder
            this_file_list = []
            n_dcm_files_here = 0
            n_raw_files_here = 0
            for this_filename in files:
                _, this_ext = os.path.splitext(this_filename.lower())
                if(this_ext in accepted_dcm_file_ext_list):
                    n_dcm_files_here += 1
                    this_file_list.append(this_filename)
                elif(this_ext in accepted_raw_file_ext_list):
                    n_raw_files_here += 1
                    this_file_list.append(this_filename)
                else:
                    # weird extension
                    pass

            # check for weird exception: a folder containing a IMA and a DICOM file representing the same data...
            if((n_dcm_files_here == 2) and (n_raw_files_here == 0)):
                _, this_ext0 = os.path.splitext(this_file_list[0].lower())
                _, this_ext1 = os.path.splitext(this_file_list[1].lower())
                if(this_ext0 != this_ext1):
                    # keep only first file
                    this_file_list = [this_file_list[0]]
                    # and pretent there is only one file in there
                    n_dcm_files_here = 1

            if( ((n_dcm_files_here == 0) and (n_raw_files_here > 0)) or (n_dcm_files_here == 1) ):
                # we reach here if (no dicom files but at least on .dat file) or (1 dicom file) is in folder
                # to accelerate search!

                for this_filename in this_file_list:
                    # try to read file if not already in df (accelerate!)
                    this_filename_fullpath = os.path.join(root, this_filename)
                    if(this_filename_fullpath not in df['filename'].to_list()):
                        log.debug("reading header from %s" % this_filename_fullpath)
                        log.pause()
                        try:
                            mfr = io.get_data_file_reader(this_filename_fullpath)
                            pid = mfr.read_param_num("lProtID")
                            file_name_list.append(this_filename_fullpath)
                            file_lprotid_list.append(pid)

                            # file type
                            _, this_ext = os.path.splitext(this_filename.lower())
                            if(this_ext in accepted_dcm_file_ext_list):
                                file_type_list.append('dcm')
                            elif(this_ext in accepted_raw_file_ext_list):
                                file_type_list.append('raw')
                            else:
                                # weird extension
                                pass

                            n_files_read += 1
                        except:
                            pass
                        log.resume()

        log.debug("found %d new files!" % n_files_read)

        # build dataframe and append
        new_df = pd.DataFrame(list(zip(file_name_list, file_type_list, file_lprotid_list)), columns=['filename', 'dataType', 'lProtID'])
        df = pd.concat([df, new_df])

        # now browse dataset list and fix missing files
        for i, this_dataset in enumerate(self.dataset):
            for a, b in zip(["dcm", 'raw'], ["raw", "dcm"]):
                # let's fix only when both files are missing
                if((this_dataset[a]["files"] == [None, None]) and (this_dataset[b]["files"] != [None, None])):
                    # originaly, this code was made to find dicoms knowing twix files: a/b are there to do both
                    # get corresponding raw data files
                    raw_file_list = this_dataset[b]["files"]
                    # find corresponding PIDs in df
                    raw_file_pid_list = df.loc[(df["filename"].isin(raw_file_list) ) &
                                           (df["dataType"] == b)]["lProtID"].tolist()
                    # find corresponding dicom files in df
                    dcm_file_list = df.loc[(df["lProtID"].isin(raw_file_pid_list) ) &
                                           (df["dataType"] == a)]["filename"].tolist()
                    # complete with Nones if needed
                    dcm_file_list = dcm_file_list + [None] * (2 - len(dcm_file_list))
                    # and rewrite dicom file list
                    self.dataset[i][a]["files"] = dcm_file_list
                    log.info("added following files to dataset list: " + str(dcm_file_list))

        # convert all paths to relative
        df["filename"] = [os.path.relpath(f, self.settings["folder_additional_datasets"]) for f in df["filename"]]

        # store dataframe in a .hidden file if possible
        try:
            df.to_csv(csv_filename, index=False)
            log.debug("updated cache file [%s] with a total of %d file IDs :)" % (csv_filename, len(df)))
        except:
            log.debug("could not store cache file in folder [%s]. Probably no write permissions?" % self.settings["folder_additional_datasets"])

        # done: dataset list completed

    def get_te_list(self):
        """
        Return the TEs for all the data signals in this pipeline.

        Returns
        -------
        [s.te for s in self._data_list] : list
            TEs for all signals stored in here
        """
        return([d["raw"]["data"].te for d in self.dataset])

    def run(self):
        """
        Run the pipeline for data MRS reconstruction.

        Returns
        -------
        self._data_list : list of MRSData2 objects
            Final MRS data signals obtained from reconstruction pipeline stored in a MRSData2 object
        """
        # --- reading and checking dataset list ---
        self._check_datasets()

        # try to fix missing files in dataset list
        self._complete_missing_datasets()

        # filter datasets
        if(self.settings["datasets_indexes"] is not None):
            # if only one index, convert to list
            if(type(self.settings["datasets_indexes"]) == int):
                self.settings["datasets_indexes"] = [self.settings["datasets_indexes"]]

            # keep only data to process
            self.dataset = [self.dataset[i] for i in self.settings["datasets_indexes"]]

        # --- reading data ---
        # now let's read the data files
        log.info("reading data files...")
        for i, d in enumerate(self.dataset):
            this_legend = d["legend"]
            this_physio_filename = d["physio_file"]
            this_imaging_filename = d["imaging_file"]

            for dtype in ["raw", "dcm"]:
                # check if any data file to read
                this_data_filename = d[dtype]["files"][0]
                if(this_data_filename is None):
                    # no? so pass
                    continue

                # reading WS data
                log.info("reading data [" + this_legend + "]")
                s = MRSData2(this_data_filename, this_physio_filename, this_imaging_filename)

                # reading noWS data
                this_data_ref_filename = d[dtype]["files"][1]
                if(this_data_ref_filename is not None):
                    log.info_line_break()
                    log.info("reading data [" + this_legend + "]")
                    s_ref = MRSData2(this_data_ref_filename, None)

                    if(self.settings["auto_detect_ref_scan"]):
                        # find ref data (highest snr) and swap if needed
                        s, s_ref = self._detect_data_ref(s, s_ref)
                else:
                    s_ref = None

                log.info("data : got a " + str(s.shape) + " vector")
                # set ppm reference
                s.ppm0 = self.settings["ppm0"]
                # set legend & offset
                this_new_legend = ("#%d " % i) + this_legend
                if(s.is_rawdata):
                    s.set_display_label(this_new_legend + " (RAW)")
                else:
                    s.set_display_label(this_new_legend + " (DCM)")

                s.set_display_offset(self.settings["display_offset"])
                # store
                self.dataset[i][dtype]["data"] = s
                self.dataset[i]["legend"] = this_new_legend

                if(s_ref is not None):
                    log.info_line_break()
                    log.info("ref. data: got a " + str(s_ref.shape) + " vector")
                    # if several averages, mean now (that could be a problem!?)
                    s_ref = s_ref.mean(axis=0)
                    # add 1 dimension
                    s_ref = s_ref.reshape((1,) + s_ref.shape)
                    log.debug("reshaped to a " + str(s_ref.shape) + " vector")
                    # set ppm reference
                    s_ref.ppm0 = self.settings["ppm0"]
                    # set legend & offset
                    this_new_legend = ("#%d " % i) + this_legend
                    if(s_ref.is_rawdata):
                        s_ref.set_display_label(this_new_legend + " [REF] [RAW]")
                    else:
                        s_ref.set_display_label(this_new_legend + " [REF] [DCM]")
                    s_ref.set_display_offset(self.settings["display_offset"])
                    # store
                    self.dataset[i][dtype]["data"].data_ref = s_ref

        # --- applying global settings ---
        log.info_line________________________()
        log.info("reading your job list...")
        # applying some global settings to jobs
        default_display_value = self._get_setting("display")
        for this_job_name in self.job:
            for this_setting_name, this_setting_value in self.job[this_job_name].items():
                # should we use a default value for this parameter?
                if(this_setting_value == pipeline._get_setting):
                    self.job[this_job_name][this_setting_name] = self._get_setting(this_setting_name)
                # with exception for display setting (mask) to control display per job
                if(this_setting_name == "display" and default_display_value is not None):
                    # False False => False
                    # False True => False
                    # True True => True
                    # True False => False
                    # that's a "and mask"
                    self.job[this_job_name][this_setting_name] &= default_display_value

        # remove display job if we don't want to display
        if(self.settings["display"] is False):
            while(self.job["displaying"] in self.job_list):
                self.job_list.remove(self.job["displaying"])
            while(self.job["displaying_anatomy"] in self.job_list):
                self.job_list.remove(self.job["displaying_anatomy"])

        # for each job
        for k, this_job in enumerate(self.job_list):
            this_job_name = this_job["job_name"]
            log.info("#%d %s" % (k, this_job_name))
        log.info_line________________________()

        # --- running job list ---
        # exception here: if any concatenate, we need to run the job list in two parts (before and after concatenation in order to reload the processed data)
        log.info("running job list...")
        log.info_line_break()

        # init job stacks
        concatenate_loop = int(self.job["concatenate"] in self.job_list) + 1
        jobs_stack_init = self.job_list[::-1].copy()
        jobs_stack = jobs_stack_init.copy()
        jobs_done_stack_init = []
        jobs_done_stack = jobs_done_stack_init.copy()
        for k in range(concatenate_loop):
            # resuming job stack after concatenate
            jobs_stack_init = jobs_stack.copy()
            jobs_done_stack_init = jobs_done_stack.copy()
            # for each dataset, raw and dicom, the list of data processing functions with the right arguments
            for i, d in enumerate(self.dataset):
                for dtype in ["raw", "dcm"]:
                    this_data = d[dtype]["data"]
                    # if no data to process (= no DCM provided)
                    if(this_data is None):
                        continue

                    log.info_line________________________()
                    log.info("processing [" + this_data.display_label + "]")
                    log.info_line________________________()

                    # run job list for this dataset
                    jobs_stack = jobs_stack_init.copy()
                    jobs_done_stack = jobs_done_stack_init.copy()
                    while(len(jobs_stack) > 0):
                        # job pop
                        job = jobs_stack.pop()
                        # if concatenate, get out of this dataset job loop
                        if(job == self.job["concatenate"]):
                            break

                        # run job on this dataset

                        # nasty exception: phasing using ref on dicoms
                        # dicom data is already phased
                        # it only makes things worse to use ref. data, just disable it
                        if((dtype == "dcm") and (job == self.job["phasing"])):
                            log.warning("changing phasing setting [using_ref_data] to %d for this dataset!" % (not self.settings["no_phasing_using_ref_for_dcm_data"]))
                            job["using_ref_data"] = not self.settings["no_phasing_using_ref_for_dcm_data"]

                        job_result = self._run_job(job, this_data)

                        # push job in the stack of jobs done
                        jobs_done_stack.append(job)

                        # replace with processed signal
                        if(type(job_result) == MRSData2):
                            this_data = job_result
                            # measure SNR/LW after this process?
                            if(self.analyze_enable):
                                # prepare storage
                                if(self.dataset[i][dtype]["analysis_results"] is None):
                                    self.dataset[i][dtype]["analysis_results"] = {}

                                # get job name, renaming if needed
                                this_job_name = job["job_name"]
                                # get list of job names already in analysis list
                                analysis_results_job_names = list(self.dataset[i][dtype]["analysis_results"].keys())
                                # check if any duplicates, meaning a job than was applied multiple times
                                analysis_results_job_names_mask = [j.startswith(this_job_name) for j in analysis_results_job_names]
                                this_job_already_applied_count = analysis_results_job_names_mask.count(True)
                                if(this_job_already_applied_count > 0):
                                    # not the first time we apply this job?
                                    this_job_name = this_job_name + " (#%d)" % (this_job_already_applied_count + 1)

                                # get SNR and LW estimations
                                this_data_snr, this_data_lw, _, _ = self._analyze(this_data, job, jobs_done_stack)

                                # store analysis results
                                self.dataset[i][dtype]["analysis_results"][this_job_name] = {"snr": this_data_snr, "lw": this_data_lw}

                    # we finish running all jobs on a dataset, storing
                    self.dataset[i][dtype]["data"] = this_data

            # if last job was concatenate, that means all datasets were processed and are ready for it
            if(job == self.job["concatenate"]):
                job_name = self.job["concatenate"][0]["name"]
                log.info("%s..." % job_name)
                # attention, concatenate is only done for raw data
                s_concatenated = self.dataset[0]["raw"]["data"]
                for i in range(1, len(self.dataset)):
                    s_concatenated = s_concatenated.concatenate_2d(self.dataset[i]["raw"]["data"])

                # empty and store the single concatenated signal in the dataset dict
                s_concatenated.set_display_label(s_concatenated.display_label + " [CONCATENATED]")
                self.dataset[0]["legend"] = s_concatenated.display_label
                self.dataset[0]["raw"]["data"] = s_concatenated

        # if we did a concatenate job, then keep only the first concatenated dataset at index 0
        if(self.job["concatenate"] in self.job_list):
            self.dataset = [self.dataset[0]]

        # before leaving, analyze ref data if available
        for i, d in enumerate(self.dataset):
            for dtype in ["raw", "dcm"]:
                this_data = d[dtype]["data"]
                # if no data to process (= no DCM provided)
                if(this_data is None):
                    continue

                _, _, this_data_ref_snr, this_data_ref_lw = self._analyze(this_data, job, jobs_done_stack)
                self.dataset[i][dtype]["ref_data_analysis_results"] = {"snr": this_data_ref_snr, "lw": this_data_ref_lw}

        # --- summary final linewidths ---
        if(self.analyze_enable):
            self.display_analyze_results()

        log.info("pipeline terminated!")
        return(self.dataset)

    def get_analyze_results(self):
        """Return analyze results as several numpy vectors eassy to plot.Dataset and job names will be padded to ease terminal output.

        Returns
        -------
        data_label_list : list (n)
            Data legends
        job_label_list : list (n)
            Job names
        snr_raw_list : numpy array (n x m)
            SNR estimated on raw data
        lw_raw_list : numpy array (n x m)
            Peak linewidth estimated on raw data
        snr_dcm_list : numpy array (n x m)
            SNR estimated on reconstructed data (dicom)
        lw_dcm_list : numpy array (n x m)
            Peak linewidth estimated on reconstructed data (dicom)
        snr_ref_raw_list : numpy array (n)
            SNR estimated on raw ref. data
        lw_ref_raw_list : numpy array (n)
            Peak linewidth estimated on raw ref. data
        snr_ref_dcm_list : numpy array (n)
            SNR estimated on reconstructed ref. data (dicom)
        lw_ref_dcm_list : numpy array (n)
            Peak linewidth estimated on reconstructed ref. data (dicom)
        """
        log.debug("converting analyze results to arrays...")

        # build label list
        data_labels = [d["legend"] for d in self.dataset]
        data_label_nchar = len(max(data_labels, key=len))
        data_label_list = [d.ljust(data_label_nchar) for d in data_labels]

        # build job label list
        if(self.dataset[0]["raw"]["analysis_results"] is not None):
            job_labels = self.dataset[0]["raw"]["analysis_results"].keys()
        elif(self.dataset[0]["dcm"]["analysis_results"] is not None):
            job_labels = self.dataset[0]["dcm"]["analysis_results"].keys()
        else:
            job_labels = []

        job_label_nchar = len(max(job_labels, key=len))
        job_label_list = [d.ljust(job_label_nchar) for d in job_labels]

        # init
        snr_raw_list = []
        lw_raw_list = []
        snr_dcm_list = []
        lw_dcm_list = []

        snr_ref_raw_list = []
        lw_ref_raw_list = []
        snr_ref_dcm_list = []
        lw_ref_dcm_list = []

        # for each dataset
        for d in self.dataset:
            # if raw data, find snr and lw estimations and store it
            if(d["raw"]["analysis_results"] is not None):
                snr_raw_list.append([d["raw"]["analysis_results"][j]["snr"] for j in d["raw"]["analysis_results"].keys()])
                lw_raw_list.append([d["raw"]["analysis_results"][j]["lw"] for j in d["raw"]["analysis_results"].keys()])
            else:
                snr_raw_list.append([np.nan] * len(job_label_list))
                lw_raw_list.append([np.nan] * len(job_label_list))

            # if dcm data, find snr and lw estimations and store it
            if(d["dcm"]["analysis_results"] is not None):
                snr_dcm_list.append([d["dcm"]["analysis_results"][j]["snr"] for j in d["dcm"]["analysis_results"].keys()])
                lw_dcm_list.append([d["dcm"]["analysis_results"][j]["lw"] for j in d["dcm"]["analysis_results"].keys()])
            else:
                snr_dcm_list.append([np.nan] * len(job_label_list))
                lw_dcm_list.append([np.nan] * len(job_label_list))

            # if raw ref data, find snr and lw estimations and store it
            if(d["raw"]["ref_data_analysis_results"] is not None):
                snr_ref_raw_list.append([d["raw"]["ref_data_analysis_results"]["snr"]])
                lw_ref_raw_list.append([d["raw"]["ref_data_analysis_results"]["lw"]])
            else:
                snr_ref_raw_list.append([np.nan])
                lw_ref_raw_list.append([np.nan])

            # if dcm ref data, find snr and lw estimations and store it
            if(d["dcm"]["ref_data_analysis_results"] is not None):
                snr_ref_dcm_list.append([d["dcm"]["ref_data_analysis_results"]["snr"]])
                lw_ref_dcm_list.append([d["dcm"]["ref_data_analysis_results"]["lw"]])
            else:
                snr_ref_dcm_list.append([np.nan])
                lw_ref_dcm_list.append([np.nan])

        snr_raw_list = np.array(snr_raw_list)
        lw_raw_list = np.array(lw_raw_list)
        snr_dcm_list = np.array(snr_dcm_list)
        lw_dcm_list = np.array(lw_dcm_list)

        snr_ref_raw_list = np.array(snr_ref_raw_list)
        lw_ref_raw_list = np.array(lw_ref_raw_list)
        snr_ref_dcm_list = np.array(snr_ref_dcm_list)
        lw_ref_dcm_list = np.array(lw_ref_dcm_list)

        return(data_label_list, job_label_list, snr_raw_list, lw_raw_list, snr_dcm_list, lw_dcm_list, snr_ref_raw_list, lw_ref_raw_list, snr_ref_dcm_list, lw_ref_dcm_list)

    def get_final_analyze_results(self):
        """Return final analyze results for each dataset.

        Returns
        -------
        data_label_list : list (n)
            Data legends
        snr_raw_list : numpy array (n)
            SNR estimated on data
        lw_raw_list : numpy array (n)
            Peak linewidth estimated on data
        lw_ref_raw_list : numpy array (n)
            Peak linewidth estimated on ref. data
        """
        log.debug("returning final analyze results...")

        data_label_list, job_label_list, snr_raw_list, lw_raw_list, snr_dcm_list, lw_dcm_list, _, lw_ref_raw_list, _, lw_ref_dcm_list = self.get_analyze_results()

        # only keep final estimations from last job
        snr_raw_list = snr_raw_list[:, -1]
        lw_raw_list = lw_raw_list[:, -1]
        snr_dcm_list = snr_dcm_list[:, -1]
        lw_dcm_list = lw_dcm_list[:, -1]
        lw_ref_raw_list = lw_ref_raw_list[:, -1]
        lw_ref_dcm_list = lw_ref_dcm_list[:, -1]

        # if any raw data estimations are None, replace them by dcm estimations
        snr_raw_list[np.isnan(snr_raw_list)] = snr_dcm_list[np.isnan(snr_raw_list)]
        lw_raw_list[np.isnan(lw_raw_list)] = lw_dcm_list[np.isnan(lw_raw_list)]
        lw_ref_raw_list[np.isnan(lw_ref_raw_list)] = lw_ref_dcm_list[np.isnan(lw_ref_raw_list)]

        return(data_label_list, snr_raw_list, lw_raw_list, lw_ref_raw_list)

    def check_analyze_results(self, ignore_error=False):
        """Check for each dataset that we got a beter reconstruction with raw data than dicom.

        Returns
        -------
        ignore_error : boolean
            Do not raise error if raw data reconstruction gave worst results than DICOM in terms of SNR and /or LW. If (False), just raise a warning.
        """
        log.debug("quality checking analyze results...")

        data_label_list, job_label_list, snr_raw_list, lw_raw_list, snr_dcm_list, lw_dcm_list, _, _, _, _ = self.get_analyze_results()

        # only keep final estimations from last job
        snr_raw_list = snr_raw_list[:, -1]
        lw_raw_list = lw_raw_list[:, -1]
        snr_dcm_list = snr_dcm_list[:, -1]
        lw_dcm_list = lw_dcm_list[:, -1]

        error_str = "final SNR and/or LW is worst for raw data than for dicom for dataset(s): "
        n_bad_reco = 0
        # compare dcm analyze results to raw data results, when possible
        for this_dataset_label, this_snr_raw, this_lw_raw, this_snr_dcm, this_lw_dcm in zip(data_label_list, snr_raw_list, lw_raw_list, snr_dcm_list, lw_dcm_list):
            if(not np.isnan(this_snr_raw) and not np.isnan(this_snr_dcm) and not np.isnan(this_lw_raw) and not np.isnan(this_lw_dcm)):
                if((this_snr_raw < this_snr_dcm) or (this_lw_raw > this_lw_dcm)):
                    # the raw data reconstruction was worst than the dicom :(
                    n_bad_reco = n_bad_reco + 1
                    error_str = error_str + ("[%s], " % this_dataset_label)

        if(n_bad_reco > 0):
            if(self.settings["raise_error_on_badreco"] and not ignore_error):
                log.error(error_str[:-2])
            else:
                log.warning(error_str[:-2])

    def display_analyze_results(self, fig_index="Quality check results", save2file=False):
        """Print final SNR and peak-linewidth for each dataset. Plot bargraph showing evolution of SNR and linewidth during data processing (to check that a job did not destroy the data for example!) and compare with dicom when possible.

        Parameters
        ----------
        fig_index: int or string
            Figure handle
        save2file: boolean
            Save figure to png file in current folder
        """
        log.info("displaying SNR and linewidth final results...")

        # get analyze data as list and np arrays
        data_label_list, job_label_list, snr_raw_list, lw_raw_list, snr_dcm_list, lw_dcm_list, snr_ref_raw_list, lw_ref_raw_list, snr_ref_dcm_list, lw_ref_dcm_list = self.get_analyze_results()

        # terminal output

        # for each dataset (line), print final SNR/LW columns
        # if DCM dataset was included, add SNR/LW for DCM too
        log.info_line________________________()
        print("RAW / DCM dataset ".ljust(len(data_label_list[0])) + "\t" + "SNR (u.a)".ljust(20) + "\t" + "LW (Hz)".ljust(20) + "\t" + "ref. data LW (Hz)".ljust(20), flush=True)
        print("", flush=True)
        for i, this_dataset_label in enumerate(data_label_list):
            this_dataset_lastjob_snr_raw = snr_raw_list[i][-1]
            this_dataset_lastjob_snr_dcm = snr_dcm_list[i][-1]
            this_dataset_lastjob_snr_str = ("%.2f / %.2f" % (this_dataset_lastjob_snr_raw, this_dataset_lastjob_snr_dcm)).ljust(20)

            this_dataset_lastjob_lw_raw = lw_raw_list[i][-1]
            this_dataset_lastjob_lw_dcm = lw_dcm_list[i][-1]
            this_dataset_lastjob_lw_str = ("%.2f / %.2f" % (this_dataset_lastjob_lw_raw, this_dataset_lastjob_lw_dcm)).ljust(20)

            this_dataset_ref_lw_raw = lw_ref_raw_list[i]
            this_dataset_ref_lw_dcm = lw_ref_dcm_list[i]
            this_dataset_ref_lw_str = ("%.2f / %.2f" % (this_dataset_ref_lw_raw, this_dataset_ref_lw_dcm)).ljust(20)
            print("%s\t%s\t%s\t%s" % (this_dataset_label, this_dataset_lastjob_snr_str, this_dataset_lastjob_lw_str, this_dataset_ref_lw_str), flush=True)
        log.info_line________________________()

        # display SNR/LW evolution for all data and all jobs
        if(self.settings["display"]):
            fig = plt.figure(fig_index)
            fig.clf()
            fig.suptitle(fig_index)
            axs = fig.subplots(2, 1, sharex='row')

            # prepare bars
            nBars = len(data_label_list)
            width = max(-0.15 * nBars + 0.6, 0.1)  # bars get thinner if more data
            pos_bars = np.arange(len(job_label_list))
            pos_shift = np.linspace(-width * (nBars - 1) / 2.0, +width * (nBars - 1) / 2.0, nBars)
            snr_ylim = [+np.inf, -np.inf]
            lw_ylim = [+np.inf, -np.inf]

            # iterate in lists and plot bars
            for i, (this_dataset_label, pos) in enumerate(zip(data_label_list, pos_shift)):
                # snr list for this dataset, raw data
                this_dataset_snr_list = snr_raw_list[i][:]
                snr_ylim[0] = min(snr_ylim[0], np.min(this_dataset_snr_list))
                snr_ylim[1] = max(snr_ylim[1], np.max(this_dataset_snr_list))
                axs[0].bar(pos_bars + pos, this_dataset_snr_list, width, label=this_dataset_label)
                # snr list for this dataset, dcm data
                this_dataset_snr_list = snr_dcm_list[i][:]
                snr_ylim[0] = min(snr_ylim[0], np.min(this_dataset_snr_list))
                snr_ylim[1] = max(snr_ylim[1], np.max(this_dataset_snr_list))
                axs[0].bar(pos_bars + pos, this_dataset_snr_list, width, fill=False)

                # lw list for this dataset, raw data
                this_dataset_lw_list = lw_raw_list[i][:]
                lw_ylim[0] = min(lw_ylim[0], np.min(this_dataset_lw_list))
                lw_ylim[1] = max(lw_ylim[1], np.max(this_dataset_lw_list))
                axs[1].bar(pos_bars + pos, this_dataset_lw_list, width, label=this_dataset_label)
                # lw list for this dataset, dcm data
                this_dataset_lw_list = lw_dcm_list[i][:]
                lw_ylim[0] = min(lw_ylim[0], np.min(this_dataset_lw_list))
                lw_ylim[1] = max(lw_ylim[1], np.max(this_dataset_lw_list))
                axs[1].bar(pos_bars + pos, this_dataset_lw_list, width, fill=False)

            # add fake DCM bars to have it in the legend
            axs[0].bar(pos_bars, this_dataset_lw_list * 0.0, width, fill=False, label="DCM")

            # snr bargraph
            axs[0].set_ylabel("SNR (u.a)")
            axs[0].set_xticks(pos_bars)
            axs[0].set_xticklabels([])
            axs[0].set_ylim(np.add(snr_ylim, [-5, +5]))
            axs[0].grid('on')
            axs[0].legend()

            # lw bargraph
            axs[1].set_ylabel("Estimated linewidth (Hz)")
            axs[1].set_xticks(pos_bars)
            axs[1].set_xticklabels(job_label_list, rotation=45)
            axs[1].set_ylim(np.add(lw_ylim, [-3, +3]))
            axs[1].grid('on')

            fig.subplots_adjust(bottom=0.2)
            fig.show()

            # save figure as png if needed
            if(save2file):
                this_filename = re.sub('[^\w\-_\. ]', '_', fig_index)
                this_filename = this_filename.replace(" ", "_")
                while(this_filename[0] == "_"):
                    this_filename = this_filename[1:]

                this_filename += ".png"
                fig.savefig(this_filename)

    def display_final_data(self, fig_index="Final spectra display", save2file=False):
        """
        Plot final processed datasets.

        Parameters
        ----------
        fig_index: int or string
            Figure handle
        save2file: boolean
            Save figure to png file in current folder
        """
        log.debug("displaying final processed data...")

        # get display job
        disp_job = self.job["displaying"]

        # change fig index
        disp_job["fig_index"] = fig_index

        # run a diplay job for each dataset
        for d in self.dataset:
            for dtype in ["raw", "dcm"]:
                this_data = d[dtype]["data"]
                # if no data
                if(this_data is None):
                    continue
                # run job on this dataset
                fig = self._run_job(disp_job, this_data)

        # save figure as png if needed
        if(save2file):
            this_filename = re.sub('[^\w\-_\. ]', '_', fig_index)
            this_filename = this_filename.replace(" ", "_")
            while(this_filename[0] == "_"):
                this_filename = this_filename[1:]

            this_filename += ".png"
            fig.savefig(this_filename)

    def load_template(self, template_filename):
        """
        Load reco pipeline attributes from template previously saved.

        Parameters
        ----------
        template_filename: string
            Path to template file
        """
        log.debug("load pipeline from template file [%s]..." % template_filename)
        # check if exist
        if(os.path.isfile(template_filename)):
            # now open pkl file
            with open(template_filename, 'rb') as f:
                [reco_pipe_template] = pickle.load(f)

            # copy attributes
            self.settings = reco_pipe_template.settings
            self.job = reco_pipe_template.job
            self.job_list = reco_pipe_template.job_list
            self.analyze_job_list = reco_pipe_template.analyze_job_list
            self.analyze_enable = reco_pipe_template.analyze_enable
        else:
            log.error("file not found [%s]!" % template_filename)

    def save_template(self, template_filename):
        """
        Save current pipeline object as template.

        Parameters
        ----------
        template_filename: string
            Pth to template file
        """
        log.debug("saving current pipeline to [%s]..." % template_filename)
        # save to pkl file
        with open(template_filename, 'wb') as f:
            pickle.dump([self], f)

    def _to_dataframe(self, prefix_str="reco_"):
        """
        Convert the object's attributes to dataframe.

        Parameters
        ----------
        prefix_str : string
            Prefix string to add to column names

        Returns
        -------
        df : Dataframe
            With all the datasets and parameters from this pipeline
        """
        log.debug("converting to dataframe...")

        # use pandas json_normalize function to flatten all nested stuff (magic!)
        # this takes care of the job attribute too
        df_p = pd.json_normalize(vars(self), sep='_')

        # need a specific call for the dataset attribute
        df_dataset = pd.json_normalize(self.dataset, sep='_')
        df_dataset = df_dataset.add_prefix("dataset_")
        del df_p["dataset"]

        # and need a specific call for the MRSData2 types
        df_dataset_raw_data_list = [s.to_dataframe(True, "dataset_raw_data_")for s in df_dataset["dataset_raw_data"] if(s is not None)]
        df_dataset_dcm_data_list = [s.to_dataframe(True, "dataset_dcm_data_") for s in df_dataset["dataset_dcm_data"] if(s is not None)]

        df = df_dataset.reset_index()

        if(len(df_dataset_raw_data_list) > 0):
            df_dataset_raw_data = pd.concat(df_dataset_raw_data_list, axis=0)
            del df_dataset["dataset_raw_data"]
            # append columns
            df = pd.concat([df, df_dataset_raw_data.reset_index()], axis=1)
        if(len(df_dataset_dcm_data_list) > 0):
            df_dataset_dcm_data = pd.concat(df_dataset_dcm_data_list, axis=0)
            del df_dataset["dataset_dcm_data"]
            # append columns
            df = pd.concat([df, df_dataset_dcm_data.reset_index()], axis=1)

        # append columns
        df = pd.concat([df, df_p.loc[df_p.index.repeat(len(df_dataset))].reset_index()], axis=1)

        # set index (prefer raw data if available)
        if("dataset_raw_data_file_hash" in df):
            df = df.set_index("dataset_raw_data_file_hash")
        elif("dataset_dcm_data_file_hash" in df):
            df = df.set_index("dataset_dcm_data_file_hash")
        else:
            log.error("no data file hash found in datasets! :(")

        del df["index"]

        # prefix
        df = df.add_prefix(prefix_str)

        return(df)

    def save_datasets(self):
        """Save corresponding this pipeline, all its parameters and its datasets to a dataframe stored on the disk."""
        log.debug("saving to disk...")

        # check it we have specified a file path
        if(self.settings["storage_file"] is None):
            log.error("no pkl file path specified to save to disk!")
            return()

        # first, convert this pipeline and all to a df
        this_df = self._to_dataframe()
        # add timestamp
        this_df["timestamp"] = datetime.now()

        # load current df from file, if exist
        pkl_filepath = self.settings["storage_file"]
        if(os.path.isfile(pkl_filepath)):
            # load
            df = pd.read_pickle(pkl_filepath)
            # remove rows for which we have new data
            df = df.drop(df.loc[df.index.isin(this_df.index)].index)
            # append new rows
            df = pd.concat([df, this_df])

            # save back
            df.to_pickle(pkl_filepath)
        else:
            # no file
            this_df.to_pickle(pkl_filepath)


def remove_grids_from_all_figs():
    """Remove the grid in all axes for all open figures."""
    figs = [manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]

    for f in figs:
        for a in f.axes:
            a.grid(False)


def reco_spatial_select_profile(dcm_folders_list, legends_list, analyze_selectivity_range_list=[[800, 3550], [-10600, -7800], [-3650, 1850]]):
    """
    Reconstruct the VOI selection profile from data acquired using the 'VOI trace' mode (available in CMRR's MRS sequences).

    Parameters
    ----------
    dcm_folders_list: list
        List of folder paths containing DICOM files
    legends_list : list
        List of legends
    analyze_selectivity_range_list : list
        List of ranges (point index) in the X, Y and Z direction

    Returns
    -------
    analyze_selectivity_list : list
        Analysis results (amount of signal IN and OUT of the ranges)
    """
    # init
    data_list = []
    xdata_list = []
    ydata_list = []
    zdata_list = []
    xdata_axis_list = []
    ydata_axis_list = []
    zdata_axis_list = []

    # read the 3 VOI traces per dataset
    log.info_line________________________()
    log.info("reading data files...")
    log.info_line________________________()
    for f in dcm_folders_list:
        # read data
        log.info("looking in folder: ")
        log.info(f)
        log.info("reading the 3 dicom files...")
        sx = suspect.io.load_siemens_dicom(f + "/original-primary_e09_0001.dcm")
        sy = suspect.io.load_siemens_dicom(f + "/original-primary_e09_0002.dcm")
        sz = suspect.io.load_siemens_dicom(f + "/original-primary_e09_0003.dcm")
        data_list.append([sx, sy, sz])
        # normalize
        sx_spectrum = np.abs(sx.spectrum()) / np.abs(sx.spectrum()).max()
        sy_spectrum = np.abs(sy.spectrum()) / np.abs(sy.spectrum()).max()
        sz_spectrum = np.abs(sz.spectrum()) / np.abs(sz.spectrum()).max()
        # store
        xdata_list.append(sx_spectrum)
        ydata_list.append(sy_spectrum)
        zdata_list.append(sz_spectrum)
        xdata_axis_list.append(sx.frequency_axis())
        ydata_axis_list.append(sy.frequency_axis())
        zdata_axis_list.append(sz.frequency_axis())

    log.info_line________________________()

    # analysis
    log.info_line________________________()
    log.info("evaluating selectivity using a " + str(analyze_selectivity_range_list) + "ppm ranges...")
    log.info_line________________________()
    analyze_selectivity_list = np.zeros([len(dcm_folders_list), 3, 2])
    for k, (sx, sy, sz, sx_ax, sy_ax, sz_ax, leg) in enumerate(zip(xdata_list, ydata_list, zdata_list, xdata_axis_list, ydata_axis_list, zdata_axis_list, legends_list)):

        # find ppm indexes corresponding to ppm values
        sx_ithreshold_l = np.argmin(np.abs(sx_ax - analyze_selectivity_range_list[0][0]))
        sx_ithreshold_r = np.argmin(np.abs(sx_ax - analyze_selectivity_range_list[0][1]))
        sy_ithreshold_l = np.argmin(np.abs(sy_ax - analyze_selectivity_range_list[1][0]))
        sy_ithreshold_r = np.argmin(np.abs(sy_ax - analyze_selectivity_range_list[1][1]))
        sz_ithreshold_l = np.argmin(np.abs(sz_ax - analyze_selectivity_range_list[2][0]))
        sz_ithreshold_r = np.argmin(np.abs(sz_ax - analyze_selectivity_range_list[2][1]))

        # evaluate signal inside voxel
        sx_in = np.trapz(sx[sx_ithreshold_l:sx_ithreshold_r], sx_ax[sx_ithreshold_l:sx_ithreshold_r])
        sy_in = np.trapz(sy[sy_ithreshold_l:sy_ithreshold_r], sy_ax[sy_ithreshold_l:sy_ithreshold_r])
        sz_in = np.trapz(sz[sz_ithreshold_l:sz_ithreshold_r], sz_ax[sz_ithreshold_l:sz_ithreshold_r])

        # put to zero the inside
        sx_masked = sx.copy()
        sx_masked[sx_ithreshold_l:sx_ithreshold_r] = 0
        sy_masked = sy.copy()
        sy_masked[sy_ithreshold_l:sy_ithreshold_r] = 0
        sz_masked = sz.copy()
        sz_masked[sz_ithreshold_l:sz_ithreshold_r] = 0

        # evaluate signal outside voxel
        sx_out = np.trapz(sx_masked, sx_ax)
        sy_out = np.trapz(sy_masked, sy_ax)
        sz_out = np.trapz(sz_masked, sz_ax)

        log.info("selectivity results for [" + leg + "]...")
        log.info(" [X] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sx_in, sx_out))
        log.info(" [Y] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sy_in, sy_out))
        log.info(" [Z] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sz_in, sz_out))

        # store
        analyze_selectivity_list[k, :, 0] = [sx_in, sy_in, sz_in]
        analyze_selectivity_list[k, :, 1] = [sx_out, sy_out, sz_out]

    log.info_line________________________()

    # display
    log.info_line________________________()
    log.info("displaying the 3 spatial profiles...")
    log.info_line________________________()
    fig = plt.figure("Spatial selection profiles")
    fig.clf()
    axs = fig.subplots(2, 3)
    fig.suptitle("Spatial selection profiles")
    for (sx, sy, sz, sx_ax, sy_ax, sz_ax, leg) in zip(xdata_list, ydata_list, zdata_list, xdata_axis_list, ydata_axis_list, zdata_axis_list, legends_list):

        axs[0, 0].plot(sx_ax, sx, linewidth=1, label=leg)
        axs[0, 0].set_xlabel('frequency dispersion (Hz)')
        axs[0, 0].grid('on')

        axs[1, 0].plot(sx_ax, sx, linewidth=1, label=leg)
        axs[1, 0].set_xlabel('frequency dispersion (Hz)')
        axs[1, 0].set_ylabel('X spatial profile')
        axs[1, 0].grid('on')

        axs[0, 1].plot(sy_ax, sy, linewidth=1, label=leg)
        axs[0, 1].set_xlabel('frequency dispersion (Hz)')
        axs[0, 1].grid('on')

        axs[1, 1].plot(sy_ax, sy, linewidth=1, label=leg)
        axs[1, 1].set_xlabel('frequency dispersion (Hz)')
        axs[1, 1].set_ylabel('Y spatial profile')
        axs[1, 1].grid('on')

        axs[0, 2].plot(sz_ax, sz, linewidth=1, label=leg)
        axs[0, 2].set_xlabel('frequency dispersion (Hz)')
        axs[0, 2].grid('on')

        axs[1, 2].plot(sz_ax, sz, linewidth=1, label=leg)
        axs[1, 2].set_xlabel('frequency dispersion (Hz)')
        axs[1, 2].set_ylabel('Z spatial profile')
        axs[1, 2].grid('on')

    # display ppm range lines
    for i in range(3):
        axs[0, i].axvline(x=analyze_selectivity_range_list[i][0], linestyle='--')
        axs[0, i].axvline(x=analyze_selectivity_range_list[i][1], linestyle='--')
        axs[1, i].axvline(x=analyze_selectivity_range_list[i][0], linestyle='--')
        axs[1, i].axvline(x=analyze_selectivity_range_list[i][1], linestyle='--')

    axs[1, 2].legend()
    fig.subplots_adjust()
    fig.show()

    return(analyze_selectivity_list)


def get_dw_reco_pipeline(data, template_name=None, legend=""):
    """
    Generate a reconstruction pipeline object to process this DW-MRS data. In short, it reshapes and split the current dataset into several datasets corresponding to different directions and bvalues, in order to process each of them in a pipeline. This method is quite weird and experimental for now. A lot of stuff missing. And that is why it is here, out of any class.

    Parameters
    ----------
    data : MRSData2 object
        DW data to process
    template_name : string
        Reco pipeline template file to use
    legend : string
        Some legend to use for this dataset

    Returns
    -------
    p : reco.pipeline object
        Ready to run reconstruction pipeline for DW-MRS data
    """
    log.debug("generating a reconstruction pipeline to process DW-MRS data...")

    if(not hasattr(data.sequence, "directions") or not hasattr(data.sequence, "bvalues") or not hasattr(data.sequence, "n_b0")):
        log.error("sorry but this dataset was not acquired with a DW-MRS sequence.")

    # get the dw parameters from the sequence
    n_averages = data.sequence.na
    n_directions = len(data.sequence.directions)
    n_bvalues = len(data.sequence.bvalues)
    n_b0 = data.sequence.n_b0

    # number of dw scans including b0
    n_diff = (n_directions * n_bvalues) + n_b0
    # number of channels
    n_chan = data.shape[1]
    # number of time points
    n_pts = data.shape[2]

    # 1st reshape
    s4d = np.reshape(data, [n_averages, n_diff, n_chan, n_pts])

    # extract b0 data temporarly
    s4d_b0 = s4d[:, 0, :, :]
    # and remove it from dataset (considering it was acquired first!)
    s4d = s4d[:, 1:, :, :]

    # 2nd reshape: separate bval from dir
    # TODO: this 2nd reshape is a shot in the dark!!!
    s5d = np.reshape(s4d, [n_averages, n_directions, n_bvalues, n_chan, n_pts])
    s5d_b0 = np.reshape(s4d_b0, [n_averages, 1, 1, n_chan, n_pts])

    # now prepape reco pipeline
    p = pipeline(template_name)

    # edit the dataset list
    i_dataset = 0

    # store b0
    s = data.inherit(np.squeeze(s5d_b0))
    # edit hash
    s._data_file_hash = s._data_file_hash + "_b0"
    # edit sequence parameters
    s.sequence.directions = [0, 0, 0]
    s.sequence.bvalues = 0
    # edit label
    s.set_display_label(legend + " b0")
    p.dataset[i_dataset]["legend"] = legend + " b0"
    # store
    p.dataset[i_dataset]["raw"]["data"] = s

    i_dataset = 1
    for this_dir_index, this_dir in enumerate(data.sequence.directions):
        for this_bval_index, this_bval in enumerate(data.sequence.bvalues):
            # build MRSData2 object
            s_extracted = s5d[:, this_dir_index, this_bval_index, :, :]
            s = data.inherit(s_extracted)
            # edit hash
            s._data_file_hash = s._data_file_hash + "_d%db%d" % (this_dir_index + 1, this_bval_index + 1)
            # edit sequence parameters
            s.sequence.directions = this_dir
            s.sequence.bvalues = this_bval
            # edit label
            p.dataset[i_dataset]["legend"] = legend + " d%d b%d" % (this_dir_index + 1, this_bval_index + 1)
            s.set_display_label(legend + " d%d b%d" % (this_dir_index, this_bval_index + 1))
            # store
            p.dataset[i_dataset]["raw"]["data"] = s

            # next
            i_dataset = i_dataset + 1

    return(p)
