#! / usr / bin / env python3
# -*- coding: utf-8 -*-
"""
Three classes' definition in here.

    * a MRSData2 class with a bunch of methods based on the suspect module to deal with SIEMENS 7T MRS data
    * a pipeline class to run the reconstruction process on a bunch of acquired data
    * a voi_pipeline class to reconstruct VOI profile data acquired using the debug mode on CMRR's MRS sequences

@author: Tangi Roussel
"""

import suspect
import numpy as np
from scipy import signal
import scipy.io as sio
import matplotlib.pylab as plt
import matplotlib._pylab_helpers
from datetime import datetime
import pickle
import os
from enum import Enum
from mrs import io
from mrs import sim
from mrs import log
from mrs import paths as default_paths

import pdb

max_number_datasets_pipeline = 1000


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

    def __new__(cls, data_filepath, physio_log_file=None, obj=None, dt=None, f0=None, te=None, tr=None, ppm0=None, voxel_dimensions=None, transform=None, metadata=None, data_ref=None, label="", offset_display=0.0, patient={}, sequence_obj=None, noise_level=None, data_rejection=None, data_file_hash=None, is_concatenated=None, is_rawdata=None):
        """
        Construct a MRSData2 object that inherits of Suspect's MRSData class. In short, the MRSData2 class is a copy of MRSData + my custom methods for post-processing. To create a MRSData2 object, you need give a path that points to a SIEMENS DICOM or a SIEMENS TWIX file.

        Parameters
        ----------
        data_filepath: string
            Full absolute file path pointing to the stored signal (DCM or TWIX file) or the folder assuming that a dcm file named "original-primary_e09_0001.dcm" is stored inside.
        physio_log_file : string
            Full absolute file path pointing to a IDEA VB17 respiratory log file
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
        # read header
        mfr = io.SIEMENS_data_file_reader(data_filepath)
        # read data and get a suspect MRSData object
        MRSData_obj = mfr.read_data()
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
        self._display_label = getattr(obj, 'display_label', None)
        self._display_offset = getattr(obj, 'display_offset', 0.0)
        self._patient = getattr(obj, 'patient', None)
        self._physio_file = getattr(obj, 'physio_file', None)
        self._sequence = getattr(obj, 'sequence', None)
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
        obj2._display_label = getattr(self, 'display_label', None)
        obj2._display_offset = getattr(self, 'display_offset', 0.0)
        obj2._patient = getattr(self, 'patient', None)
        obj2._physio_file = getattr(self, 'physio_file', None)
        obj2._sequence = getattr(self, 'sequence', None)
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

    def _analyze_peak_1d(self, ppm_range):
        """
        Find peak in specific ppm range using magnitude mode and return stuff.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        ppm_range : list [2]
            Range in ppm used for peak searching

        Returns
        -------
        peak_index : float
            Index position of the peak
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
        ppm = self.frequency_axis_ppm()
        sf = self.spectrum()
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
        peak_lw = dppm * self.f0

        # peak segment
        ippm_half_peak = np.arange(ippm_min, ippm_max)
        ppm_seg = ppm[ippm_half_peak]
        peak_seg = sf[ippm_half_peak]

        return(peak_index, peak_ppm, peak_val, peak_lw, ppm_seg, peak_seg)

    def _analyze_peak_2d(self, peak_range=[4.5, 5]):
        """
        Analyze a peak in the spectrum by estimating its amplitude, linewidth, frequency shift and phase for each average.

        * Works only with a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : array [2]
            Range in ppm used to analyze peak phase when no reference signal is specified

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

        # first, find peak of interest in range, just to check
        s_avg = np.mean(self, axis=0)
        _, peak_ppm, _, _, _, _ = s_avg._analyze_peak_1d(peak_range)
        log.debug("found peak of interest at %0.2fppm!" % peak_ppm)

        # for each average in moving averaged data
        peak_trace = np.zeros([self.shape[0], 4])
        pbar = log.progressbar("analyzing", self.shape[0])
        for a in range(0, self.shape[0]):

            # call 1D peak analysis
            peak_index, peak_ppm, peak_val, peak_lw, _, _ = self[a, :]._analyze_peak_1d(peak_range)

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
        peak_trace_rel2mean = np.zeros([self.shape[0], 4])
        peak_trace_rel2mean[:, 0] = peak_trace[:, 0] / peak_trace[:, 0].mean() * 100 - 100
        peak_trace_rel2mean[:, 1] = peak_trace[:, 1] - peak_trace[:, 1].mean()
        peak_trace_rel2mean[:, 2] = peak_trace[:, 2] - peak_trace[:, 2].mean()
        peak_trace_rel2mean[:, 3] = peak_trace[:, 3] - peak_trace[:, 3].mean()

        # relative to 1st pt
        peak_trace_rel2firstpt = np.zeros([self.shape[0], 4])
        peak_trace_rel2firstpt[:, 0] = peak_trace[:, 0] / peak_trace[0, 0] * 100 - 100
        peak_trace_rel2firstpt[:, 1] = peak_trace[:, 1] - peak_trace[0, 1]
        peak_trace_rel2firstpt[:, 2] = peak_trace[:, 2] - peak_trace[0, 2]
        peak_trace_rel2firstpt[:, 3] = peak_trace[:, 3] - peak_trace[0, 3]

        pbar.finish("done")
        return(peak_trace, peak_trace_rel2mean, peak_trace_rel2firstpt)

    def analyze_cs_displacement(self):
        """
        Calculate chemical shift displacement in three directions.

        Returns
        -------
        csd : numpy array of floats
            Estimated CS displacement in %/ppm
        """
        # init
        log.debug("estimating chemical shift displacement error for [%s]..." % self.display_label)
        csd = [None, None, None]
        df_abs_Hz = self.f0  # yes, 1ppm==f0[MHz]/1e6=297Hz at 7T

        if(type(self.sequence) == sim.mrs_seq_eja_svs_slaser):
            # assuming X-Y-Z is done with 90-180-180
            log.info("estimating CS displacement for semiLASER: assuming (90x)-(180y)-(180z)!...")

            # X selection done with asymmetric 90° pulse
            # we do not know much about this pulse. We can only say it is 3.4kHz large if the duration is 2ms
            log.debug("estimating CS displacement for (90x): this pulse is the weird asymmetric one, we have no idea what it is exactly!")
            if(self.pulse_laser_exc_length == 2000.0):
                log.debug("since its duration is 2ms here, we assume, according to Oz & Tkac, MRM 65:901-910 (2011), that its bandwidth is 3.4kHz.")
                bw_x_Hz = 3400.0
                grad_x_Hz_m = bw_x_Hz / (self.voxel_size[0] * 0.001)
                d_x_m = 1000.0 * df_abs_Hz / grad_x_Hz_m
                d_x_prct = d_x_m / self.voxel_size[0] * 100.0
            else:
                log.debug("since its duration is not 2ms here, we do not know its bandwith. Therefore, no way to calculate the CS displacement for this axis, sorry ;)")
                d_x_m = None

            # Y selection done with 180°
            bw_y_Hz = self.pulse_laser_rfc_r / (self.pulse_laser_rfc_length / 1000000.0)
            grad_y_Hz_m = bw_y_Hz / (self.voxel_size[1] * 0.001)
            d_y_m = 1000.0 * df_abs_Hz / grad_y_Hz_m
            d_y_prct = d_y_m / self.voxel_size[1] * 100.0

            # Z selection done with 180°
            bw_z_Hz = self.pulse_laser_rfc_r / (self.pulse_laser_rfc_length / 1000000.0)
            grad_z_Hz_m = bw_z_Hz / (self.voxel_size[2] * 0.001)
            d_z_m = 1000.0 * df_abs_Hz / grad_z_Hz_m
            d_z_prct = d_z_m / self.voxel_size[2] * 100.0

            csd = np.array([d_x_prct, d_y_prct, d_z_prct])
        else:
            log.warning("no idea how to calculate CS displacement for this sequence...")

        return(csd)

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
        log.debug("zero-filling [%s]..." % self.display_label)
        s = self.copy()
        nZeros = nPoints_final - s.np
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

            fig = plt.figure(100)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='row', sharey='row')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_zerofill_nd")
            fig.suptitle("zero-filling [%s]" % self.display_label)

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

    def analyze_noise_nd(self, n_pts=100):
        """
        Measure noise level in time domain and store it in the "noise_level" attribute. This is usefull to keep track of the original noise level for later use, CRB normalization durnig quantification for example.

        * Works with multi-dimensional signals.
        * Returns a multi-dimensional signal.

        Parameters
        ----------
        n_pts : int
            Apodization factor in Hz

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
        # we really want real noise, not zeros from zero-filling
        s_nonzero_mask = (s_real != 0.0)
        s_analyze = s_real[s_nonzero_mask]
        # now take the last 100 points
        noise_lev = np.std(s_analyze[-n_pts:-1])
        log.info("noise level = %.2E" % noise_lev)

        # changing noise level attribute
        log.debug("updating noise level...")
        s._noise_level = noise_lev

        # if any ref data available, we crop it too (silently)
        if(s.data_ref is not None):
            s.data_ref = s.data_ref.analyze_noise_nd(n_pts)

        return(s)

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
            # prepare subplots
            fig = plt.figure(110)
            fig.clf()
            axs = fig.subplots(2, 3, sharex='col')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_phase_3d")
            fig.suptitle("phasing [%s]" % self.display_label)

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
                _, peak_ppm, peak_val, _, _, _ = s_avg._analyze_peak_1d(peak_range)
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
                    _, peak_ppm, peak_val, _, _, _ = this_s._analyze_peak_1d(peak_range)
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
        moving_Naverages_half = int((nAvgWindow - 1) / 2)
        for a in range(0, s.shape[0]):
            ia = max(0, a - moving_Naverages_half)
            ib = min(s.shape[0], a + moving_Naverages_half + 1)
            s_ma[a, :] = np.mean(s[ia:ib, :], axis=0)

        return(s_ma)

    def analyze_physio_2d(self, peak_range=[4.5, 5], delta_time_range=1000.0, display=False):
        """
        Analyze the physiological signal and try to correlate it to a peak amplitude, linewidth, frequency shift and phase variations.

        * Works only with a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : array [2]
            Range in ppm used to analyze peak phase when no reference signal is specified
        delta_time_range : float
            Range in ms used to correlate / match the NMR and the physiological signal. Yes, since we are not really sure of the start timestamp we found in the TWIX header, we try to match perfectly the two signals.
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
        peak_prop_abs, _, _ = self.correct_zerofill_nd()._analyze_peak_2d(peak_range)

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
            fig = plt.figure(120)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_physio_2d_1")
            fig.suptitle("analyzing physiological signals for [%s]" % self.display_label)

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
            fig = plt.figure(121)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_physio_2d_2")
            fig.suptitle("analyzing physiological signals for [%s]" % self.display_label)

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
            fig = plt.figure(122)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_physio_2d_3")
            fig.suptitle("analyzing physiological signals for [%s]" % self.display_label)

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
            fig = plt.figure(123)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_physio_2d_4")
            fig.suptitle("analyzing physiological signals for [%s]" % self.display_label)

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

    def correct_analyze_and_reject_2d(self, peak_range=[4.5, 5], moving_Naverages=1, peak_properties_ranges={"amplitude (%)": None, "linewidth (Hz)": [5.0, 30.0], "chemical shift (ppm)": 0.5, "phase std. factor (%)": 60.0}, peak_properties_rel2mean=True, auto_method_list=None, auto_adjust_allowed_snr_change=0.0, display=False, display_range=[1, 6]):
        """
        Analyze peak in each average in terms intensity, linewidth, chemical shift and phase and reject data if one of these parameters goes out of the min / max bounds. Usefull to understand what the hell went wrong during your acquisition when you have the raw data (TWIX) and to try to improve things a little. You can choose to set the bounds manually or automatically based on a peak property (amplitude, linewidth, frequency, phase). And you can run several automatic adjusment methods, the one giving the highest SNR and/or the lowest peak linewidth will be selected.

        * Works only with a 2D [averages,timepoints] signal.
        * Returns a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : list
            Range in ppm used to analyze peak phase when no reference signal is specified
        moving_Naverages : int
            Number of averages to perform when using moving average, need to be an odd number
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
            log.info("This is the 1st time we perform data rejection on this signal!")
            iround_data_rej = 1
        else:
            iround_data_rej = len(self.data_rejection) + 1
            log.info("This is the %dth time we perform data rejection on this signal!" % iround_data_rej)

        # estimate initial SNR and linewidth
        old_level = log.getLevel()
        log.setLevel(log.ERROR)
        initial_snr, _, _ = s.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd().analyze_snr_1d(peak_range)
        initial_lw = s.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd().analyze_linewidth_1d(peak_range)
        log.setLevel(old_level)
        log.info("* Pre-data-rejection SNR = %.2f" % initial_snr)
        log.info("* Pre-data-rejection linewidth = %.2f Hz" % initial_lw)

        # build moving averaged data
        s_ma = self.correct_zerofill_nd()._build_moving_average_data_2d(moving_Naverages)

        # perform peak analysis
        peak_prop_abs, peak_prop_rel2mean, peak_prop_rel2firstpt = s_ma._analyze_peak_2d(peak_range)

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

            display_axes_ready = [False, False]
            properties_names = list(peak_properties_ranges.keys())
            auto_method_final_snr_list = np.array([0.0] * 4)
            auto_method_final_lw_list = np.array([np.inf] * 4)
            peak_prop_min_auto_res = peak_prop_min.copy()
            peak_prop_max_auto_res = peak_prop_max.copy()

            for this_auto_method in auto_method_list:

                if(this_auto_method == data_rejection_method.AUTO_LINEWIDTH):
                    this_prop_min = max(peak_prop_analyze[:, 1].min(), peak_properties_ranges_list[1][0])
                    this_prop_max = peak_prop_analyze[:, 1].max()
                else:
                    this_prop_min = np.abs(peak_prop_analyze[:, this_auto_method.value]).min()
                    this_prop_max = np.abs(peak_prop_analyze[:, this_auto_method.value]).max()

                # adjust the number of tries / number of steps depending on criteria
                if(this_auto_method == data_rejection_method.AUTO_AMPLITUDE):
                    this_prop_step = 1.0
                elif(this_auto_method == data_rejection_method.AUTO_LINEWIDTH):
                    this_prop_step = 1.0
                elif(this_auto_method == data_rejection_method.AUTO_FREQUENCY):
                    this_prop_step = 0.001
                elif(this_auto_method == data_rejection_method.AUTO_PHASE):
                    this_prop_step = 0.01
                else:
                    log.error("upsyy! I am not aware of this automatic data rejection method: " + str(this_auto_method))

                # generate a range of criteria value to test, using the step
                this_prop_range = np.arange(this_prop_min, this_prop_max, this_prop_step)

                # checking that there is actually a variation and a range to test
                if(this_prop_range.size == 0):
                    # let's skip this method
                    this_prop_range = np.array([this_prop_min])
                # if that is too short/long, just generate a list of 50
                if(this_prop_range.size > 50):
                    this_prop_range = np.linspace(this_prop_min, this_prop_max, 50)

                # iterate between max and min for linewidth, and test the resulting data
                pbar = log.progressbar("adjusting rejection threshold for [" + properties_names[this_auto_method.value] + "]", this_prop_range.shape[0])

                test_snr_list = np.zeros(this_prop_range.shape)
                test_lw_list = np.zeros(this_prop_range.shape)
                test_nrej_list = np.zeros(this_prop_range.shape)

                # test each criteria bound
                for (i_prop_val, this_prop_val) in enumerate(this_prop_range):
                    # rebuild min/max rejection bounds
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
                    # by default, apply silently zerofilling, realigning, averaging and apodization
                    # TODO: maybe need to make this customizable?
                    if(this_mask_reject_data_sumup.sum() < s_ma.shape[0]):
                        old_level = log.getLevel()
                        log.setLevel(log.ERROR)
                        test_snr_list[i_prop_val], _, _ = this_s_cor.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd().analyze_snr_1d(peak_range)
                        test_lw_list[i_prop_val] = this_s_cor.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd().analyze_linewidth_1d(peak_range)
                        log.setLevel(old_level)
                        test_nrej_list[i_prop_val] = this_mask_reject_data_sumup.sum()

                    # progression
                    pbar.update(i_prop_val)

                pbar.finish("done")

                # analyze SNR curve
                test_snr_initial = test_snr_list[-1]
                test_snr_threshold = test_snr_initial + test_snr_initial * auto_adjust_allowed_snr_change / 100.0
                test_snr_list_rel = test_snr_list / test_snr_initial * 100.0 - 100.0

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
                        # that was a bit ambitious
                        log.info("sorry, this is only making your SNR worse...")
                        log.debug("the best SNR change we found was %.2f%% compared to initial :(" % test_snr_list_rel.max())
                        # set optimal LW to max
                        optim_prop = this_prop_max
                        optim_res_snr = test_snr_list[-1]
                        optim_res_lw = test_lw_list[-1]
                    else:
                        # we found SNR values which matches our request
                        # let's choose the one with the lowest LW
                        log.info("could not improve SNR above threshold but will reduce peak linewidth! :)")
                        ind_min_lw_snr_masked = np.argmin(test_lw_list[test_snr_list_mask])
                        optim_prop = this_prop_range[test_snr_list_mask][ind_min_lw_snr_masked]
                        optim_res_snr = test_snr_list[test_snr_list_mask][ind_min_lw_snr_masked]
                        optim_res_lw = test_lw_list[test_snr_list_mask][ind_min_lw_snr_masked]
                        log.info("optimal [" + properties_names[this_auto_method.value] + "] = %.1f" % optim_prop)

                # display and save the final snr and lw
                log.info("* Post-data-rejection based on [" + properties_names[this_auto_method.value] + "] SNR = %.2f" % optim_res_snr)
                log.info("* Post-data-rejection based on [" + properties_names[this_auto_method.value] + "] linewidth = %.2f Hz" % optim_res_lw)
                auto_method_final_snr_list[this_auto_method.value] = optim_res_snr
                auto_method_final_lw_list[this_auto_method.value] = optim_res_lw

                # plot SNR / LW combinaisons and optimal choice
                if(display):
                    # plot the SNRs versus LWs
                    fig = plt.figure(129 + (iround_data_rej - 1) * 3 + 1)
                    if(not display_axes_ready[0]):
                        # we just created the figure, let's create the axes
                        fig.clf()
                        fig.canvas.set_window_title("mrs.reco.MRSData2.correct_analyze_and_reject_2d (auto 1/2)")
                        fig.suptitle("adjusting data rejection criteria for [%s] (round #%d)" % (self.display_label, iround_data_rej))
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

                    # plot the data rejection percentage
                    fig = plt.figure(129 + (iround_data_rej - 1) * 3 + 2)
                    if(not display_axes_ready[1]):
                        # we just created the figure, let's create the axes
                        fig.clf()
                        fig.canvas.set_window_title("mrs.reco.MRSData2.correct_analyze_and_reject_2d (auto 2/2)")
                        fig.suptitle("estimating data rejection rates for [%s] (round #%d)" % (self.display_label, iround_data_rej))
                        fig.subplots(2, 2)
                        for a in fig.axes:
                            a.twinx()
                        display_axes_ready[1] = True

                    fig.axes[this_auto_method.value].plot(this_prop_range, test_nrej_list, 'ko-', label='Number of scans rejected')
                    fig.axes[this_auto_method.value].set_xlabel(properties_names[this_auto_method.value][0].upper() + properties_names[this_auto_method.value][1:])
                    fig.axes[this_auto_method.value].set_ylabel('Number of scans rejected')
                    fig.axes[this_auto_method.value].grid('on')

                    fig.axes[this_auto_method.value + 4].plot(this_prop_range, test_nrej_list / s_ma.shape[0] * 100, 'ko-', label='Total percentage of scans rejected')
                    fig.axes[this_auto_method.value + 4].set_ylabel('Rejection percentage (%)')

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
            auto_method_final_snr_list_rel = auto_method_final_snr_list / test_snr_initial * 100.0 - 100.0

            # is this higher than the initial snr?
            if(auto_method_final_snr_list_rel.max() > auto_adjust_allowed_snr_change):
                # apply this method
                ind_max_snr_auto_method = np.argmax(auto_method_final_snr_list)
                optim_auto_method = data_rejection_method(ind_max_snr_auto_method)
                log.info("best adjustment done with " + str(optim_auto_method) + " regarding SNR! (round #%d)" % iround_data_rej)
            else:
                # check that we have a segment of the curve above the initial SNR
                auto_method_final_snr_list_mask = (auto_method_final_snr_list_rel >= 0.0)
                if(not auto_method_final_snr_list_mask.any()):
                    # that was a bit ambitious
                    log.info("sorry but only made your SNR worse...")
                    log.debug("the best SNR change we found was %.2f%% compared to initial :(" % auto_method_final_snr_list_rel.max())
                    log.info("automatic data rejection failed, no optimal method found, sorry! :(")
                    # no optimal method!
                    optim_auto_method = None
                else:
                    # we found SNR values which matches our request
                    # let's choose the one with the lowest LW
                    log.info("could not find a method that improves SNR above threshold but will reduce the peak linewidth! :)")

                    # could not find a method giving a higher SNR than the initial SNR
                    # let's find a method that gives a lower linewidth ?
                    ind_min_lw_auto_method = np.argmin(auto_method_final_lw_list[auto_method_final_snr_list_mask])
                    optim_auto_method = data_rejection_method(ind_min_lw_auto_method)
                    log.info("best adjustment done with " + str(optim_auto_method) + " regarding linewidth! (round #%d)" % iround_data_rej)

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
        s_rej_avg = s_rej.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd()

        log.info("TOTAL data rejection = %d / %d (%.0f%%)" % (mask_reject_data_sumup.sum(), s_ma.shape[0], (mask_reject_data_sumup.sum() / s_ma.shape[0] * 100)))

        # perform post-correction measurements
        peak_prop_abs, peak_prop_rel2mean, peak_prop_rel2firstpt = s_ma._analyze_peak_2d(peak_range)

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

            fig = plt.figure(129 + (iround_data_rej - 1) * 3 + 3)
            fig.clf()
            axs = fig.subplots(2, 3, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_analyze_and_reject_2d")
            fig.suptitle("analyzing data and rejecting some for [%s] (round #%d)" % (self.display_label, iround_data_rej))

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
            plt.subplot(1, 3, 3)
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
                _, _, _, _, peak_seg_ppm, peak_seg_val = s_ma[k, :]._analyze_peak_1d(peak_range)
                plt.plot(peak_seg_ppm, np.real(peak_seg_val) * ampfactor + ystep * k, 'k-', linewidth=1)

            plt.xlim(peak_range[1], peak_range[0])
            plt.xlabel('chemical shift (ppm)')
            plt.ylabel('individual spectra')
            plt.yticks([])
            plt.grid('on')
            fig.subplots_adjust()
            fig.show()

        # wait, are we removing all data ???
        if(mask_reject_data_sumup.sum() == s.shape[0]):
            log.error("all data is rejected! You need to readjust your rejection bounds...")

        # display corected and rejected spectra
        if(display):
            s_avg = s.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd()
            s_cor_avg = s_cor.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd()

            # change legends
            s_avg.set_display_label("original spectrum")
            s_cor_avg.set_display_label("corrected spectrum")
            s_rej_avg.set_display_label("rejected spectrum")

            fig = plt.figure(129 + (iround_data_rej - 1) * 3 + 4)
            fig.clf()
            ax = fig.subplots()
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_analyze_and_reject_2d")
            fig.suptitle("original, corrected and rejected spectra for [%s] (round #%d)" % (self.display_label, iround_data_rej))

            ax.plot(s_rej_avg.frequency_axis_ppm(), s_rej_avg.spectrum().real, 'r-', linewidth=1, label=s_rej_avg.display_label)
            ax.plot(s_avg.frequency_axis_ppm(), s_avg.spectrum().real, 'k-', linewidth=1, label=s_avg.display_label)
            ax.plot(s_cor_avg.frequency_axis_ppm(), s_cor_avg.spectrum().real, 'b-', linewidth=1, label=s_cor_avg.display_label)

            if any(display_range):
                ax.set_xlim(display_range[1], display_range[0])

            ax.set_xlabel('chemical shift (ppm)')
            ax.set_ylabel('spectrum')
            ax.grid('on')
            ax.legend()

            fig.subplots_adjust()
            fig.show()

        # estimate final SNR and linewidth
        old_level = log.getLevel()
        log.setLevel(log.ERROR)
        final_snr, _, _ = s_cor.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd().analyze_snr_1d(peak_range)
        final_lw = s_cor.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_nd().analyze_linewidth_1d(peak_range)
        log.setLevel(old_level)
        log.info("* Final post-data-rejection SNR = %.2f" % final_snr)
        log.info("* Final post-data-rejection linewidth = %.2f Hz" % final_lw)

        # fill up dict about this data rejection
        data_rej_dict = {}
        data_rej_dict["Pre-rejection"] = {}
        data_rej_dict["Pre-rejection"]["snr"] = initial_snr
        data_rej_dict["Pre-rejection"]["lw"] = initial_lw
        data_rej_dict["Pre-rejection"]["na"] = s.shape[0]
        data_rej_dict["Pre-rejection"]["useful scantime"] = s.shape[0] * s._tr
        data_rej_dict["Pre-rejection"]["measurements"] = peak_prop_analyze
        data_rej_dict["Post-rejection"] = {}
        data_rej_dict["Post-rejection"]["snr"] = final_snr
        data_rej_dict["Post-rejection"]["lw"] = final_lw
        data_rej_dict["Post-rejection"]["na"] = s_cor.shape[0]
        data_rej_dict["Post-rejection"]["measurements"] = peak_prop_analyze_postcor
        data_rej_dict["Pre-rejection"]["useful scantime"] = s_cor.shape[0] * s_cor._tr
        # rejected spectrum
        data_rej_dict["Rejected spectrum"] = s_rej_avg
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

    def correct_realign_2d(self, peak_range=[4.5, 5], moving_Naverages=1, inter_corr_mode=False, display=False, display_range=[1, 6]):
        """
        Realign each signal of interest in frequency by taking as a reference the first spectra in absolute mode.

        * Works only with a 2D [averages,timepoints] signal.
        * Returns a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to analyze peak phase
        moving_Naverages : int
            Number of averages to perform when using moving average, need to be an odd number
        inter_corr_mode : boolean
            Use inter-correlation technique to adjust frequency shifts. Could be more robust when SNR is low.
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
            s_ma = self._build_moving_average_data_2d(moving_Naverages)

            # init
            s_avg = np.mean(s, axis=0)
            if(inter_corr_mode):
                # let's fix a +/-0.5ppm range
                f_shifts_min = - np.abs(peak_range[1] - peak_range[0]) * s_ma.f0
                f_shifts_max = + np.abs(peak_range[1] - peak_range[0]) * s_ma.f0
                # let's fix a 0.1ppm resolution here for inter-correlation tests
                f_shifts_step = 0.1 * s_ma.f0
                f_shifts_list = np.arange(f_shifts_min, f_shifts_max, f_shifts_step)
            else:
                # find peak in average spectrum absolute mode
                ippm_peak_avg, ppm_peak_avg, peak_val, _, _, _ = s_avg._analyze_peak_1d(peak_range)
                log.debug("measuring peak properties at %0.2fppm!" % ppm_peak_avg)

            # for each average in moving averaged data
            s_realigned_ma = s_ma.copy()
            df_trace = np.zeros(s_ma.shape[0])
            pbar = log.progressbar("realigning", s_ma.shape[0])
            for a in range(0, s_ma.shape[0]):

                if(inter_corr_mode):
                    # compare this individual spectrum with the first,  using inter-correlation

                    # use the peak_range as a range for inter-corr tests
                    cc_2d = f_shifts_list * 0.0
                    for ifs, fs in enumerate(f_shifts_list):
                        s_ma_ref = np.abs(s_ma[0, :].spectrum())
                        s_ma_shifted = np.abs(s_ma[a, :].adjust_frequency(fs).spectrum())
                        cc = np.corrcoef(s_ma_ref, s_ma_shifted)
                        cc_2d[ifs] = np.abs(cc[0, 1])

                    # find max correlation
                    optimal_fs_ind = np.argmax(cc_2d)
                    optimal_fs = f_shifts_list[optimal_fs_ind]
                    df_trace[a] = optimal_fs
                else:
                    # measure shift on moving average data
                    ippm_peak, ppm_peak, _, _, _, _ = s_ma[a, :]._analyze_peak_1d(peak_range)

                    # estimate frequency shift in Hz compared to average spectrum
                    dppm = -(ppm_peak_avg - ppm_peak)
                    df_trace[a] = dppm * s_ma.f0

                # correct moving averaged data
                s_realigned_ma[a, :] = s_ma[a, :].adjust_frequency(df_trace[a])

                # correct original data
                s_realigned[a, :] = s[a, :].adjust_frequency(df_trace[a])

                pbar.update(a)

            pbar.finish("done")

            # final display
            if(display):

                fig = plt.figure(140)
                fig.clf()
                axs = fig.subplots(2, 3, sharex='all', sharey='all')
                fig.canvas.set_window_title("mrs.reco.MRSData2.correct_realign_2d")
                fig.suptitle("frequency realigning [%s]" % self.display_label)

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

                fig = plt.figure(150)
                fig.clf()
                axs = fig.subplots(2, 1, sharex='all', sharey='all')
                fig.canvas.set_window_title("mrs.reco.MRSData2.correct_average_2d")
                fig.suptitle("averaging [%s]" % self.display_label)

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

    def correct_phase_1d(self, suspect_method=suspect_phasing_method.MATCH_MAGNITUDE_REAL, ppm_range=[0, 6], display=False, display_range=[1, 6]):
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

        # estimate phases
        if(suspect_method == suspect_phasing_method.MATCH_MAGNITUDE_REAL):
            phi0, phi1 = suspect.processing.phase.mag_real(s, range_ppm=ppm_range)
        elif(suspect_method == suspect_phasing_method.MIN_IMAG_INTEGRAL):
            phi0, phi1 = suspect.processing.phase.ernst(s)
        elif(suspect_method == suspect_phasing_method.ACME):
            phi0, phi1 = suspect.processing.phase.acme(s, range_ppm=ppm_range)
        else:
            log.error("hey, I do not know this suspect phasing method!?")

        # apply phase corrections
        s_phased = s.adjust_phase(phi0, phi1)

        # convert back to MRSData2
        s_phased = self.inherit(s_phased)

        if(display):
            fig = plt.figure(160)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_phase_1d")
            fig.suptitle("phasing (suspect) [%s]" % self.display_label)

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

            fig = plt.figure(170)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='row', sharey='row')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_apodization")
            fig.suptitle("apodizing [%s]" % self.display_label)

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
            if(s_crop.sequence is not None):
                log.debug("updating sequence.npts...")
                s_crop.sequence.npts = nPoints_final
                s_crop.sequence._ready = False
        else:
            s_crop = self.copy()
            log.debug("no cropping needed, getting bored...")

        if(display):
            t = s.time_axis()
            t_crop = s_crop.time_axis()
            ppm = s.frequency_axis_ppm()
            ppm_crop = s_crop.frequency_axis_ppm()

            fig = plt.figure(180)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='row', sharey='row')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_crop_1d")
            fig.suptitle("cropping [%s]" % self.display_label)

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

        # if any ref data available, we crop it too (silently)
        if(s_crop.data_ref is not None):
            s_crop.data_ref = s_crop.data_ref.correct_crop_1d(nPoints_final, False)

        return(s_crop)

    def correct_water_removal_1d(self, hsvd_nComponents=5, hsvd_range=[4.6, 4.8], display=False, display_range=[1, 6]):
        """
        Remove any water residual peak(s) within a ppm range using HSVD.

        * Works only with a 1D [timepoints] signal.
        * Returns a 1D [timepoints] signal.

        Parameters
        ----------
        hsvd_nComponents : int
            Number of components for HSVD water residue removal
        hsvd_range : list [2]
            Range in ppm of HSVD components to keep for the water peak removal
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_water_removed : MRSData2 numpy array [timepoints]
            Resulting water HSVD suppressed data stored in a MRSData2 object
        """
        log.debug("removing water peak for [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()
        pbar = log.progressbar("removing residual water peak with HSVD", 5)

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
        s_water_removed = s - hsvd_fid
        pbar.update(5)

        # display this over the data
        if(display):
            fig = plt.figure(190)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_water_removal_1d")
            fig.suptitle("removing water peak for [%s]" % self.display_label)

            # original spectrum
            axs[0].plot(ppm, s.spectrum().real, 'k-', linewidth=1, label='original data (real part)')
            axs[0].plot(ppm, hsvd_fid.spectrum().real, 'r-', linewidth=1, label='estimated HSVD water peak')
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('original spectrum')
            axs[0].grid('on')
            axs[0].legend()

            # water removed spectrum
            axs[1].plot(ppm, s_water_removed.spectrum().real, 'b-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('water removed spectrum')
            axs[1].grid('on')
            axs[1].legend()

            fig.subplots_adjust()
            fig.show()

        pbar.finish("done")

        # convert back to MRSData2
        s_water_removed = self.inherit(s_water_removed)

        return(s_water_removed)

    def correct_freqshift_1d(self, peak_range=[4.5, 5], peak_real_ppm=4.7, display=False, display_range=[1, 6]):
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
        _, ppm_peak, peak_val, _, _, _ = s._analyze_peak_1d(peak_range)
        log.debug("peak detected at %0.2fppm -> %0.2fppm!" % (ppm_peak, peak_real_ppm))

        # estimate frequency shift in Hz
        log.debug("frequency shifting data...")
        dppm = (peak_real_ppm - ppm_peak)
        df = dppm * s.f0
        s_shifted = s.adjust_frequency(-df)

        if(display):
            fig = plt.figure(200)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_freqshift_1d")
            fig.suptitle("calibrating [%s]" % self.display_label)

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
            fig = plt.figure(210)
            fig.clf()
            axs = fig.subplots(2, 2)
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_bandpass_filtering_1d")
            fig.suptitle("filtering [%s]" % self.display_label)

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

    def analyze_snr_1d(self, peak_range, noise_range=[-1, 0], magnitude_mode=False, display=False, display_range=[1, 6]):
        """
        Estimate the SNR of a peak in the spectrum ; chemical shift ranges for the peak and the noise regions are specified by the user. Can also look at time-domain SNR. Works only for a 1D MRSData2 objects.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to find a peak of interest
        noise_range : list [2]
            Range in ppm used to estimate noise
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
            fig = plt.figure(220)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_snr_1d")
            fig.suptitle("analyzing SNR for [%s]" % self.display_label)

        # find maximum peak in range
        sf = s.spectrum()
        _, ppm_peak, peak_val, _, _, _ = s._analyze_peak_1d(peak_range)
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

    def analyze_linewidth_1d(self, POI_range_ppm, magnitude_mode=False, display=False, display_range=[1, 6]):
        """
        Estimate the linewidth of a peak in the spectrum ; chemical shift ranges for the peak and the noise regions are specified by the user.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        POI_range_ppm : list [2]
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

        # call 1D peak analysis
        _, ppm_peak, _, lw, peak_seg_ppm, peak_seg_val = s._analyze_peak_1d(POI_range_ppm)
        if(magnitude_mode):
            log.debug("estimating the MAGNITUDE peak linewidth at %0.2fppm!" % ppm_peak)
        else:
            log.debug("estimating the REAL peak linewidth at %0.2fppm!" % ppm_peak)

        log.info("results for [" + s.display_label + "] coming...")
        log.info("LW = %0.2f Hz!" % lw)

        if(display):
            ppm = s.frequency_axis_ppm()
            sf = s.spectrum()

            fig = plt.figure(230)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_linewidth_1d")
            fig.suptitle("analyzing peak linewidth for [%s]" % self.display_label)

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

    def display_spectrum_1d(self, ifig=1, display_range=[1, 6], magnitude_mode=False):
        """
        Display spectrum in figure 'ifig', overlaying if needed.

        * Works only with a 1D [timepoints] signal.

        Parameters
        ----------
        ifig: int
            The figure index that shoud host the plot
        s : MRSData2 numpy array [timepoints]
            MRS data to display
        display_range : list [2]
            Range in ppm used for display
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
        s = self.copy()
        log.debug("displaying stuff!")

        plt.figure(ifig).canvas.set_window_title("mrs.reco.MRSData2.display_spectrum_1d")
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

        return(plt)

    def save_ismrmd(self, h5_filepath):
        """
        Save the MRSData2 object to a ISMRMRD format file. This function depends on ismrmrd-python, available at https://github.com/ismrmrd/ismrmrd-python. General info about this open file format is available here: https://ismrmrd.github.io/. Got inspired by this example: https://github.com/ismrmrd/ismrmrd-python-tools/blob/master/generate_cartesian_shepp_logan_dataset.py .

        Parameters
        ----------
        h5_filepath: string
            Full absolute file path pointing to h5 file
        """
        log.debug("saving MRS signal to " + h5_filepath + "...")
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 3D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # try importing here in order not to break everything because you do not have this dependency
        import ismrmrd

        # TODO: need to add more info to the header: voxel size/position, etc.

        # open  dataset
        dset = ismrmrd.Dataset(h5_filepath, "dataset", create_if_needed=True)

        # create the XML header
        header = ismrmrd.xsd.ismrmrdHeader()

        # subject stuff
        subj = ismrmrd.xsd.subjectInformationType()
        subj.patientName = self.patient["name"]
        subj.patientID = self.patient["name"]
        # bug here...
        # subj.patientBirthdate = self.patient["birthday"].strftime('%Y-%m-%d')
        if(self.patient["sex"] == 0):
            subj.patientGender = "M"
        elif(self.patient["sex"] == 1):
            subj.patientGender = "F"
        elif(self.patient["sex"] == 2):
            subj.patientGender = "O"
        else:
            log.error("patient gender unknown!")
        subj.patientWeight_kg = self.patient["weight"]
        # no patient height in the ismrmrd format!?
        # add to header
        header.subjectInformation = subj

        # sequence stuff
        seq = ismrmrd.xsd.sequenceParametersType()
        seq.TR = self.sequence.tr
        seq.TE = [self.sequence.te]
        seq.sequence_type = self.sequence.name
        # add to header
        header.sequenceParameters = seq

        # experimental conditions
        exp = ismrmrd.xsd.experimentalConditionsType()
        exp.H1resonanceFrequency_Hz = self.f0 * 1e6
        header.experimentalConditions = exp

        # dummy encoding
        encoding = ismrmrd.xsd.encoding()
        encoding.trajectory = ismrmrd.xsd.trajectoryType.cartesian
        efov = ismrmrd.xsd.fieldOfView_mm()
        efov.x = 1
        efov.y = 1
        efov.z = 1
        rfov = ismrmrd.xsd.fieldOfView_mm()
        rfov.x = 1
        rfov.y = 1
        rfov.z = 1
        ematrix = ismrmrd.xsd.matrixSize()
        ematrix.x = 1
        ematrix.y = 1
        ematrix.z = 1
        rmatrix = ismrmrd.xsd.matrixSize()
        rmatrix.x = 1
        rmatrix.y = 1
        rmatrix.z = 1
        espace = ismrmrd.xsd.encodingSpaceType()
        espace.matrixSize = ematrix
        espace.fieldOfView_mm = efov
        rspace = ismrmrd.xsd.encodingSpaceType()
        rspace.matrixSize = rmatrix
        rspace.fieldOfView_mm = rfov
        encoding.encodedSpace = espace
        encoding.reconSpace = rspace
        limits = ismrmrd.xsd.encodingLimitsType()
        encoding.encodingLimits = limits
        header.encoding.append(encoding)

        # write header
        dset.write_xml_header(header.toxml('utf-8'))

        # data
        acq = ismrmrd.Acquisition()
        acq.resize(self.shape[0], 0)
        acq.data[:] = self[:]
        dset.append_acquisition(acq)

        # done
        dset.close()

    def save_mat(self, mat_filepath):
        """
        Save the numpy array content to a MATLAB mat file.

        Parameters
        ----------
        mat_filepath: string
            Full absolute file path pointing to mat file
        """
        # TODO: need to add some header info
        log.debug("saving MRS signal to " + mat_filepath + "...")
        sio.savemat(mat_filepath, {'MRSdata': self})

    def save_pkl(self, pkl_filepath):
        """
        Save the whole object to a pickle file.

        Parameters
        ----------
        pkl_filepath: string
            Full absolute file path pointing to pkl file
        """
        log.debug("saving MRS signal to " + pkl_filepath + "...")
        with open(pkl_filepath, 'wb') as f:
            pickle.dump(self, f)

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
        self.dataset = [{}] * max_number_datasets_pipeline
        for i in range(max_number_datasets_pipeline):
            self.dataset[i] = {"legend": None,
                               "raw": {"files": [None, None], "data": None, "analysis-results": None, "ref-data-analysis-results": None},
                               "dcm": {"files": [None, None], "data": None, "analysis-results": None, "ref-data-analysis-results": None},
                               "physio-file": None,
                               "imaging-file": None}

        # --- global settings ---
        self.settings = {   # option to process only a set of datasets: list of indexes
                            "datasets_indexes": None,
                            # ppm scale reference
                            "ppm0": 4.7,
                            # ppm range to search for peak used for phasing, etc.
                            "POI_range_ppm": [4.5, 5.2],
                            # ppm range to search for ppm scale calibration
                            "POI_shift_range_ppm": [4.5, 5.2],
                            # real ppm value the above peak
                            "POI_shift_true_ppm": 4.7,
                            # ppm range to search for peak for SNR estimation
                            "POI_SNR_range_ppm": [1.8, 2.1],
                            # ppm range to search for peak for FWHM estimation
                            "POI_LW_range_ppm": [4.5, 5.2],
                            # ppm range used for display
                            "display_range_ppm": [1, 6],
                            # y offset used for display
                            "display_offset": 0.0}

        # --- available jobs and their parameters ---
        self.job = {}
        # --- job: spectrum final display ---
        self.job["displaying"] = {0:
                                  {"func": MRSData2.display_spectrum_1d, "name": "displaying"},
                                  # figure index
                                  "fig_index": 1,
                                  # ppm range used for display
                                  "range_ppm": self.settings["display_range_ppm"],
                                  # display spectrum in magnitude mode?
                                  "magnitude_mode": False
                                  }

        # --- job: automatic rephasing ---
        self.job["phasing"] = {0:
                               {"func": MRSData2.correct_phase_3d, "name": "phasing"},
                               # use reference data is available?
                               "using_ref_data": True,
                               # ppm range to look fo peak used to estimate phase
                               "POI_range_ppm": self.settings["POI_range_ppm"],
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
                               "display_range_ppm": self.job["displaying"]["range_ppm"]
                               }

        # --- job: amplification ---
        self.job["scaling"] = {0:
                               {"func": MRSData2.correct_intensity_scaling_nd, "name": "scaling intensity"},
                               "scaling_factor_rawdata": 1e8,
                               "scaling_factor_dcm": 1.0
                               }

        # --- job: FID modulus ---
        self.job["FID modulus"] = {0:
                                   {"func": MRSData2.correct_fidmodulus_nd, "name": "FID modulus"}
                                   }

        # --- job: channel combination ---
        self.job["channel-combining"] = {0:
                                         {"func": MRSData2.correct_combine_channels_3d, "name": "channel-combining"},
                                         # use non water-suppressed data to recombine and rephase channels
                                         "using_ref_data": True,
                                         # should we rephase (0th order) data while combining?
                                         "phasing": False,
                                         # boolean mask to switch on/off some Rx channels
                                         "weights": [True]
                                         }

        # --- job: concatenate ---
        self.job["concatenate"] = {0:
                                   {"func": MRSData2.concatenate_2d, "name": "concatenate"}
                                   }

        # --- job: zero-filling ---
        self.job["zero-filling"] = {0:
                                    {"func": MRSData2.correct_zerofill_nd, "name": "zero-filling"},
                                    # number of signal points after zf
                                    "npts": 8192 * 2,
                                    # display all this process to check what the hell is going on
                                    "display": True,
                                    "display_range_ppm": self.job["displaying"]["range_ppm"]
                                    }

        # --- job: analyze physio signal ---
        self.job["physio-analysis"] = {0:
                                       {"func": MRSData2.analyze_physio_2d, "name": "analyzing physio. signals"},
                                       # ppm range to look for a peak to analyze
                                       "POI_range_ppm": self.settings["POI_range_ppm"],
                                       # time range in (ms) to look around timestamp for correlation physio/MRS
                                       "delta_time_ms": 1000.0,
                                       # display all this process to check what the hell is going on
                                       "display": True
                                       }

        # --- job: automatic data rejection based on criterias ---
        self.job["data-rejecting"] = {0:
                                      {"func": MRSData2.correct_analyze_and_reject_2d, "name": "data rejecting"},
                                      # ppm range to look for a peak to analyze
                                      "POI_range_ppm": self.settings["POI_range_ppm"],
                                      # size of moving average window
                                      "moving_averages": 1,
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
                                      "auto_method_list": None,
                                      # minimum allowed SNR change (%) when adjusting the linewidth criteria, this can be positive (we want to increase SNR +10% by rejecting crappy dat) or negative (we are ok in decreasing the SNR -10% in order to get better resolved spectra)
                                      "auto_allowed_snr_change": 1.0,
                                      # display all this process to check what the hell is going on
                                      "display": True
                                      }

        # --- job: automatic data frequency realignment ---
        self.job["realigning"] = {0:
                                  {"func": MRSData2.correct_realign_2d, "name": "frequency realigning"},
                                  # ppm range to look for a peak to analyze
                                  "POI_range_ppm": self.settings["POI_range_ppm"],
                                  # size of moving average window
                                  "moving_averages": 1,
                                  # use correlation mode
                                  "inter_corr_mode": False,
                                  # display all this process to check what the hell is going on
                                  "display": True,
                                  "display_range_ppm": self.job["displaying"]["range_ppm"]
                                  }

        # --- job: spectral filtering ---
        self.job["filtering"] = {0:
                                 {"func": MRSData2.correct_bandpass_filtering_1d, "name": "FFT filtering"},
                                 # ppm range to keep
                                 "range_ppm": [1, 6],
                                 # type of apodization window (take it from numpy/scipy)
                                 "window_func": signal.tukey,
                                 # display all this process to check what the hell is going on
                                 "display": True,
                                 "display_range_ppm": self.job["displaying"]["range_ppm"]
                                 }

        # --- job: data averaging ---
        self.job["averaging"] = {0:
                                 {"func": MRSData2.correct_average_2d, "name": "averaging"},
                                 # number of averages to mean (None = all)
                                 "na": None,
                                 # display all this process to check what the hell is going on
                                 "display": True,
                                 "display_range_ppm": self.job["displaying"]["range_ppm"]
                                 }

        # --- job: phasing using suspect ---
        self.job["phasing (suspect)"] = {0:
                                         {"func": MRSData2.correct_phase_1d, "name": "phasing (suspect)"},
                                         # phasing method
                                         "suspect_method": suspect_phasing_method.MATCH_MAGNITUDE_REAL,
                                         # ppm range to analyze phase
                                         "range_ppm": [1, 6],
                                         # display all this process to check what the hell is going on
                                         "display": True,
                                         "display_range_ppm": self.job["displaying"]["range_ppm"]
                                         }

        # --- job: noise level analysis ---
        self.job["noise-estimation"] = {0:
                                        {"func": MRSData2.analyze_noise_nd, "name": "estimating noise level"},
                                        # estimate noise std time-domain on the last 100 pts of the FID
                                        "npts": 100
                                        }

        # --- job: data apodization ---
        self.job["apodizing"] = {0:
                                 {"func": MRSData2.correct_apodization_nd, "name": "apodizing"},
                                 # exponential damping factor for apodization (Hz)
                                 "damping_hz": 5,
                                 # display all this process to check what the hell is going on
                                 "display": True,
                                 "display_range_ppm": self.job["displaying"]["range_ppm"]
                                 }

        # --- job: data cropping ---
        self.job["cropping"] = {0:
                                {"func": MRSData2.correct_crop_1d, "name": "cropping"},
                                # final number of signal points after crop
                                "final_npts": 6144,
                                # display all this process to check what the hell is going on
                                "display": True,
                                "display_range_ppm": self.job["displaying"]["range_ppm"]
                                }

        # --- job: water post-acquisition removal ---
        self.job["water-removal"] = {0:
                                     {"func": MRSData2.correct_water_removal_1d, "name": "removing water peak"},
                                     # number of components when running HSVD
                                     "hsvd_components": 5,
                                     # ppm range where all components will be remove
                                     "hsvd_range": self.settings["POI_range_ppm"],
                                     # display all this process to check what the hell is going on
                                     "display": True,
                                     "display_range_ppm": self.job["displaying"]["range_ppm"]
                                     }

        # --- job: spectrum chemical shift calibration ---
        self.job["calibrating"] = {0:
                                   {"func": MRSData2.correct_freqshift_1d, "name": "frequency shifting"},
                                   # ppm range to look for the peak of interest (NAA by default)
                                   "POI_shift_range_ppm": self.settings["POI_shift_range_ppm"],
                                   # real ppm value for this peak
                                   "POI_shift_true_ppm": self.settings["POI_shift_true_ppm"],
                                   # display all this process to check what the hell is going on
                                   "display": True,
                                   "display_range_ppm": self.job["displaying"]["range_ppm"]
                                   }

        # --- job: SNR analysis ---
        self.job["analyzing-snr"] = {0:
                                     {"func": MRSData2.analyze_snr_1d, "name": "analyzing SNR"},
                                     # ppm range to look for a peak to analyze
                                     "POI_SNR_range_ppm": self.settings["POI_SNR_range_ppm"],
                                     # ppm range to look for pure noise
                                     "n_range_ppm": [-1, 0],
                                     # should we look at the magnitude or real spectrum?
                                     "magnitude_mode": False,
                                     # display all this process to check what the hell is going on
                                     "display": True,
                                     "display_range_ppm": self.job["displaying"]["range_ppm"]
                                     }

        # --- job: LW analysis ---
        self.job["analyzing-lw"] = {0:
                                    {"func": MRSData2.analyze_linewidth_1d, "name": "analyzing peak-linewidth"},
                                    # ppm range to look for a peak to analyze
                                    "POI_range_ppm": self.settings["POI_range_ppm"],
                                    # should we look at the magnitude or real spectrum?
                                    "magnitude_mode": False,
                                    # display all this process to check what the hell is going on
                                    "display": True,
                                    "display_range_ppm": self.job["displaying"]["range_ppm"]
                                    }

        # --- job list ---
        # list of data processing to apply to the data
        # beware, you need to know what you are doing here
        # also, be careful with the dimensionality of data 3D, 2D, 1D along the data processing
        # order is important!
        self.job_list = [self.job["phasing"],
                         self.job["scaling"],
                         self.job["FID modulus"],
                         self.job["channel-combining"],
                         self.job["concatenate"],
                         self.job["zero-filling"],
                         self.job["physio-analysis"],
                         self.job["data-rejecting"],
                         self.job["realigning"],
                         self.job["averaging"],
                         self.job["noise-estimation"],
                         self.job["apodizing"],
                         self.job["cropping"],
                         self.job["water-removal"],
                         self.job["calibrating"],
                         self.job["displaying"]]

        # --- analyze job list ---
        # SNR/LW analysis job list
        self.analyze_job_list = [self.job["channel-combining"],
                                 self.job["zero-filling"],
                                 self.job["averaging"],
                                 self.job["calibrating"]]

        # --- SNR/LW analysis ---
        self.analyze_enable = True
        self.analyze_display = True

        # --- template loading if needed ---
        if(template_name is not None):
            # overwrite everything above with the template
            self.load_template(template_name)

        # freeze the object and prevent the creation of new attributes
        self.__isfrozen = True

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
        log.info_line_break()
        log.info_line________________________()
        job_name = job[0]["name"]
        log.info("%s on [%s]..." % (job_name, data.display_label))
        # get function
        job_func = job[0]["func"]
        if(default_args):
            job_args = [data]
        else:
            # get arguments
            job_args = job.copy()
            del job_args[0]
            job_args = [data] + list(job_args.values())

        # call job on data
        job_result = job_func(*job_args)
        log.info_line________________________()

        # return
        return(job_result)

    def _analyze(self, data, current_job, already_done_jobs):
        """
        Estimate SNR and/or peak linewidth for this dataset. Values are stored. A mini default pipeline is applied before SNR/LW measurements and can be set with self.analyze_job_list.

        Parameters
        ----------
        data : MRSData2 object [whatever,...,timepoints]
            One of the processing function of the MRSData2 class
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
        """
        log.debug("estimating SNR and peak-linewidth for [%s]..." % data.display_label)

        # init job list
        this_analyze_job_list = [j for j in self.analyze_job_list + already_done_jobs if (j in self.analyze_job_list) and (j not in already_done_jobs)]

        # run mini-pipeline with default arguments (= no display)
        # and with no log output
        old_level = log.getLevel()
        log.setLevel(log.ERROR)
        for j in this_analyze_job_list:
            # run job on this dataset with default arguments
            data = self._run_job(j, data, True)

        # measure snr
        data_snr, _, _ = self._run_job(self.job["analyzing-snr"], data)
        data_lw = self._run_job(self.job["analyzing-lw"], data)

        # allow outputs
        log.setLevel(old_level)

        # output
        job_label = "post-" + current_job[0]["name"]
        log.info(job_label + " SNR of [%s] = %.2f" % (data.display_label, data_snr))
        log.info(job_label + " LW of [%s] = %.2f" % (data.display_label, data_lw))

        return(data_snr, data_lw)

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
        Run the pipeline for MRS reco. It includes.

            * reading the data
            * phasing
            * combining channels
            * analyze real noise level before going on with processing
            * zero-filling
            * peak analysis & data rejection
            * realigning
            * averaging
            * checking original noise level
            * apodizing
            * cropping
            * removing water reidue with HSVD
            * shifting the spectrum
            * analyzing the SNR
            * analyzing the linewidth
            * displaying the spectrum
            +
            * analyzing the SNR evolution during data processing
            * analyzing the linewidth evolution during data processing

        Returns
        -------
        self._data_list : list of MRSData2 objects
            Final MRS data signals obtained from reconstruction pipeline stored in a MRSData2 object
        """
        # --- reading dataset dict
        log.info_line________________________()
        log.info("checking datasets...")
        log.info_line________________________()

        # for each data set, check dict fields
        ind_dataset_ok = []
        for i, d in enumerate(self.dataset):
            legend_ok = (d["legend"] is not None)
            rawfiles_ok = (d["raw"]["files"][0] is not None)
            dicomfiles_ok = (d["dcm"]["files"][0] is not None)
            physiofile_ok = (d["physio-file"] is not None)
            imagingfile_ok = (d["imaging-file"] is not None)

            # checking that this is an entry
            if(legend_ok or rawfiles_ok or dicomfiles_ok or physiofile_ok or imagingfile_ok):
                # seems like a entry

                # checking legend
                if(not legend_ok):
                    log.error("dataset[%d] is missing a legend! :(" % i)
                # checking datafiles
                if(not rawfiles_ok and not dicomfiles_ok):
                    log.error("dataset[%d] is missing data, no raw or dicom files were set! :(" % i)
                # checking number of datafiles
                if(rawfiles_ok and len(d["raw"]["files"]) > 2):
                    log.error("dataset[%d] has %d raw-files! There should be 1 or 2. :(" % (i, len(d["raw"]["files"])))
                if(dicomfiles_ok and len(d["dcm"]["files"]) > 2):
                    log.error("dataset[%d] has %d dicom-files! There should be 1 or 2. :(" % (i, len(d["dcm"]["files"])))

                # clean and print strings
                log.info("dataset[%d]:" % i)
                # legend
                if(legend_ok):
                    self.dataset[i]["legend"] = d["legend"].strip()
                log.info("  legend = " + self.dataset[i]["legend"])
                # raw data
                if(self.dataset[i]["raw"]["files"][0] is not None):
                    self.dataset[i]["raw"]["files"][0] = self.dataset[i]["raw"]["files"][0].strip()
                    log.info("  raw #0 = " + self.dataset[i]["raw"]["files"][0])
                if(len(self.dataset[i]["raw"]["files"]) == 1):
                    self.dataset[i]["raw"]["files"].append(None)
                if(self.dataset[i]["raw"]["files"][1] is not None):
                    self.dataset[i]["raw"]["files"][1] = self.dataset[i]["raw"]["files"][1].strip()
                    log.info("  raw #1 [REF] = " + self.dataset[i]["raw"]["files"][1])
                else:
                    log.info("  raw #1 [REF] = None")
                # dicom
                if(self.dataset[i]["dcm"]["files"][0] is not None):
                    self.dataset[i]["dcm"]["files"][0] = self.dataset[i]["dcm"]["files"][0].strip()
                    log.info("  dicom #0 = " + self.dataset[i]["dcm"]["files"][0])
                if(len(self.dataset[i]["dcm"]["files"]) == 1):
                    self.dataset[i]["dcm"]["files"].append(None)
                if(self.dataset[i]["dcm"]["files"][1] is not None):
                    self.dataset[i]["dcm"]["files"][1] = self.dataset[i]["dcm"]["files"][1].strip()
                    log.info("  dicom #1 [REF] = " + self.dataset[i]["dcm"]["files"][1])
                else:
                    log.info("  dicom #1 [REF] = None")
                # physio
                if(physiofile_ok):
                    self.dataset[i]["physio-file"] = d["physio-file"].strip()
                    log.info("  physio = " + self.dataset[i]["physio-file"])
                # imaging
                if(imagingfile_ok):
                    self.dataset[i]["imaging-file"] = d["imaging-file"].strip()
                    log.info("  imaging = " + self.dataset[i]["imaging-file"])

                # index of non-empty datasets
                ind_dataset_ok.append(i)

        # filter datasets
        if(self.settings["datasets_indexes"] is not None):
            # if only one index, convert to list
            if(type(self.settings["datasets_indexes"]) == int):
                self.settings["datasets_indexes"] = [self.settings["datasets_indexes"]]

            # keep only data to process
            self.dataset = [self.dataset[i] for i in self.settings["datasets_indexes"]]
        else:
            # keep only the non-empty datasets
            self.dataset = [self.dataset[i] for i in ind_dataset_ok]
        log.info_line________________________()

        # --- reading data ---
        # now let's read the data files
        log.info_line_break()
        log.info_line________________________()
        log.info("reading data files...")
        log.info_line________________________()
        for i, d in enumerate(self.dataset):
            this_legend = d["legend"]
            this_physio_filename = d["physio-file"]

            for dtype in ["raw", "dcm"]:
                # check if any data file to read
                this_data_filename = d[dtype]["files"][0]
                if(this_data_filename is None):
                    # no? so pass
                    continue

                this_data_ref_filename = d[dtype]["files"][1]
                log.info("reading data [" + this_legend + "]")
                s = MRSData2(this_data_filename, this_physio_filename)
                log.debug("got a " + str(s.shape) + " vector")
                # set ppm reference
                s.ppm0 = self.settings["ppm0"]
                # set legend & offset
                this_new_legend = ("#%d " % i) + this_legend
                if(s.is_rawdata):
                    s.set_display_label(this_new_legend + " [RAW]")
                else:
                    s.set_display_label(this_new_legend + " [DCM]")

                s.set_display_offset(self.settings["display_offset"])
                # store
                self.dataset[i][dtype]["data"] = s
                self.dataset[i]["legend"] = this_new_legend

                if(this_data_ref_filename is not None):
                    log.info_line_break()
                    log.info("reading ref. data [" + this_legend + "]")
                    s = MRSData2(this_data_ref_filename, None)
                    log.info("got a " + str(s.shape) + " vector")
                    # if several averages, mean now (that could be a problem!?)
                    s = s.mean(axis=0)
                    # add 1 dimension
                    s = s.reshape((1,) + s.shape)
                    log.debug("reshaped to a " + str(s.shape) + " vector")
                    # set ppm reference
                    s.ppm0 = self.settings["ppm0"]
                    # set legend & offset
                    this_new_legend = ("#%d " % i) + this_legend
                    if(s.is_rawdata):
                        s.set_display_label(this_new_legend + " [REF] [RAW]")
                    else:
                        s.set_display_label(this_new_legend + " [REF] [DCM]")
                    s.set_display_offset(self.settings["display_offset"])
                    # store
                    self.dataset[i][dtype]["data"].data_ref = s

                log.info_line________________________()

        # --- reading job list ---
        log.info_line________________________()
        log.info("reading your job list...")
        log.info_line________________________()
        # applying some global settings to jobs
        # (BTW, PYTHON CANNOT HANDLE POINTERS -_-)
        for this_setting in list(self.settings.keys()):
            for this_job in list(self.job.keys()):
                if(this_setting in self.job[this_job]):
                    self.job[this_job][this_setting] = self.settings[this_setting]

        # for each job
        for k, this_job in enumerate(self.job_list):
            this_job_name = this_job[0]["name"]
            log.info("#%d %s" % (k, this_job_name))
        log.info_line________________________()

        # --- running job list ---
        # exception here: if any concatenate, we need to run the job list in two parts (before and after concatenation in order to reload the processed data)
        log.info_line_break()
        log.info_line________________________()
        log.info("running job list...")
        log.info_line________________________()
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

                    log.info_line_break()
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
                        job_result = self._run_job(job, this_data)

                        # push job in the stack of jobs done
                        jobs_done_stack.append(job)

                        # replace with processed signal
                        if(type(job_result) == MRSData2):
                            this_data = job_result
                            # measure SNR/LW after this process?
                            if(self.analyze_enable):
                                # prepare storage
                                if(self.dataset[i][dtype]["analysis-results"] is None):
                                    self.dataset[i][dtype]["analysis-results"] = {}

                                # get job name, renaming if needed
                                this_job_name = job[0]["name"]
                                if(this_job_name in list(self.dataset[i][dtype]["analysis-results"].keys())):
                                    # not the first time we apply this job?
                                    this_job_name = this_job_name + " (#2)"

                                # get SNR and LW estimations
                                this_data_snr, this_data_lw = self._analyze(this_data, job, jobs_done_stack)

                                # store analysis results
                                self.dataset[i][dtype]["analysis-results"][this_job_name] = {"snr": this_data_snr, "lw": this_data_lw}
                                log.info_line________________________()

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

        # before leaving, analyse ref data if available
        for i, d in enumerate(self.dataset):
            for dtype in ["raw", "dcm"]:
                this_data = d[dtype]["data"]
                # if no data to process (= no DCM provided)
                if(this_data is None):
                    continue

                if(this_data.data_ref is not None):
                    this_data_ref_snr = self._run_job(self.job["analyzing-snr"], this_data.data_ref)
                    this_data_ref_lw = self._run_job(self.job["analyzing-lw"], this_data.data_ref)
                    self.dataset[i][dtype]["ref-data-analysis-results"] = {"snr": this_data_ref_snr, "lw": this_data_ref_lw}

        # --- summary final linewidths ---
        if(self.analyze_enable):
            self.display_analyze_results()

        log.info("pipeline terminated!")
        log.info("returning processed signals!")
        return(self.dataset)

    def get_analyze_results(self):
        """Return analyse results as several numpy vectors eassy to plot.Dataset and job names will be padded to ease terminal output.

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
        log.info("converting analyze results to arrays...")

        # build label list
        data_labels = [d["legend"] for d in self.dataset]
        data_label_nchar = len(max(data_labels, key=len))
        data_label_list = [d.ljust(data_label_nchar) for d in data_labels]

        # build job label list
        if(self.dataset[0]["raw"]["analysis-results"] is not None):
            job_labels = self.dataset[0]["raw"]["analysis-results"].keys()
        elif(self.dataset[0]["dcm"]["analysis-results"] is not None):
            job_labels = self.dataset[0]["dcm"]["analysis-results"].keys()
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
            if(d["raw"]["analysis-results"] is not None):
                snr_raw_list.append([d["raw"]["analysis-results"][j]["snr"] for j in d["raw"]["analysis-results"].keys()])
                lw_raw_list.append([d["raw"]["analysis-results"][j]["lw"] for j in d["raw"]["analysis-results"].keys()])
            else:
                snr_raw_list.append([np.nan] * len(job_label_list))
                lw_raw_list.append([np.nan] * len(job_label_list))

            # if dcm data, find snr and lw estimations and store it
            if(d["dcm"]["analysis-results"] is not None):
                snr_dcm_list.append([d["dcm"]["analysis-results"][j]["snr"] for j in d["dcm"]["analysis-results"].keys()])
                lw_dcm_list.append([d["dcm"]["analysis-results"][j]["lw"] for j in d["dcm"]["analysis-results"].keys()])
            else:
                snr_dcm_list.append([np.nan] * len(job_label_list))
                lw_dcm_list.append([np.nan] * len(job_label_list))

            # if raw ref data, find snr and lw estimations and store it
            if(d["raw"]["ref-data-analysis-results"] is not None):
                snr_ref_raw_list.append([d["raw"]["ref-data-analysis-results"]["snr"]])
                lw_ref_raw_list.append([d["raw"]["ref-data-analysis-results"]["lw"]])
            else:
                snr_ref_raw_list.append([np.nan])
                lw_ref_raw_list.append([np.nan])

            # if dcm ref data, find snr and lw estimations and store it
            if(d["dcm"]["ref-data-analysis-results"] is not None):
                snr_ref_dcm_list.append([d["dcm"]["ref-data-analysis-results"]["snr"]])
                lw_ref_dcm_list.append([d["dcm"]["ref-data-analysis-results"]["lw"]])
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
        """Return final analyse results for each dataset.

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
        log.info("returning final analyze results...")

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

    def display_analyze_results(self):
        """Print final SNR and peak-linewidth for each dataset. Plot bargraph showing evolution of SNR and linewidth during data processing (to check that a job did not destroy the data for example!) and compare with dicom when possible."""
        log.info("displaying SNR and linewidth final results...")

        # get analyse data as list and np arrays
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
        if(self.analyze_display):
            fig = plt.figure(300)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='row')
            fig.canvas.set_window_title("mrs.reco.pipeline.display_analyze_results")
            fig.suptitle("SNR/peak-linewidth results")

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

    def display_final_data(self):
        """Plot final processed datasets."""
        log.info("displaying final processed data...")

        # get display job
        disp_job = self.job["displaying"]

        # run a diplay job for each dataset
        for d in self.dataset:
            for dtype in ["raw", "dcm"]:
                this_data = d[dtype]["data"]
                # if no data
                if(this_data is None):
                    continue
                # run job on this dataset
                self._run_job(disp_job, this_data)

    def save_datasets(self, rdb):
        """
        Save each data set and the current pipeline to the PKL data storage file.

        Parameters
        ----------
        rdb: data_db object
            Data storage object
        """
        # for each dataset
        log.info("saving all datasets to file [%s]..." % rdb.db_file)
        for d in self.dataset:
            rdb.save_dataset(d, self)

    def load_template(self, template_name):
        """
        Load reco pipeline attributes from template previously saved.

        Parameters
        ----------
        template_name: string
            Template name
        """
        template_filename = template_name + ".pkl"
        log.info("load pipeline from template file [%s]..." % template_filename)
        # check if exist
        if(os.path.isfile(default_paths.DEFAULT_RECO_TEMPLATE_FOLDER + template_filename)):
            # now open pkl file
            with open(default_paths.DEFAULT_RECO_TEMPLATE_FOLDER + template_filename, 'rb') as f:
                [reco_pipe_template] = pickle.load(f)

            # copy attributes
            self.settings = reco_pipe_template.settings
            self.job = reco_pipe_template.job
            self.job_list = reco_pipe_template.job_list
            self.analyze_job_list = reco_pipe_template.analyze_job_list
            self.analyze_enable = reco_pipe_template.analyze_enable
            self.analyze_display = reco_pipe_template.analyze_display

    def save_template(self, template_name):
        """
        Save current pipeline object as template.

        Parameters
        ----------
        template_name: string
            Template name
        """
        template_filename = template_name + ".pkl"
        log.info("saving current pipeline to template file [%s]..." % template_filename)
        # save to pkl file
        with open(default_paths.DEFAULT_RECO_TEMPLATE_FOLDER + template_filename, 'wb') as f:
            pickle.dump([self], f)


def remove_grids_from_all_figs():
    """Remove the grid in all axes for all open figures."""
    figs = [manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]

    for f in figs:
        for a in f.axes:
            a.grid(False)


def reco_spatial_select_profile(dcm_folders_list, legends_list, analyze_selectivity_range_list = [[800, 3550], [-10600, -7800], [-3650, 1850]]):
    """
    Reconstruct the VOI selection profile from data acquired using the 'VOI trace' mode (available in CMRR's MRS sequences)

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
    fig = plt.figure()
    fig.clf()
    axs = fig.subplots(2, 3)
    fig.canvas.set_window_title("mrs.reco.voi_pipeline.run")
    fig.suptitle("spatial selection profiles")
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
