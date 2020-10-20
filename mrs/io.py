#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stuff related to MRS data file reading.

@author: tangir
"""

import suspect
import suspect.io.twix as sit
import mapvbvd
import numpy as np
from mrs import sim
from mrs import log
import struct
from abc import ABCMeta
from datetime import datetime
import os
import pydicom
import hashlib

import pdb

class data_file_reader(metaclass=ABCMeta):
    """A virtual class used to read data and acquisition parameters from files. Subclasses depending on MRI constructors."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, data_fullfilepath):
        pass

    def is_rawdata(self):
        """Return true if data is read from a raw data file."""
        pass

    def get_number_rx_channels(self):
        """Return the number of channels."""
        pass

    def read_data(self):
        """Read MRS data from file and return a suspect's MRSData object."""
        pass

    def read_param_num(self, param_name, file_index=0):
        """Look for parameter in TWIX/DICOM data headers and return its float value."""
        pass

    def get_md5_hash(self):
        """Return a MD5 hash code of the whole binary data file content."""
        pass

    def get_nucleus(self):
        """Return the nucleus setting used to acquire the MRS data."""
        pass

    def get_patient_name(self):
        """Return the patient name field."""
        pass

    def get_patient_birthday(self):
        """Return the patient birthday field."""
        pass

    def get_patient_sex(self):
        """Return the patient sex field."""
        pass

    def get_patient_weight(self):
        """Return the patient weight field."""
        pass

    def get_patient_height(self):
        """Return the patient height field."""
        pass

    def get_sequence_name(self, file_index=0):
        """Return the acquisition sequence name."""
        pass

    def get_gating_mode(self):
        """Return the gating signal source."""
        pass

    def get_timestamp(self):
        """Return the acquisition start timestamp."""
        pass


class SIEMENS_data_file_reader(data_file_reader):
    """A class used to scrap parameters out of SIEMENS DICOM and TWIX files, sometimes in a very dirty way."""

    def __init__(self, data_fullfilepath):
        """
        Initialize by reading a DICOM or TWIX file in text mode and extracting the DICOM header with pydicom, if required.

        Parameters
        ----------
        data_fullfilepath : string
            Full path to the DICOM/TWIX file
        """
        super().__init__(data_fullfilepath)

        log.debug("checking data file path and extension...")
        data_filename, data_file_extension = os.path.splitext(data_fullfilepath)
        if(len(data_file_extension) == 0):
            # if empty extension, assuming the filename is not present in the path
            # lasy-mode where I copy-pasted only the folder paths
            # therefore, I will complete the path with the dicom name
            # which has always been "original-primary_e09_0001.dcm" for me up to now
            self.fullfilepath = data_fullfilepath + "/original-primary_e09_0001.dcm"
        else:
            self.fullfilepath = data_fullfilepath

        # detect DICOM or TWIX?
        data_filename, data_file_extension = os.path.splitext(self.fullfilepath)
        self.file_ext = data_file_extension.lower()

        # open file in text mode and save content
        log.debug("dumping file in ASCII mode...")
        f = open(self.fullfilepath, "r", encoding="ISO-8859-1")
        self.file_content_str = f.read()
        f.close()

        # if DICOM, use pydicom to extract basic header (SIEMENS hidden header not extracted)
        if(self.is_rawdata()):
            self.dcm_header = None
        else:
            log.debug("reading DICOM header...")
            self.dcm_header = pydicom.dcmread(self.fullfilepath)

        # freeze
        self.__isfrozen = True

    def is_rawdata(self):
        """
        Return true if data is read from a TWIX file.

        Returns
        -------
        self.file_ext == '.dat' : boolean
            True if TWIX file
        """
        return(self.file_ext == '.dat')

    def get_number_rx_channels(self):
        """
        Return the number of channels.

        Returns
        -------
        nchan : int
            Number of channels
        """
        nchan = self.read_param_num("lMaximumNofRxReceiverChannels")
        return(int(nchan))

    def read_data(self):
        """
        Read MRS data from file and return a suspect's MRSData object.

        Returns
        -------
        MRSData_obj : MRSData object
            MRS data read from DICOM/TWIX file
        """
        if(self.is_rawdata()):
            # try and read this TWIX file
            try:
                log.debug("reading TWIX file...")
                MRSData_obj = suspect.io.load_twix(self.fullfilepath)
            except:
                # well maybe it is broken, maybe the acquisition was interrupted
                # let's try to read it using this modified verion of suspect.io.load_twix
                log.debug("reading broken TWIX file...")
                MRSData_obj = self._load_broken_twix()
        elif(self.file_ext == '.dcm'):
            log.debug("reading DICOM file...")
            MRSData_obj = suspect.io.load_siemens_dicom(self.fullfilepath)
        else:
            log.error("unknown data file format!?")

        return(MRSData_obj)

    def read_param_num(self, param_name, file_index=0):
        """
        Look for parameter in TWIX/DICOM data headers and return its float value.

        Parameters
        ----------
        param_name : string
            Name of parameter to extract
        file_index : int
            Index in file from where we should start searching

        Returns
        -------
        param_val : float
            Value of parameter
        """
        # scrap out parameter value
        a = self.file_content_str.find(param_name, file_index)
        # could not find it?
        if(a < 0):
            log.warning("could not find parameter [%s]! :(" % param_name)
            return(np.nan)

        a = self.file_content_str.find("=", a + 1)
        b = self.file_content_str.find("\n", a + 1)
        param_val_str = self.file_content_str[(a + 1):b]
        param_val_str = param_val_str.strip()

        # try to convert to float
        try:
            param_val_float = float(param_val_str)
        except:
            # that did not work, look for next occurence
            param_val_float = self.read_param_num(param_name, b)

        # return it
        return(param_val_float)

    def get_md5_hash(self):
        """
        Return a MD5 hash code of the whole binary data file content.

        Returns
        -------
        h : string
            MD5 hexadecimal hash code
        """
        hasher = hashlib.md5()
        with open(self.fullfilepath, 'rb') as f:
            buf = f.read()
            hasher.update(buf)
        return(hasher.hexdigest())

    def get_nucleus(self):
        """
        Return the nucleus setting used to acquire the MRS data.

        Returns
        -------
        nucleus_str : string
            String representing nucleus. Example: "1H"
        """
        a = self.file_content_str.find("Info[0].tNucleus")
        a = self.file_content_str.find("=", a + 1)
        b = self.file_content_str.find("\n", a + 1)
        nucleus_str = self.file_content_str[(a + 1):b]
        nucleus_str = nucleus_str.replace('"', '')
        nucleus_str = nucleus_str.strip()

        return(nucleus_str)

    def get_patient_name(self):
        """
        Return the patient name field.

        Returns
        -------
        patient_name_str : string
            Patient name
        """
        if(self.is_rawdata()):
            # find patient name dirty way
            a = self.file_content_str.find("tPatientName")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find("\"", a)
            b = self.file_content_str.find("\"", a + 1)
            patient_name_str = self.file_content_str[(a + 1):b]
            patient_name_str = patient_name_str.strip()
        elif(self.file_ext == '.dcm'):
            patient_name_str = str(self.dcm_header.PatientName)
        else:
            log.error("unknown file extension!")

        return(patient_name_str)

    def get_patient_birthday(self):
        """
        Return the patient birthday field.

        Returns
        -------
        patient_birthday_datetime : datetime
            Patient birthday
        """
        if(self.is_rawdata()):
            # find patient birthday dirty way
            a = self.file_content_str.find("PatientBirthDay")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find("\"", a)
            birthday_str = self.file_content_str[(a + 1):(a + 9)]
            birthday_str = birthday_str.strip()
            if(birthday_str):
                patient_birthday_datetime = datetime.strptime(birthday_str, '%Y%m%d')
        elif(self.file_ext == '.dcm'):
            patient_birthday_datetime = datetime.strptime(str(self.dcm_header.PatientBirthDate), '%Y%m%d')
        else:
            log.error("unknown file extension!")

        return(patient_birthday_datetime)

    def get_patient_sex(self):
        """
        Return the patient sex field.

        Returns
        -------
        patient_sex_str : string
            Patient sex ('M', 'F' or 'O')
        """
        if(self.is_rawdata()):
            # find patient sex dirty way
            a = self.file_content_str.find("PatientSex")
            a = self.file_content_str.find("{", a)
            patient_sex_str = self.file_content_str[(a + 2):(a + 3)]
            patient_sex_int = int(patient_sex_str.strip())
            if(patient_sex_int == 2):
                patient_sex_str = 'M'
            elif(patient_sex_int == 1):
                patient_sex_str = 'F'
            elif(patient_sex_int == 3):
                patient_sex_str = 'O'
            else:
                log.error("unknown patient sex!")
        elif(self.file_ext == '.dcm'):
            patient_sex_str = str(self.dcm_header.PatientSex)
        else:
            log.error("unknown file extension!")

        return(patient_sex_str)

    def get_patient_weight(self):
        """
        Return the patient weight field.

        Returns
        -------
        patient_weight_kgs : float
            Patient weight in kg
        """
        if(self.is_rawdata()):
            # find patient weight dirty way
            a = self.file_content_str.find("flUsedPatientWeight")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find(">", a)
            b = self.file_content_str.find("}", a + 1)
            patient_weight_str = self.file_content_str[(a + 5):(b - 2)]
            patient_weight_kgs = float(patient_weight_str.strip())
        elif(self.file_ext == '.dcm'):
            patient_weight_kgs = float(self.dcm_header.PatientWeight)
        else:
            log.error("unknown file extension!")

        return(patient_weight_kgs)

    def get_patient_height(self):
        """
        Return the patient height field.

        Returns
        -------
        patient_height_m : float
            Patient height in meters
        """
        if(self.is_rawdata()):
            # find patient height dirty way
            a = self.file_content_str.find("flPatientHeight")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find(">", a)
            a = self.file_content_str.find(">", a + 1)
            b = self.file_content_str.find("}", a + 1)
            patient_height_str = self.file_content_str[(a + 6):(b - 3)]
            try:
                patient_height_m = float(patient_height_str.strip()) / 1000.0
            except:
                patient_height_m = np.nan
        elif(self.file_ext == '.dcm'):
            try:
                patient_height_m = float(self.dcm_header.PatientSize)
            except:
                patient_height_m = np.nan
        else:
            log.error("unknown file extension!")

        return(patient_height_m)

    def get_sequence(self):
        """
        Extract all parameters related to acquisition and sequence and return it as a sequence object.

        Returns
        -------
        seq : mrs.sim.mrs_sequence
            Sequence object
        """
        # echo time (btw already extracted within suspect, this is a bit redandent)
        te = self.read_param_num("alTE") / 1000.0
        log.debug("extracted echo time (%.2f)" % te)

        # repetion time (btw already extracted within suspect, this is a bit redandent)
        tr = self.read_param_num("alTR") / 1000.0
        log.debug("extracted repetition time (%.2f)" % tr)

        # averages
        na = int(self.read_param_num("lAverages"))
        log.debug("extracted number of averages (%d)" % na)

        # dummy scans
        # for some reason, lPreparingScans is sometimes not set in the header...
        ds = self.read_param_num("lPreparingScans")
        if(np.isnan(ds)):
            ds = 0 # I know, that is bad
        else:
            ds = int(ds)
        log.debug("extracted number of dummy scans (%d)" % ds)

        # reference voltage
        vref = self.read_param_num("flReferenceAmplitude")
        log.debug("extracted reference voltage (%.2fV)" % vref)

        # 1st order shim X
        shim_1st_X = self.read_param_num("lOffsetX")
        log.debug("extracted 1st order shim value for X (%.2f)" % shim_1st_X)

        # 1st order shim Y
        shim_1st_Y = self.read_param_num("lOffsetY")
        log.debug("extracted 1st order shim value for Y (%.2f)" % shim_1st_Y)

        # 1st order shim Z
        shim_1st_Z = self.read_param_num("lOffsetZ")
        log.debug("extracted 1st order shim value for Z (%.2f)" % shim_1st_Z)

        # 2nd / 3rd shims
        shims_2nd_3rd_list = []
        for i_shim in range(5):
            this_shim_val = self.read_param_num("alShimCurrent[%d]" % i_shim)
            shims_2nd_3rd_list.append(this_shim_val)
        log.debug("extracted 2nd/3rd shims (" + str(shims_2nd_3rd_list) + ")")

        # merge shims values into one vector
        shims_values = [shim_1st_X, shim_1st_Y, shim_1st_Z] + shims_2nd_3rd_list

        # nucleus (used for sequence object)
        nucleus = self.get_nucleus()
        log.debug("extracted nuclei (%s)" % nucleus)

        # number of points
        npts = int(self.read_param_num("lVectorSize"))
        log.debug("extracted number of points (%d)" % npts)

        # dwell time (btw already extracted within suspect, this is a bit redandent)
        dt = self.read_param_num("DwellTime") * 1e-9
        log.debug("extracted dwell time (%.2f)" % dt)

        # f0 frequency (btw already extracted within suspect, this is a bit redandent)
        f0 = self.read_param_num("Frequency") * 1e-6
        log.debug("extracted dwell time (%.6f)" % f0)

        # special timestamp
        ulTimeStamp_ms = self._get_sequence_timestamp()
        log.debug("found a timestamp (%.0f)" % ulTimeStamp_ms)

        # gating mode
        gss = self._get_gating_mode()
        log.debug("extracted gating mode (%d)" % gss.value)

        # effective acquisition time
        eff_acq_time = self._get_sequence_acquisition_time_effective()
        log.debug("extracted effective acquisition time (%.0f)" % eff_acq_time)

        # --- sequence-specific parameters ---
        sequence_name = self._get_sequence_name()
        log.debug("extracted sequence name (%s)" % sequence_name)

        if(sequence_name == "eja_svs_slaser"):
            log.debug("this a semi-LASER acquisition, let's extract some specific parameters!")

            # afp pulse (fake) flip angle
            pulse_laser_rfc_fa = self.read_param_num("adFlipAngleDegree[1]")
            log.debug("extracted LASER AFP refocussing pulse flip angle (%.2f)" % pulse_laser_rfc_fa)

            # afp pulse length
            pulse_laser_rfc_length = self.read_param_num("sWiPMemBlock.alFree[1]")
            log.debug("extracted LASER AFP refocussing pulse duration (%.2f)" % pulse_laser_rfc_length)

            # afp pulse R
            pulse_laser_rfc_r = self.read_param_num("sWiPMemBlock.alFree[49]")
            log.debug("extracted LASER AFP refocussing pulse R (%.2f)" % pulse_laser_rfc_r)

            # afp pulse N
            pulse_laser_rfc_n = self.read_param_num("sWiPMemBlock.alFree[48]")
            log.debug("extracted LASER AFP refocussing pulse N (%.2f)" % pulse_laser_rfc_n)

            # afp pulse voltage
            pulse_laser_rfc_voltage = self.read_param_num("aRFPULSE[1].flAmplitude")
            log.debug("extracted LASER AFP refocussing pulse voltage (%.2f)" % pulse_laser_rfc_voltage)

            # exc pulse duration
            pulse_laser_exc_length = self.read_param_num("sWiPMemBlock.alFree[24]")
            log.debug("extracted LASER exc. pulse length (%.2f)" % pulse_laser_exc_length)

            # exc pulse voltage
            pulse_laser_exc_voltage = self.read_param_num("aRFPULSE[0].flAmplitude")
            log.debug("extracted LASER excitation pulse voltage (%.2f)" % pulse_laser_exc_voltage)

            # spoiler length
            spoiler_length = self.read_param_num("sWiPMemBlock.alFree[12]")
            log.debug("extracted LASER spoiler length (%.2f)" % spoiler_length)

            # build sequence object
            sequence_obj = sim.mrs_seq_eja_svs_slaser(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0, pulse_laser_exc_length / 1000.0, pulse_laser_exc_voltage, pulse_laser_rfc_length / 1000.0, pulse_laser_rfc_fa, pulse_laser_rfc_r, pulse_laser_rfc_n, pulse_laser_rfc_voltage, spoiler_length / 1000.0)

        elif(sequence_name == "eja_svs_press"):
            sequence_obj = sim.mrs_seq_eja_svs_press(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0)

        elif(sequence_name == "eja_svs_steam"):
            sequence_obj = sim.mrs_seq_eja_svs_steam(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0)

        elif(sequence_name == "fid"):
            sequence_obj = sim.mrs_seq_fid(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0)

        elif(sequence_name == "svs_se"):
            sequence_obj = sim.mrs_seq_svs_se(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0)

        elif(sequence_name == "svs_st"):
            sequence_obj = sim.mrs_seq_svs_st(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0)

        elif(sequence_name == "svs_st_vapor_643"):
            log.debug("this is CMRR's STEAM sequence, let's extract some specific parameters!")

            # TM value
            TM_ms = self.read_param_num("alTD[0]") / 1000.0
            log.debug("extracted TM value (%.0fms)" % TM_ms)

            # build sequence object
            sequence_obj = sim.mrs_seq_svs_st_vapor_643(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0, TM_ms)

        elif(sequence_name == "bow_isis_15"):
            # TODO : create a sequence implementation for ISIS?
            sequence_obj = sim.mrs_seq_fid(te, tr, na, ds, nucleus, npts, 1.0 / dt, f0, vref, shims_values, ulTimeStamp_ms, gss, eff_acq_time, 1.0)

        elif(sequence_name is None):
            sequence_obj = None

        else:
            # unknown!
            log.error("ouch unknown sequence!")

        return(sequence_obj)

    def _get_sequence_name(self, file_index=0):
        """
        Return the acquisition sequence name.

        Parameters
        ----------
        file_index : int
            Index in file from where we should start searching

        Returns
        -------
        sequence_name : string
            Sequence used for the acquisition
        """
        # now some sequence-specific parameters
        a = self.file_content_str.find('%' + "CustomerSeq" + '%', file_index)
        # could not find it?
        if(a < 0):
            # maybe a siemens seq then ?
            a = self.file_content_str.find('%' + "SiemensSeq" + '%', file_index)
            # could not find it?
            if(a < 0):
                log.warning("could not find sequence name! :(")
                return(None)

        a = self.file_content_str.find('\\', a)
        b = self.file_content_str.find('"', a + 1)
        sequence_name = self.file_content_str[(a + 1):b]
        sequence_name = sequence_name.strip()

        # check we did not find a long garbage string
        if(len(sequence_name) > 16):
            # if so, go on searching in file
            sequence_name = self.get_sequence_name(b)

        return(sequence_name)

    def _get_sequence_timestamp(self):
        """
        Return the acquisition start timestamp.

        Returns
        -------
        ulTimeStamp_ms : float
            Acquisition start timestamp
        """
        # open TWIX file in binary mode to get MDH header
        # we do it here and not in the __init__ because we only do that once
        f = open(self.fullfilepath, "rb")
        binaryDump = f.read()
        hdr_len = struct.unpack("i", binaryDump[:4])
        sMDH = struct.unpack("iiiii", binaryDump[hdr_len[0]:hdr_len[0] + 20])
        # and extract timestamp
        ulTimeStamp = sMDH[3]
        ulTimeStamp_ms = float(ulTimeStamp * 2.5)
        return(ulTimeStamp_ms)

    def _get_gating_mode(self):
        """
        Return the gating signal source.

        Returns
        -------
        gts : gating_signal_source
            Type of signal used during gated acquisition
        """
        lsig1 = self.read_param_num("lSignal1")
        if(lsig1 == 1):
            return(sim.gating_signal_source.NO_GATING)
        elif(lsig1 == 16):
            return(sim.gating_signal_source.RESP_GATING)
        elif(lsig1 == 4):
            return(sim.gating_signal_source.CARDIAC_GATING)
        elif(lsig1 == 2):
            return(sim.gating_signal_source.CARDIAC_ECG)
        else:
            log.error("Sorry, I don't recognize what type of gating you used! :(")

    def _get_sequence_acquisition_time_effective(self):
        """
        Return the real acquisition duration, usefull especially when using gating. Works only on dicom files.

        Returns
        -------
        acq_time : float
            Acquisition effective duration (s)
        """
        if(self.is_rawdata()):
            log.warning("cannot extract effective acquisition time from a TWIX file :(")
            acq_time = np.nan
        elif(self.file_ext == '.dcm'):
            # the dicom header includes a SeriesInstanceUID string that contains the date and time of the start of the acquisition
            SeriesInstanceUID_splitted = self.dcm_header.SeriesInstanceUID.split('.')
            SeriesInstanceUID_str = max(SeriesInstanceUID_splitted, key=len)
            try:
                # so for some reason, I get very very rarely malformed SeriesInstanceUID where the data/time is placed at a different index...
                SeriesInstanceUID_datetime = datetime.strptime(SeriesInstanceUID_str[0:13], '%Y%m%d%H%M%S')
            except:
                return(np.nan)
            # the dicom header includes a Acquisition Time string that contains the time of the end of the acquisition
            AcquisitionTime_datetime = datetime.strptime(self.dcm_header.AcquisitionTime[0:6], '%H%M%S')
            # this should work except if we were scanning for several days or over midnight, ahah
            AcquisitionTime_datetime = AcquisitionTime_datetime.replace(day=SeriesInstanceUID_datetime.day)
            AcquisitionTime_datetime = AcquisitionTime_datetime.replace(month=SeriesInstanceUID_datetime.month)
            AcquisitionTime_datetime = AcquisitionTime_datetime.replace(year=SeriesInstanceUID_datetime.year)

            # et voila, difference to get real acquisition time
            acq_time = float((AcquisitionTime_datetime - SeriesInstanceUID_datetime).seconds)
        else:
            log.error("unknown file extension!")

        return(acq_time)

    def _load_broken_twix(self):
        """
        Read broken TWIX files.

        Returns
        -------
        broken_data : MRSData oject
            Data read from file
        """
        # try suspect lib
        # probably reads the header ok but the data will be corrupted
        broken_raw_data = self._load_broken_twix_suspect()

        # get the actual data
        broken_raw_data_numpy = self._load_broken_twix_vbvd()
        broken_raw_data_numpy = np.transpose(broken_raw_data_numpy, (2, 1, 0))

        # doubling dt: somehow vbvd reads twice less points...
        broken_raw_data.inherit(broken_raw_data_numpy)
        broken_raw_data._dt = broken_raw_data._dt * 2

        return(broken_raw_data)

    def _load_broken_twix_vbvd(self):
        """
        Read broken TWIX files using mapVBVD lib (python >=3.7 only).

        Returns
        -------
        broken_data : MRSData oject
            Data read from file
        """
        vbvd_rawdata = mapvbvd.mapVBVD(self.fullfilepath)
        vbvd_rawdata_sq = np.squeeze(vbvd_rawdata.image[:, :, 0, 0, 0, 0, 0, 0, 0, :, 0, 0, 0, 0, 0, 0])

        return(vbvd_rawdata_sq)

    def _load_broken_twix_suspect(self):
        """
        Read broken TWIX files using mapVBVD lib (python >=3.7 only).

        Returns
        -------
        builder.build_mrsdata() : MRSData object
            Resulting constructed MRSData object
        """
        with open(self.fullfilepath, 'rb') as fin:

            # we can tell the type of file from the first two uints in the header
            first_uint, second_uint = struct.unpack("II", fin.read(8))

            # reset the file pointer before giving to specific function
            fin.seek(0)

            # create a TwixBuilder object for the actual loader function to use
            builder = sit.TwixBuilder()

            # first four bytes are the size of the header
            header_size = struct.unpack("I", fin.read(4))[0]

            # read the rest of the header minus the four bytes we already read
            header = fin.read(header_size - 4)
            # for some reason the last 24 bytes of the header contain some junk that is not a string
            header = header[:-24].decode('latin-1')
            builder.set_header_string(header)

            # the way that vb files are set up we just keep reading scans until the acq_end flag is set

            while True:
                # start by keeping track of where in the file this scan started
                # this will be used to jump to the start of the next scan
                start_position = fin.tell()
                acq_end = False

                try:
                    # the first four bytes contain composite information
                    temp = struct.unpack("I", fin.read(4))[0]
                except:
                    acq_end = True

                if acq_end:
                    break

                # 25 LSBs contain DMA length (size of this scan)
                DMA_length = temp & (2 ** 26 - 1)
                # next we have the "pack" flag bit and the rest is PCI_rx
                # not sure what either of these are for but break them out in case
                # pack_flag = (temp >> 25) & 1
                # PCI_rx = temp >> 26

                meas_uid, scan_counter, time_stamp, pmu_time_stamp = struct.unpack("IIII", fin.read(16))

                # next long int is actually a lot of bit flags
                # a lot of them don't seem to be relevant for spectroscopy
                eval_info_mask = struct.unpack("Q", fin.read(8))[0]
                acq_end = eval_info_mask & 1
                rt_feedback = eval_info_mask >> 1 & 1
                hp_feedback = eval_info_mask >> 2 & 1
                sync_data = eval_info_mask >> 5 & 1
                # raw_data_correction = eval_info_mask >> 10 & 1
                # ref_phase_stab_scan = eval_info_mask >> 14 & 1
                # phase_stab_scan = eval_info_mask >> 15 & 1
                # sign_rev = eval_info_mask >> 17 & 1
                phase_correction = eval_info_mask >> 21 & 1
                # pat_ref_scan = eval_info_mask >> 22 & 1
                # pat_ref_ima_scan = eval_info_mask >> 23 & 1
                # reflect = eval_info_mask >> 24 & 1
                noise_adj_scan = eval_info_mask >> 25 & 1

                if acq_end:
                    break

                # if any of these flags are set then we should ignore the scan data
                if rt_feedback or hp_feedback or phase_correction or noise_adj_scan or sync_data:
                    fin.seek(start_position + DMA_length)
                    continue

                # now come the actual parameters of the scan
                num_samples, num_channels = struct.unpack("HH", fin.read(4))
                builder.set_num_channels(num_channels)

                # the loop counters are a set of 14 shorts which are used as indices
                # for the parameters an acquisition might loop over, including
                # averaging repetitions, COSY echo time increments and CSI phase
                # encoding steps
                # we have no prior knowledge about which counters might loop in a given
                # scan so we have to read in all scans and then sort out the data shape
                loop_counters = struct.unpack("14H", fin.read(28))

                cut_off_data, kspace_centre_column, coil_select, readout_offcentre = struct.unpack("IHHI", fin.read(12))
                time_since_rf, kspace_centre_line_num, kspace_centre_partition_num = struct.unpack("IHH", fin.read(8))

                # ice_program_params = struct.unpack("4H", fin.read(8))
                free_params = struct.unpack("4H", fin.read(8))

                # there are some dummy points before the data starts
                num_dummy_points = free_params[0]

                # we want our np to be the largest power of two within the num_samples - num_dummy_points
                npp = int(2 ** np.floor(np.log2(num_samples - num_dummy_points)))

                # slice_position = struct.unpack("7f", fin.read(28))

                # construct a numpy ndarray to hold the data from all the channels in this scan
                scan_data = np.zeros((num_channels, npp), dtype='complex')

                # loop over all the channels and extract data
                for channel_index in range(num_channels):
                    channel_id, ptab_pos_neg = struct.unpack("Hh", fin.read(4))
                    raw_data = struct.unpack("<{}f".format(num_samples * 2), fin.read(num_samples * 4 * 2))
                    # turn the raw data into complex pairs
                    data_iter = iter(raw_data)
                    complex_iter = (complex(r, -i) for r, i in zip(data_iter, data_iter))
                    try:
                        scan_data[channel_index, :] = np.fromiter(complex_iter, "complex64", num_samples)[num_dummy_points:(num_dummy_points + npp)]
                    except:
                        acq_end = True
                        break

                    # the vb format repeats all the header data for each channel in
                    # turn, obviously this is redundant so we read all but the channel
                    # index from the next header here
                    fin.read(124)

                if acq_end:
                    break

                builder.set_np(npp)
                # pass the data from this scan to the builder
                builder.add_scan(loop_counters, scan_data)

                # go to the next scan and the top of the loop
                fin.seek(start_position + DMA_length)

        return(builder.build_mrsdata())


def read_physio_file(physio_filepath):
    """
    Read physiological data signal from SIEMENS IDEA VB17 .resp log files and stores the resp trace, the trigger flags in the object.

    Parameters
    ----------
    physio_filepath: string
        Full absolute file path pointing to physio log file

    Returns
    -------
    resp_t : numpy array [n]
        Timestamp vector in ms since midnight
    resp_trigUp : numpy array [n]
        Logical vector indicating trigger on rising edges
    resp_trigDown : numpy array [n]
        Logical vector indicating trigger on falling edges
    resp_data : numpy array [n]
        Respiratory trace
    """
    log.debug("reading physio data...")

    # check file extension
    resp_log_filename_short, resp_log_file_extension = os.path.splitext(physio_filepath)
    if(resp_log_file_extension == ".resp"):
        # resp trace are sampled with a 50Hz frequency
        fs = 50.0  # Hz
    else:
        log.error("no data scan files specified!")

    # open the log file and the read the trace
    resp_file = open(physio_filepath, 'r')
    one_big_line = resp_file.readline()
    resp_data = np.fromstring(one_big_line, dtype=np.float, sep=' ')

    # jump some lines
    for i in range(9):
        this_line = resp_file.readline()

    # and read interesting timestamp parameters
    this_line = resp_file.readline()
    a = this_line.find("LogStartMDHTime:")
    b = this_line.find("\n", a)
    LogStartMDHTime_str = this_line[(a + len("LogStartMDHTime:")):b]
    LogStartMDHTime = float(LogStartMDHTime_str)

    '''
    this_line=resp_file.readline()
    a=this_line.find("LogStopMDHTime:")
    b=this_line.find("\n",a)
    LogStopMDHTime_str=this_line[(a + len("LogStopMDHTime:")):b]
    LogStopMDHTime=float(LogStopMDHTime_str)

    this_line=resp_file.readline()
    a=this_line.find("LogStartMPCUTime:")
    b=this_line.find("\n",a)
    LogStartMPCUTime_str=this_line[(a + len("LogStartMPCUTime:")):b]
    LogStartMPCUTime=float(LogStartMPCUTime_str)

    this_line=resp_file.readline()
    a=this_line.find("LogStopMPCUTime:")
    b=this_line.find("\n",a)
    LogStopMPCUTime_str=this_line[(a + len("LogStopMPCUTime:")):b]
    LogStopMPCUTime=float(LogStopMPCUTime_str)
    '''

    # let's clean the data

    # remove some useless points at the beginning
    resp_data = resp_data[4:]
    # remove the trigger up and down flags
    # keeping the trigger timestamps in separate logical vectors
    ts = 1 / fs
    resp_t = np.zeros(resp_data.shape)
    resp_trig5000 = np.full(resp_data.shape, False, dtype=bool)
    resp_trig6000 = np.full(resp_data.shape, False, dtype=bool)
    for it, t in enumerate(resp_data):
        if(t == 5000):
            resp_t[it] = resp_t[it - 1]
            resp_trig5000[it - 1] = True
        elif(t == 6000):
            resp_t[it] = resp_t[it - 1]
            resp_trig6000[it - 1] = True
        else:
            resp_t[it] = it * ts * 1000.0 + LogStartMDHTime  # ms

    # filter out all the trigger flags
    resp_t = resp_t[np.logical_and(resp_data != 5000, resp_data != 6000)]
    resp_trigUp = resp_trig5000[np.logical_and(resp_data != 5000, resp_data != 6000)]
    resp_trigDown = resp_trig6000[np.logical_and(resp_data != 5000, resp_data != 6000)]
    resp_data = resp_data[np.logical_and(resp_data != 5000, resp_data != 6000)]

    # filter out some weird flags
    resp_t = resp_t[resp_data != 5003]
    resp_trigUp = resp_trigUp[resp_data != 5003]
    resp_trigDown = resp_trigDown[resp_data != 5003]
    resp_data = resp_data[resp_data != 5003]

    # normalize data 0-1
    resp_data = resp_data - resp_data.min()
    resp_data = resp_data / resp_data.max()

    '''
    #remove extra first samples
    n_first_useless_samples=int(np.floor(((LogStartMPCUTime-LogStartMDHTime) / 1000.0)*fs));
    resp_t=resp_t[n_first_useless_samples:]
    resp_trigUp=resp_trigUp[n_first_useless_samples:]
    resp_trigDown=resp_trigDown[n_first_useless_samples:]
    resp_data=resp_data[n_first_useless_samples:]
    '''

    # return all this stuff
    return(resp_t, resp_trigUp, resp_trigDown, resp_data)
