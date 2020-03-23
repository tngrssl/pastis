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
import os
import datetime
import struct
import warnings
from termcolor import cprint
import suspect.io.twix as sit
import mrs.sim as sim

import pdb


class MRSData2(suspect.mrsobjects.MRSData):
    """A class based on suspect's MRSData to store MRS data."""

    def __new__(cls, data_filepath, coil_nChannels, physio_log_file, obj=[], dt=[], f0=[], te=[], ppm0=[], voxel_dimensions=[], transform=[], metadata=[], tr=[], timestamp=[], patient_birthyear=[], patient_sex=[], patient_name=[], patient_weight=[], patient_height=[], vref=[], shims=[], sequence_obj=[], display_apo_hz=[]):
        """
        Construct a MRSData2 object that inherits of Suspect's MRSData class. In short, the MRSData2 class is a copy of MRSData + my custom methods for post-processing. To create a MRSData2 object, you need give a path that points to a SIEMENS DICOM or a SIEMENS TWIX file.

        Parameters
        ----------
        data_filepath: string
            Full absolute file path pointing to the stored signal (DCM or TWIX file) or the folder assuming that a dcm file named "original-primary_e09_0001.dcm" is stored inside.
        coil_nChannels : int
            Number of RX channels
        physio_log_file : string
            Full absolute file path pointing to a IDEA VB17 respiratory log file
        obj,dt,f0,te,ppm0,voxel_dimensions,transform,metadata
            Please check suspect's MRSData class for those arguments
        tr : float
            TR in ms
        timestamp : float
            timestamp in ms
        patient_birthyear : int
            birthyear of patient
        patient_sex : int
            sex of patient (0:M 1:F)
        patient_name : string
            patient name
        patient_weight : float
            patient weight in kgs
        patient_height : float
            patient high in meters
        vref : float
            reference voltage (V)
        shims : list of floats
            list of shim voltages in volts
        sequence_obj : sim.mrs_sequence object
            sequence object
        display_apo_hz : float
            apodization factor applied when displaying the spectrum (Hz)

        Returns
        -------
        obj : MRSData2 numpy array [averages,channels,timepoints]
            Resulting constructed MRSData2 object
        """
        if(data_filepath == [] and coil_nChannels == []):
            # calling the parent class' constructor
            obj = super(suspect.mrsobjects.MRSData, cls).__new__(cls, obj, dt, f0, te, ppm0, voxel_dimensions, transform, metadata)
            # adding attributes
            obj._tr = tr
            obj._timestamp = timestamp
            obj._patient_birthyear = patient_birthyear
            obj._patient_sex = patient_sex
            obj._patient_name = patient_name
            obj._patient_weight = patient_weight
            obj._patient_height = patient_height
            obj._vref = vref
            obj._shims = shims
            obj._sequence = sequence_obj
            obj.display_apo = display_apo_hz
            # bye
            return(obj)

        # hello
        cprint(">> mrs.reco.MRSData2.__new__:", 'green')

        # init
        TR_ms = None
        year_int = None
        sex_int = None
        patient_name_str = None
        patient_weight_kgs = None
        patient_height_m = None
        vref_v = None
        shims_values = None
        sequence_obj = None
        display_apo = display_apo_hz
        ngrklmgnfkdlgnfdklmgnfdkmlgnfkdlmgnkfdlmgnkfdlmgnkl

        ulTimeStamp_ms = None

        print(" > checking data file path...")
        data_filename, data_file_extension = os.path.splitext(data_filepath)
        if(len(data_file_extension) == 0):
            # if empty extension, assuming the filename is not present in the path
            # lasy-mode where I copy-pasted only the folder paths
            # therefore, I will complete the path with the dicom name
            # which has always been "original-primary_e09_0001.dcm" for me up to now
            data_fullfilepath = data_filepath + "/original-primary_e09_0001.dcm"
        else:
            data_fullfilepath = data_filepath

        # find out if DICOM or TWIX file
        data_filename, data_file_extension = os.path.splitext(data_fullfilepath)
        if(data_file_extension.lower() == '.dat'):
            # twix!
            print(" > reading DAT / TWIX file...")
            print(data_fullfilepath)

            # try and read this TWIX file
            try:
                MRSData_obj = suspect.io.load_twix(data_fullfilepath)
            except:
                # well maybe it is broken, maybe the acquisition was interrupted
                MRSData_obj = load_broken_twix_vb(data_fullfilepath)

            # open TWIX file in text mode to scrap data out (dirty!!!)
            f = open(data_fullfilepath, "r", encoding="ISO-8859-1")
            twix_txtdata = f.read()
            f.close()

            # find TR value
            a = twix_txtdata.find("alTR[0]")
            a = twix_txtdata.find("=", a)
            b = twix_txtdata.find("\n", a)
            TR_str_us = twix_txtdata[(a + 2):b]
            if(TR_str_us.strip()):
                TR_ms = float(TR_str_us) / 1000.0
                print(" > extracted TR value (%.0fms)" % TR_ms)

            # find patient birthyear
            a = twix_txtdata.find("PatientBirthDay")
            a = twix_txtdata.find("{", a)
            a = twix_txtdata.find("\"", a)
            year_str = twix_txtdata[(a + 1):(a + 5)]
            if(year_str.strip()):
                year_int = int(year_str)
                print(" > extracted patient birthyear (%d)" % year_int)

            # find patient sex
            a = twix_txtdata.find("PatientSex")
            a = twix_txtdata.find("{", a)
            sex_str = twix_txtdata[(a + 2):(a + 3)]
            if(sex_str.strip()):
                sex_int = int(sex_str)
                print(" > extracted patient sex (%d)" % sex_int)

            # find patient name
            a = twix_txtdata.find("tPatientName")
            a = twix_txtdata.find("{", a)
            a = twix_txtdata.find("\"", a)
            b = twix_txtdata.find("\"", a + 1)
            patient_name_str = twix_txtdata[(a + 1):b]
            if(patient_name_str.strip()):
                print(" > extracted patient name (%s)" % patient_name_str)

            # find patient weight
            a = twix_txtdata.find("flUsedPatientWeight")
            a = twix_txtdata.find("{", a)
            a = twix_txtdata.find(">", a)
            b = twix_txtdata.find("}", a + 1)
            patient_weight_str = twix_txtdata[(a + 5):(b - 2)]
            if(patient_weight_str.strip()):
                patient_weight_kgs = float(patient_weight_str)
                print(" > extracted patient weight (%.2fkg)" % patient_weight_kgs)

            # find patient height
            a = twix_txtdata.find("flPatientHeight")
            a = twix_txtdata.find("{", a)
            a = twix_txtdata.find(">", a)
            a = twix_txtdata.find(">", a + 1)
            b = twix_txtdata.find("}", a + 1)
            patient_height_str = twix_txtdata[(a + 6):(b - 3)]
            if(patient_height_str.strip()):
                patient_height_m = float(patient_height_str) / 1000.0
                print(" > extracted patient height (%.2fm)" % patient_height_m)

            # find reference voltage
            a = twix_txtdata.find("TransmitterReferenceAmplitude")
            a = twix_txtdata.find("{", a)
            a = twix_txtdata.find(">", a)
            b = twix_txtdata.find("}", a + 1)
            vref_str = twix_txtdata[(a + 6):(b - 14)]
            if(vref_str.strip()):
                vref_v = float(vref_str)
                print(" > extracted reference voltage (%.2fV)" % vref_v)

            # find 1st order shim X
            a = twix_txtdata.find("lOffsetX")
            a = twix_txtdata.find("{", a)
            b = twix_txtdata.find("}", a + 1)
            shim_1st_X_str = twix_txtdata[(a + 4):b]
            if(shim_1st_X_str.strip()):
                shim_1st_X = float(shim_1st_X_str)
                print(" > extracted 1st order shim value for X (%.2f)" % shim_1st_X)
            else:
                shim_1st_X = np.nan

            # find 1st order shim Y
            a = twix_txtdata.find("lOffsetY")
            a = twix_txtdata.find("{", a)
            b = twix_txtdata.find("}", a + 1)
            shim_1st_Y_str = twix_txtdata[(a + 4):b]
            if(shim_1st_Y_str.strip()):
                shim_1st_Y = float(shim_1st_Y_str)
                print(" > extracted 1st order shim value for Y (%.2f)" % shim_1st_Y)
            else:
                shim_1st_Y = np.nan

            # find 1st order shim Z
            a = twix_txtdata.find("lOffsetZ")
            a = twix_txtdata.find("{", a)
            b = twix_txtdata.find("}", a + 1)
            shim_1st_Z_str = twix_txtdata[(a + 4):b]
            if(shim_1st_Z_str.strip()):
                shim_1st_Z = float(shim_1st_Z_str)
                print(" > extracted 1st order shim value for Z (%.2f)" % shim_1st_Z)
            else:
                shim_1st_Z = np.nan

            # find 2nd / 3rd shims
            a = twix_txtdata.find("alShimCurrent")
            a = twix_txtdata.find("{", a)
            b = twix_txtdata.find("}", a + 1)
            shims_2nd_3rd_str = twix_txtdata[(a + 12):(b - 31)]
            if(shims_2nd_3rd_str.strip()):
                shims_2nd_3rd_str = shims_2nd_3rd_str.split(" ")
                shims_2nd_3rd_list = [float(s) for s in shims_2nd_3rd_str]
                print(" > extracted 2nd / 3rd shims (" + str(shims_2nd_3rd_list) + ")")
            else:
                shims_2nd_3rd_list = np.nan

            # merge shims values into one vector
            shims_values = [shim_1st_X, shim_1st_Y, shim_1st_Z] + shims_2nd_3rd_list

            # now some sequence-specific parameters
            a = twix_txtdata.find("tSequenceFileName")
            a = twix_txtdata.find('\\', a)
            b = twix_txtdata.find('\n', a + 1)
            sequence_name = twix_txtdata[(a + 1):(b - 2)]
            print(" > extracted sequence name (%s)" % sequence_name)

            # nucleus (used for sequence object)
            a = twix_txtdata.find("Nucleus")
            a = twix_txtdata.find("{", a)
            a = twix_txtdata.find('"', a)
            b = twix_txtdata.find('"', a + 1)
            nucleus_str = twix_txtdata[(a + 1):b]
            nucleus = nucleus_str.strip()
            print(" > extracted nuclei (%s)" % nucleus)

            # find numer of points
            a = twix_txtdata.find("lVectorSize")
            a = twix_txtdata.find("{", a)
            b = twix_txtdata.find("}", a + 1)
            npts_str = twix_txtdata[(a + 4):b]
            if(npts_str.strip()):
                npts = int(npts_str)
                print(" > extracted number of points (%d)" % npts)
            else:
                npts = np.nan

            # open TWIX file in binary mode to get MDH header
            f = open(data_fullfilepath, "rb")
            binaryDump = f.read()
            hdr_len = struct.unpack("i", binaryDump[:4])
            sMDH = struct.unpack("iiiii", binaryDump[hdr_len[0]:hdr_len[0] + 20])
            # and extract timestamp
            ulTimeStamp = sMDH[3]
            ulTimeStamp_ms = float(ulTimeStamp * 2.5)
            print(" > extracted timestamp (%.0fms)" % ulTimeStamp_ms)

            if(sequence_name == "eja_svs_slaser"):
                print(" > This a semi-LASER acquisition, let's extract some specific parameters!")

                # find afp pulse (fake) flip angle
                a = twix_txtdata.find("FlipAngleDeg2")
                a = twix_txtdata.find('>', a)
                a = twix_txtdata.find('>', a + 1)
                b = twix_txtdata.find('}', a + 1)
                pulse_laser_rfc_fa_str = twix_txtdata[(a + 4):b]
                if(pulse_laser_rfc_fa_str.strip()):
                    pulse_laser_rfc_fa = float(pulse_laser_rfc_fa_str)
                    print(" > extracted LASER AFP refocussing pulse flip angle (%.2f)" % pulse_laser_rfc_fa)

                # find afp pulse length
                a = twix_txtdata.find("sWiPMemBlock.alFree[1]")
                a = twix_txtdata.find('=', a)
                b = twix_txtdata.find('\n', a + 1)
                pulse_laser_rfc_length_str = twix_txtdata[(a + 2):b]
                if(pulse_laser_rfc_length_str.strip()):
                    pulse_laser_rfc_length = float(pulse_laser_rfc_length_str)
                    print(" > extracted LASER AFP refocussing pulse duration (%.2f)" % pulse_laser_rfc_length)

                # find afp pulse R
                a = twix_txtdata.find("sWiPMemBlock.alFree[49]")
                a = twix_txtdata.find('=', a)
                b = twix_txtdata.find('\n', a + 1)
                pulse_laser_rfc_r_str = twix_txtdata[(a + 2):b]
                if(pulse_laser_rfc_r_str.strip()):
                    pulse_laser_rfc_r = float(pulse_laser_rfc_r_str)
                    print(" > extracted LASER AFP refocussing pulse R (%.2f)" % pulse_laser_rfc_r)

                # find afp pulse N
                a = twix_txtdata.find("sWiPMemBlock.alFree[48]")
                a = twix_txtdata.find('=', a)
                b = twix_txtdata.find('\n', a + 1)
                pulse_laser_rfc_n_str = twix_txtdata[(a + 2):b]
                if(pulse_laser_rfc_n_str.strip()):
                    pulse_laser_rfc_n = float(pulse_laser_rfc_n_str)
                    print(" > extracted LASER AFP refocussing pulse N (%.2f)" % pulse_laser_rfc_n)

                # find afp pulse voltage
                a = twix_txtdata.find("aRFPULSE[1].flAmplitude")
                a = twix_txtdata.find('=', a)
                b = twix_txtdata.find('\n', a + 1)
                pulse_laser_rfc_voltage_str = twix_txtdata[(a + 2):b]
                if(pulse_laser_rfc_voltage_str.strip()):
                    pulse_laser_rfc_voltage = float(pulse_laser_rfc_voltage_str)
                    print(" > extracted LASER AFP refocussing pulse voltage (%.2f)" % pulse_laser_rfc_voltage)

                # find exc pulse duration
                a = twix_txtdata.find("sWiPMemBlock.alFree[24]")
                a = twix_txtdata.find('=', a)
                b = twix_txtdata.find('\n', a + 1)
                pulse_laser_exc_length_str = twix_txtdata[(a + 2):b]
                if(pulse_laser_exc_length_str.strip()):
                    pulse_laser_exc_length = float(pulse_laser_exc_length_str)
                    print(" > extracted LASER exc. pulse length (%.2f)" % pulse_laser_exc_length)

                # find exc pulse voltage
                a = twix_txtdata.find("aRFPULSE[0].flAmplitude")
                a = twix_txtdata.find('=', a)
                b = twix_txtdata.find('\n', a + 1)
                pulse_laser_exc_voltage_str = twix_txtdata[(a + 2):b]
                if(pulse_laser_exc_voltage_str.strip()):
                    pulse_laser_exc_voltage = float(pulse_laser_exc_voltage_str)
                    print(" > extracted LASER excitation pulse voltage (%.2f)" % pulse_laser_exc_voltage)

                # find spoiler length
                a = twix_txtdata.find("sWiPMemBlock.alFree[12]")
                a = twix_txtdata.find('=', a)
                b = twix_txtdata.find('\n', a + 1)
                spoiler_length_str = twix_txtdata[(a + 2):b]
                if(spoiler_length_str.strip()):
                    spoiler_length = float(spoiler_length_str)
                    print(" > extracted spoiler length (%.2f)" % spoiler_length)

            elif(sequence_name == "svs_st_vapor_643"):
                print(" > This CMRR's STEAM sequence, let's extract some specific parameters!")

                # find TR value
                a = twix_txtdata.find("alTD[0]")
                a = twix_txtdata.find("=", a)
                b = twix_txtdata.find("\n", a)
                TM_str_us = twix_txtdata[(a + 2):b]
                if(TM_str_us.strip()):
                    TM_ms = float(TM_str_us) / 1000.0
                    print(" > extracted TM value (%.0fms)" % TM_ms)

        elif(data_file_extension.lower() == '.dcm' or data_file_extension.lower() == '.ima'):
            # dicom!
            print(" > reading DICOM file...")
            print(data_fullfilepath)
            MRSData_obj = suspect.io.load_siemens_dicom(data_fullfilepath)
        else:
            # unknown!
            raise Exception(' > ouch unknown file extension!')

        # add dimensions if needed
        if(MRSData_obj.ndim == 1):
            # this could be a single-shot single-rx signal in twix
            # or an averaged single-rx signal in dicom...

            # add 2 dimensions
            MRSData_obj = MRSData_obj.reshape((1, 1,) + MRSData_obj.shape)

        if(MRSData_obj.ndim == 2):
            # this could be a single-shot multi-rx signal
            # or a averaged single-rx signal...

            if(coil_nChannels == MRSData_obj.shape[0]):
                # beware, same number of averages / coil elements, which is which?
                print(" > WARNING: ambiguous data dimensions: " + str(MRSData_obj.shape) + "-> Assuming (" + str(MRSData_obj.shape[0]) + ") to be the coil channels!")
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

        print(" > read a " + str(MRSData_obj.shape) + " vector")

        # calling the parent class' constructor
        obj = super(suspect.mrsobjects.MRSData, cls).__new__(cls, MRSData_obj, MRSData_obj.dt, MRSData_obj.f0, MRSData_obj.te, MRSData_obj.ppm0, MRSData_obj.voxel_dimensions, MRSData_obj.transform, MRSData_obj.metadata)

        # create sequence object now we have all the parameters
        if(sequence_name == "eja_svs_slaser"):
            sequence_obj = sim.mrs_seq_eja_svs_slaser(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0, pulse_laser_exc_length / 1000.0, pulse_laser_exc_voltage, pulse_laser_rfc_length / 1000.0, pulse_laser_rfc_fa, pulse_laser_rfc_r, pulse_laser_rfc_n, pulse_laser_rfc_voltage, vref_v, spoiler_length / 1000.0)
        elif(sequence_name == "eja_svs_press"):
            sequence_obj = sim.mrs_sequence.mrs_seq_eja_svs_press(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)
        elif(sequence_name == "eja_svs_steam"):
            sequence_obj = sim.mrs_sequence.mrs_seq_eja_svs_steam(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)
        elif(sequence_name == "fid"):
            sequence_obj = sim.mrs_sequence.mrs_seq_fid(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)
        elif(sequence_name == "svs_st"):
            sequence_obj = sim.mrs_seq_svs_st(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)
        elif(sequence_name == "svs_st_vapor_643"):
            sequence_obj = sim.mrs_seq_svs_st_vapor_643(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0, TM_ms)
        else:
            # unknown!
            raise Exception(' > ouch unknown sequence!')

        # adding attributes
        obj._tr = TR_ms
        obj._patient_birthyear = year_int
        obj._patient_sex = sex_int
        obj._patient_name = patient_name_str
        obj._patient_weight = patient_weight_kgs
        obj._patient_height = patient_height_m
        obj._vref = vref_v
        obj._shims = shims_values
        obj._sequence = sequence_obj
        obj._timestamp = ulTimeStamp_ms

        # respiratory trace if any
        if(physio_log_file is None):
            obj._resp_trace = None
        else:
            # parse the physio file
            [rt, rup, rd, rr] = cls._read_physio_file(cls, physio_log_file)
            obj._resp_trace = [rt, rup, rd, rr]

        return(obj)

    def _read_physio_file(self, physio_filepath):
        """
        Read physiological data signal from SIEMENS IDEA VB17 .resp log files and stores the resp trace, the trigger flags in the object.

        Parameters
        ----------
        physio_filepath: string
            Full absolute file path pointing to phusio log file

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
        # hello
        cprint(">>> mrs.MRSData2._read_physio_file:", 'green')

        # check file extension
        resp_log_filename_short, resp_log_file_extension = os.path.splitext(
            physio_filepath)
        if(resp_log_file_extension == ".resp"):
            # resp trace are sampled with a 50Hz frequency
            fs = 50.0  # Hz
        else:
            raise Exception(" > no data scan files specified!")

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

    def __array_finalize__(self, obj):
        """
        Overload of special numpy array method called when playing around with stuff relative to object copy etc...

        Parameters
        ----------
        obj : MRSData2 numpy array [1,channels,timepoints]
        """
        super().__array_finalize__(obj)
        self._tr = getattr(obj, 'tr', None)
        self._patient_birthyear = getattr(obj, 'patient_birthyear', None)
        self._patient_sex = getattr(obj, 'patient_sex', None)
        self._resp_trace = getattr(obj, 'resp_trace', None)
        self._patient_name = getattr(obj, 'patient_name', None)
        self._patient_weight = getattr(obj, 'patient_weight', None)
        self._patient_height = getattr(obj, 'patient_height', None)
        self._vref = getattr(obj, 'vref', None)
        self._shims = getattr(obj, 'shims', None)
        self._timestamp = getattr(obj, 'timestamp', None)
        self._sequence = getattr(obj, 'sequence', None)
        self._pulse_laser_rfc_length = getattr(obj, 'pulse_laser_rfc_length', None)
        self._pulse_laser_rfc_r = getattr(obj, 'pulse_laser_rfc_r', None)
        self._pulse_laser_rfc_n = getattr(obj, 'pulse_laser_rfc_n', None)
        self._pulse_laser_rfc_voltage = getattr(obj, 'pulse_laser_rfc_voltage', None)
        self._pulse_laser_exc_length = getattr(obj, 'pulse_laser_exc_length', None)
        self._pulse_laser_exc_voltage = getattr(obj, 'pulse_laser_exc_voltage', None)

    def inherit(self, obj):
        """
        Overload of special suspect's MRSData method to import a numpy array in here.

        Parameters
        ----------
        obj : MRSData2 numpy array [1,channels,timepoints]
            Multi-channel reference signal
        """
        obj2 = super().inherit(obj)
        obj2._tr = getattr(self, 'tr', None)
        obj2._patient_birthyear = getattr(self, 'patient_birthyear', None)
        obj2._patient_sex = getattr(self, 'patient_sex', None)
        obj2._resp_trace = getattr(self, 'resp_trace', None)
        obj2._patient_name = getattr(self, 'patient_name', None)
        obj2._patient_weight = getattr(self, 'patient_weight', None)
        obj2._patient_height = getattr(self, 'patient_height', None)
        obj2._vref = getattr(self, 'vref', None)
        obj2._shims = getattr(self, 'shims', None)
        obj2._timestamp = getattr(self, 'timestamp', None)
        obj2._sequence = getattr(self, 'sequence', None)
        obj2._pulse_laser_rfc_length = getattr(self, 'pulse_laser_rfc_length', None)
        obj2._pulse_laser_rfc_r = getattr(self, 'pulse_laser_rfc_r', None)
        obj2._pulse_laser_rfc_n = getattr(self, 'pulse_laser_rfc_n', None)
        obj2._pulse_laser_rfc_voltage = getattr(self, 'pulse_laser_rfc_voltage', None)
        obj2._pulse_laser_exc_length = getattr(self, 'pulse_laser_exc_length', None)
        obj2._pulse_laser_exc_voltage = getattr(self, 'pulse_laser_exc_voltage', None)
        return(obj2)

    @property
    def tr(self):
        """
        Property get function for TR (ms).

        Returns
        -------
        self._tr : float
            TR in (ms)
        """
        return(self._tr)

    @property
    def patient_birthyear(self):
        """
        Property get function for patient_birthyear.

        Returns
        -------
        self._patient_birthyear : int
            birthyear of patient
        """
        return(self._patient_birthyear)

    @property
    def patient_sex(self):
        """
        Property get function for patient_sex.

        Returns
        -------
        self._patient_sex : int
            sex of patient (0:M 1:F)
        """
        return(self._patient_sex)

    @property
    def patient_name(self):
        """
        Property get function for patient_name.

        Returns
        -------
        self._patient_name : string
            name of patient
        """
        return(self._patient_name)

    @property
    def patient_weight(self):
        """
        Property get function for patient_weight.

        Returns
        -------
        self._patient_weight : float
            weight of patient in kgs
        """
        return(self._patient_weight)

    @property
    def patient_height(self):
        """
        Property get function for patient_height.

        Returns
        -------
        self._patient_height : float
            height of patient in meters
        """
        return(self._patient_height)

    @property
    def vref(self):
        """
        Property get function for vref.

        Returns
        -------
        self._vref : float
            reference voltage
        """
        return(self._vref)

    @property
    def shims(self):
        """
        Property get function for shims.

        Returns
        -------
        self._shims : list of floats
            shim voltages
        """
        return(self._shims)

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
    def timestamp(self):
        """
        Property get function for timestamp (ms).

        Returns
        -------
        self._timestamp : float
            timestamp in (ms)
        """
        return(self._timestamp)

    @property
    def resp_trace(self):
        """
        Property get function for resp_trace.

        Returns
        -------
        self._resp_trace : list
            Respiratory trace stuff
        """
        return(self._resp_trace)

    def _print_progress_bar(self, current_step, number_steps=None):
        """
        Print a small and pretty ASCII progress bar in the terminal. There is maybe some n / n + 1 bug but I don't care.

        Parameters
        ----------
        current_step : int
            Integer describing current status or number of achieved steps
        number_steps : int
            Integer describing maximum number of steps to achieve task
        """
        max_counter = 30
        if(current_step == 0 and number_steps is not None):
            # initialization
            self._progress_bar_counter = 0
            self._progress_bar_counter_max = number_steps
            bar_str = ""
            for i in range(max_counter):
                bar_str = bar_str + '░'
            print(bar_str, end="", flush=True)
        else:
            # growing bar?
            current_bar_counter = int(np.ceil(current_step / self._progress_bar_counter_max * max_counter))
            if(current_bar_counter > self._progress_bar_counter):
                bar_str = ""
                for i in range(max_counter):
                    bar_str = bar_str + '\b'
                for i in range(current_bar_counter):
                    bar_str = bar_str + '█'
                for i in range(max_counter - current_bar_counter):
                    bar_str = bar_str + '░'

                print(bar_str, end="", flush=True)
                self._progress_bar_counter = current_bar_counter

    def cs_displacement(self):
        """
        Calculate chemical shift displacement in three directions.

        Returns
        -------
        csd : numpy array of floats
            Estimated CS displacement in %/ppm
        """
        # hello
        cprint(">> mrs.reco.MRSData2.cs_displacement:", 'green')

        # init
        csd = [None, None, None]
        df_abs_Hz = self.f0  # yes, 1ppm==f0[MHz]/1e6=297Hz at 7T

        if(self.sequence_name == sim.mrs_sequence.SLASER):
            # assuming X-Y-Z is done with 90-180-180
            print(" > estimating CS displacement for semiLASER: assuming (90x)-(180y)-(180z)!...")

            # X selection done with asymmetric 90° pulse
            # we do not know much about this pulse. We can only say it is 3.4kHz large if the duration is 2ms
            print(" > estimating CS displacement for (90x): this pulse is the weird asymmetric one, we have no idea what it is exactly!")
            if(self.pulse_laser_exc_length == 2000.0):
                print(" > since its duration is 2ms here, we assume, according to Oz & Tkac, MRM 65:901-910 (2011), that its bandwidth is 3.4kHz.")
                bw_x_Hz = 3400.0
                grad_x_Hz_m = bw_x_Hz / (self.voxel_size[0] * 0.001)
                d_x_m = 1000.0 * df_abs_Hz / grad_x_Hz_m
                d_x_prct = d_x_m / self.voxel_size[0] * 100.0
            else:
                print(" > since its duration is not 2ms here, we do not know its bandwith. Therefore, no way to calculate the CS displacement for this axis, sorry ;)")
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
            print(" > no idea how to calculate CS displacement for this sequence...")

        return(csd)

    def correct_fidmodulus(self):
        """
        Calculate absolute mode of FID signals.

        Returns
        -------
        s : MRSData2 numpy array [averages,channels,timepoints]
            Resulting FID modulus data stored in a MRSData2 object
        """
        # hello
        cprint(">> mrs.reco.MRSData2.correct_fidmodulus:", 'green')

        # init
        s = self.copy()

        # return magnitude
        s = s.inherit(np.abs(s))
        return(s)

    def correct_phase(self, s_ref, peak_range, weak_ws, high_snr, phase_order, phase_offset=0.0, display=False, display_range=[1, 6]):
        """
        Well, that's a big one but basically it rephases the signal of interest.

        >In the case of multi-channel acquisition, the phase will be estimated for each channel, for each average using phase time evolution estimated on reference signal.
        >If the reference signal is not specified and weak water suppression was performed, 0th order phase correction will be done using the first point in the fid.
        >If strong water suppression was performed, 0th order phase correction will be done using the phase of a peak in the spectrum; the chemical shift range to find this peak is specified by the user.
        >>Note that those last two approaches can be performed for each average in the case of high SNR (rare) or by averaging all the scans for each channel in the case of lower SNR.

        Parameters
        ----------
        s_ref : MRSData2 numpy array [1,channels,timepoints]
            Multi-channel reference signal used to estimate the phase variation in time
        peak_range : list [2]
            Range in ppm used to peak-pick and estimate a phase
        weak_ws : boolean
            Weak water suppression, enough to have substantial water signal contributing to the 1st time point for phasing
        high_snr : boolean
            High SNR, enough to process each channel separately
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
        s : MRSData2 numpy array [averages,channels,timepoints]
            Resulting phased data stored in a MRSData2 object
        """
        # hello
        cprint(">> mrs.reco.MRSData2.correct_phase:", 'green')

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
        if(s_ref is None):
            # we do not have a ref scan, so method #2, #3, #4, #5
            # that depends on water intensity and water suppression
            if(weak_ws):
                phase_method = 3
            else:
                phase_method = 5

            # and on SNR
            if(high_snr):
                phase_method -= 1
        else:
            # we have a ref scan, so method #0 or #1
            phase_method = 0 + only_0th_order
            t_ref = s_ref.time_axis()

        # display chosen method
        print(" > " + list_phase_method[phase_method])
        print(" > phasing... ", end="", flush=True)

        if(display):
            # prepare subplots
            fig = plt.figure(100)
            fig.clf()
            axs = fig.subplots(2, 3, sharex='col')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_phase")

        # for each channel
        self._print_progress_bar(0, s.shape[1] * s.shape[0])
        for c in range(0, s.shape[1]):

            if(phase_method == 0 or phase_method == 1):
                # time-domain phase of reference signal for this channel
                sp_ref = np.angle(s_ref[0, c, :])

                if(display):
                    # display reference FID
                    axs[0, 0].cla()
                    axs[0, 0].plot(t_ref, np.real(s_ref[0, c, :]), linewidth=1, label='real part')
                    axs[0, 0].plot(t_ref, np.imag(s_ref[0, c, :]), linewidth=1, label='imag part')
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
                sf_avg = s[:, c, :].mean(axis=0).spectrum()

                # find maximum peak in range and its phase
                ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
                sf_avg_abs_masked = np.abs(sf_avg)
                sf_avg_abs_masked[ippm_peak_range] = 0
                ippm_peak = np.argmax(sf_avg_abs_masked)
                ppm_peak = ppm[ippm_peak]
                phase_peak_avg = np.angle(sf_avg[ippm_peak])
                if(ippm_peak == 0):
                    raise Exception(" > no peak found in specified ppm range or badly phased data!")
                if(c == 0):
                    print(" > measuring phase at %0.2fppm on 1st channel!" % ppm_peak)

            # now, for each average in meta signal acquired with this channel
            for a in range(0, s.shape[0]):

                # this spectrum
                sf = s[a, c, :].spectrum()

                if(phase_method == 0):
                    # correct phase using reference time-domain phase estimation
                    s_phased[a, c, :] = s[a, c, :] * np.exp(-1j * (sp_ref + phase_offset))
                elif(phase_method == 1):
                    # correct phase using first point of reference time-domain phase estimation
                    s_phased[a, c, :] = s[a, c, :] * np.exp(-1j * (sp_ref[0] + phase_offset))
                elif(phase_method == 2):
                    # phase of first point of this fid
                    phase_fid = np.angle(s[a, c, 0])
                    # and apply it to correct the spectrum
                    s_phased[a, c, :] = s[a, c, :] * np.exp(-1j * (phase_fid + phase_offset))
                elif(phase_method == 3):
                    # apply to correct the spectrum
                    s_phased[a, c, :] = s[a, c, :] * np.exp(-1j * (phase_fid_avg + phase_offset))
                elif(phase_method == 4):
                    # find maximum peak in range and its phase
                    ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
                    sf_abs_masked = np.abs(sf)
                    sf_abs_masked[ippm_peak_range] = 0
                    ippm_peak = np.argmax(sf_abs_masked)
                    phase_peak = np.angle(sf[ippm_peak])

                    # apply phase to spectrum
                    sf_phased = s[a, c, :].spectrum() * np.exp(-1j * (phase_peak + phase_offset))
                    # ifft back
                    s_phased[a, c, :] = np.fft.ifft(np.fft.ifftshift(sf_phased))
                elif(phase_method == 5):
                    # apply phase to spectrum
                    sf_phased = s[a, c, :].spectrum() * np.exp(-1j * (phase_peak_avg + phase_offset))
                    # ifft back
                    s_phased[a, c, :] = np.fft.ifft(np.fft.ifftshift(sf_phased))

                if(display):
                    # display original meta FID
                    axs[0, 1].cla()
                    axs[0, 1].plot(t, np.real(s[a, c, :]), linewidth=1, label='real part')
                    axs[0, 1].plot(t, np.imag(s[a, c, :]), linewidth=1, label='imag part')
                    axs[0, 1].set_xlabel('time (s)')
                    axs[0, 1].set_ylabel('original')
                    axs[0, 1].grid('on')
                    axs[0, 1].legend()
                    axs[0, 1].set_title("Meta channel #" + str(c + 1) + " average #" + str(a + 1))

                    # display original meta spectrum
                    axs[0, 2].cla()
                    axs[0, 2].plot(ppm, np.real(
                        s[a, c, :].spectrum()), linewidth=1, label='real part')
                    if(phase_method == 3):
                        axs[0, 2].plot(ppm[ippm_peak], np.real(sf[ippm_peak]), 'ro')
                    if(phase_method == 4):
                        axs[0, 2].plot(ppm[ippm_peak], np.real(sf[ippm_peak]), 'ro')
                    axs[0, 2].set_xlim(display_range[1], display_range[0])
                    axs[0, 2].set_xlabel('time (s)')
                    axs[0, 2].set_ylabel('original')
                    axs[0, 2].grid('on')
                    axs[0, 2].legend()
                    axs[0, 2].set_title("Meta channel #" + str(c + 1) + " average #" + str(a + 1))

                    # display corrected meta FID
                    axs[1, 1].cla()
                    axs[1, 1].plot(t, np.real(s_phased[a, c, :]), linewidth=1, label='real part')
                    axs[1, 1].set_xlabel('time (s)')
                    axs[1, 1].set_ylabel('corrected')
                    axs[1, 1].grid('on')
                    axs[1, 1].legend()

                    # display corrected meta spectrum
                    axs[1, 2].cla()
                    axs[1, 2].plot(ppm, s_phased[a, c, :].spectrum().real, linewidth=1, label='real part')
                    axs[1, 2].set_xlim(display_range[1], display_range[0])
                    axs[1, 2].set_xlabel('frequency (ppm)')
                    axs[1, 2].set_ylabel('corrected')
                    axs[1, 2].grid('on')
                    axs[1, 2].legend()

                    fig.tight_layout()
                    plt.pause(0.1)

                self._print_progress_bar(c * s.shape[0] + a + 1)

        print(" done.")

        return(self.inherit(s_phased))

    def correct_combine_channels(self, s_ref, phasing=True, channels_onoff=[True]):
        """
        Recombine Rx channels using SVD stuff. If no reference signal is specified, the recombination weights will be calculated from the signal of interest (not optimal).

        Parameters
        ----------
        s_ref : MRSData2 numpy array [averages,channels,timepoints]
            Multi-channel reference data (usually non water suppressed) stored in a MRSData2 object
        phasing : boolean
            Allow 0th order phasing during channel recombine or not
        channels_onoff : boolean list [nChannels]
            Binary weights to apply to each channel forexample to turn off some of them

        Returns
        -------
        s_combined : MRSData2 numpy array [averages,timepoints]
            Resulting channel combined data stored in a MRSData2 object
        """
        # hello
        cprint(">> mrs.reco.MRSData2.correct_combine_channels:", 'green')

        # init
        s = self.copy()
        if(s.shape[1] == 1):
            print(" > no need to recombine this!")
            s_combined = np.mean(s, axis=1)
            print(" > reshaped to " + str(s_combined.shape))
        else:
            if(phasing):
                if(s_ref is not None):
                    print(
                        " > WITH reference scan & phasing (original suspect code)...", end="", flush=True)
                    weights = suspect.processing.channel_combination.svd_weighting(s_ref.mean(axis=0))
                else:
                    print(
                        " > WITHOUT reference scan & phasing (original suspect code)...", end="", flush=True)
                    s_dirty_mean = np.mean(s, axis=0)
                    weights = suspect.processing.channel_combination.svd_weighting(s_dirty_mean)
            else:
                if(s_ref is not None):
                    print(" > WITH reference scan & NO phasing...",
                          end="", flush=True)
                    p, _, v = np.linalg.svd(s_ref.mean(axis=0), full_matrices=False)
                    channel_weights = p[:, 0].conjugate()
                    weights = -channel_weights / np.sum(np.abs(channel_weights))
                else:
                    print(" > WITHOUT reference scan & NO phasing...",
                          end="", flush=True)
                    s_dirty_mean = np.mean(s, axis=0)
                    p, _, v = np.linalg.svd(s_dirty_mean, full_matrices=False)
                    channel_weights = p[:, 0].conjugate()
                    weights = -channel_weights / np.sum(np.abs(channel_weights))

            # turn off some channels ?
            channels_onoff_np = np.array(channels_onoff)
            if((channels_onoff_np == False).any()):
                channels_onoff_float = channels_onoff_np.astype(float)
                print(" > playing with channel weights...")
                print(channels_onoff_float)
                weights = weights * channels_onoff_float

            s_combined = suspect.processing.channel_combination.combine_channels(
                s, weights)
            print(" done.")

        return(self.inherit(s_combined))

    def correct_zerofill(self, nPoints_final, display=True, display_range=[1, 6]):
        """
        Zero-fill MRS data signals.

        Parameters
        ----------
        nPoints_final : int
            Final number of points
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display

        Returns
        -------
        s_zf : MRSData2 numpy array [averages,timepoints]
            Resulting zero-filled data stored in a MRSData2 object
        """
        # hello
        cprint(">> mrs.reco.MRSData2.correct_zerofill:", 'green')

        # init
        s = self.copy()
        print(" > zero-filling data...", end="", flush=True)
        nZeros = nPoints_final - s.np
        s_new_shape = list(s.shape)
        s_new_shape[-1] = nZeros
        s_zf = s.inherit(np.concatenate((s, np.zeros(s_new_shape)), axis=s.ndim - 1))

        if(display):
            s_disp = s.copy()
            s_zf_disp = s_zf.copy()
            while(s_disp.ndim > 1):
                s_disp = np.mean(s_disp, axis=0)
                s_zf_disp = np.mean(s_zf_disp, axis=0)

            fig = plt.figure(110)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='row', sharey='row')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_zerofill")

            axs[0, 0].plot(s_disp.time_axis(), np.real(s_disp), 'k-', linewidth=1)
            axs[0, 0].set_xlabel('time (s)')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(s_zf_disp.time_axis(), np.real(s_zf_disp), 'k-', linewidth=1)
            axs[0, 1].set_xlabel('time (s)')
            axs[0, 1].set_ylabel('zero-filled')
            axs[0, 1].grid('on')

            axs[1, 0].plot(s_disp.frequency_axis_ppm(), s_disp.spectrum().real, 'k-', linewidth=1)
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('original')
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].grid('on')

            axs[1, 1].plot(s_zf_disp.frequency_axis_ppm(), s_zf_disp.spectrum().real, 'k-', linewidth=1)
            axs[1, 1].set_xlabel("chemical shift (ppm)")
            axs[1, 1].set_ylabel('zero-filled')
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].grid('on')

            fig.tight_layout()
            # plt.pause(0.1)

        print("done.")
        return(self.inherit(s_zf))

    def _build_moving_average_data(self, nAvgWindow, beSilent=False):
        """
        Build moving average data in the average dimension. Usefull for the analyse_peak and correct_realign functions.

        Parameters
        ----------
        nAvgWindow : int
            Size of the moving average window
        beSilent : boolean
            No output in console (True)

        Returns
        -------
        s_ma : MRSData2 numpy array [averages,timepoints]
            Resulting moving average data stored in a MRSData2 object. The number of averages is the same as the original data BUT each of those average is actually an average of nAvgWindow spectra.
        """
        # hello
        if(not beSilent):
            print(">>> mrs.reco.MRSData2._build_moving_average_data:")

        # init
        if(not beSilent):
            print("  > using a moving window of (" + str(nAvgWindow) + ") averages!")
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

    def _analyse_peak(self, peak_range, relativeToMean=False):
        """
        Analyse a peak in the spectrum by estimating its amplitude, linewidth, frequency shift and phase for each average.

        Parameters
        ----------
        peak_range : array [2]
            Range in ppm used to analyse peak phase when no reference signal is specified
        relativeToMean : boolean
            The variation calculation needs a reference. Use the average over the whole data (True) or the first shot (False)

        Returns
        -------
        peak_trace : numpy array [averages,4]
            Peak changes (amplitude, linewidth, frequency and phase) for each average in raw data
        peak_trace_rel : numpy array [averages,4]
            Peak relative changes (amplitude, linewidth, frequency and phase) for each average in raw data
        """
        # hello
        print(">>> mrs.reco.MRSData2._analyse_peak:")

        # first, find peak of interest in range
        ppm = self.frequency_axis_ppm()
        s_avg = np.mean(self, axis=0)
        sf_avg_abs = np.abs(s_avg.spectrum())
        ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
        sf_masked = sf_avg_abs.copy()
        sf_masked[ippm_peak_range] = 0
        ippm_peak_avg = np.argmax(sf_masked, axis=0)
        if(ippm_peak_avg == 0):
            raise Exception(" > no peak found in specified ppm range or badly phased data!")
        ppm_peak_avg = ppm[ippm_peak_avg]
        print("  > found peak of interest at %0.2fppm!" % ppm_peak_avg)

        # for each average in moving averaged data
        print("  > analyzing... ", end="", flush=True)
        peak_trace = np.zeros([self.shape[0], 4])
        self._print_progress_bar(0, self.shape[0])
        for a in range(0, self.shape[0]):

            # first, measure shift in ppm
            sf_masked = np.abs(self[a, :].spectrum())
            sf_masked[ippm_peak_range] = 0
            ippm_peak = np.argmax(sf_masked)
            if(ippm_peak == 0):
                raise Exception(" > no peak found in specified ppm range or badly phased data!")
            ppm_peak = ppm[ippm_peak]
            peak_trace[a, 2] = ppm_peak

            # estimate amplitude
            sf_masked = np.real(self[a, :].spectrum())
            sf_masked[ippm_peak_range] = 0
            ippm_peak = np.argmax(sf_masked)
            if(ippm_peak == 0):
                raise Exception(" > no peak found in specified ppm range or badly phased data!")
            ppm_peak = ppm[ippm_peak]
            amp_peak = sf_masked[ippm_peak]
            peak_trace[a, 0] = amp_peak

            # estimate linewidth in Hz
            ippm_half_peak = np.where(sf_masked > amp_peak / 2.0)
            ippm_min = np.min(ippm_half_peak)
            ippm_max = np.max(ippm_half_peak)
            dppm = np.abs(ppm[ippm_max] - ppm[ippm_min])
            peak_trace[a, 1] = dppm * self.f0

            # estimate phase in rad
            sf_masked = self[a, :].spectrum()
            peak_trace[a, 3] = np.angle(sf_masked[ippm_peak])

            self._print_progress_bar(a)

        # normalize stuff
        peak_trace_rel = np.zeros([self.shape[0], 4])
        if(relativeToMean):
            peak_trace_rel[:, 0] = peak_trace[:, 0] / peak_trace[:, 0].mean() * 100 - 100
            peak_trace_rel[:, 1] = peak_trace[:, 1] - peak_trace[:, 1].mean()
            peak_trace_rel[:, 2] = peak_trace[:, 2] - peak_trace[:, 2].mean()
            peak_trace_rel[:, 3] = peak_trace[:, 3] - peak_trace[:, 3].mean()
        else:
            peak_trace_rel[:, 0] = peak_trace[:, 0] / peak_trace[0, 0] * 100 - 100
            peak_trace_rel[:, 1] = peak_trace[:, 1] - peak_trace[0, 1]
            peak_trace_rel[:, 2] = peak_trace[:, 2] - peak_trace[0, 2]
            peak_trace_rel[:, 3] = peak_trace[:, 3] - peak_trace[0, 3]

        print(" done.")
        return(peak_trace, peak_trace_rel)

    def analyse_physio(self, peak_range, delta_time_range, display=True):
        """
        Analyse the physiological signal and try to correlate it to a peak amplitude, linewidth, frequency shift and phase variations.

        Parameters
        ----------
        peak_range : array [2]
            Range in ppm used to analyse peak phase when no reference signal is specified
        delta_time_range : float
            Range in ms used to correlate / match the NMR and the physiological signal. Yes, since we are not really sure of the start timestamp we found in the TWIX header, we try to match perfectly the two signals.
        display : boolean
            Display correction process (True) or not (False)
        """
        # hello
        cprint(">> mrs.reco.MRSData2.analyse_physio:", 'green')

        if(self.resp_trace is None):
            # no physio signal here, exiting
            return()

        # perform peak analysis
        params_trace, params_trace_rel = self._analyse_peak(peak_range, False)

        # physio signal
        resp_t = self.resp_trace[0]
        resp_s = self.resp_trace[3]

        # init
        mri_t = np.linspace(self.timestamp, self.timestamp + self.tr * params_trace.shape[0], params_trace.shape[0])
        dt_array = np.arange(-delta_time_range / 2.0, delta_time_range / 2.0, 1.0)
        cc_2d = np.zeros([dt_array.shape[0], 4])

        print("  > correlating signals... ", end="", flush=True)
        self._print_progress_bar(0, dt_array.shape[0])

        # shift signal and calculate corr coeff
        for idt, dt in enumerate(dt_array):
            # build time scale
            this_resp_t_interp = mri_t.copy() + dt
            this_resp_s_interp = np.interp(this_resp_t_interp, resp_t, resp_s)

            # now crop the signals to have the same length
            final_length = min(
                this_resp_s_interp.shape[0], params_trace.shape[0])
            this_mri_t = mri_t[0:final_length]
            this_params_trace = params_trace[0:final_length, :]
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

            self._print_progress_bar(idt)

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
        print(" done.")

        # some info in the term
        st_ms = self.timestamp
        st_str = datetime.datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        print(">> Data timestamp=\t" + str(st_ms) + "ms\t" + st_str)
        print(">> Best start time for...")

        st_ms = self.timestamp + best_dt_per_par[0]
        st_str = datetime.datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        print(" > amplitude=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_per_par[1]
        st_str = datetime.datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        print(" > linewidth=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_per_par[2]
        st_str = datetime.datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        print(" > frequency=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_per_par[3]
        st_str = datetime.datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        print(" > phase=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_all
        st_str = datetime.datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        print(" > total=\t\t" + str(st_ms) + "ms\t" + st_str)

        imaxR = np.argmax(best_dt_per_par)
        best_dt = best_dt_per_par[imaxR]
        st_ms = self.timestamp + best_dt
        st_str = datetime.datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        print(" > max R for=\t\t" + str(st_ms) + "ms\t" + st_str)

        # time shift the signals with optimal shift
        # build time scale
        this_resp_t_interp = mri_t.copy() + best_dt
        this_resp_s_interp = np.interp(this_resp_t_interp, resp_t, resp_s)

        # now crop the signals to have the same length
        final_length = min(this_resp_s_interp.shape[0], params_trace.shape[0])
        this_mri_t = mri_t[0:final_length]
        this_params_trace = params_trace[0:final_length, :]
        this_resp_s_interp = this_resp_s_interp[0:final_length]

        # evaluate correlation coeff.
        # for each parameter
        cc_final = np.zeros(4)
        for p in range(4):
            # estimate some R coeff
            cc = np.corrcoef(this_resp_s_interp, this_params_trace[:, p])
            cc_final[p] = cc[0, 1]

        # now let's talk about FFT
        print(">> FFT analysis...", end="", flush=True)
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

        print(" done.")

        if(display):
            # display the cross-correlation plots
            fig = plt.figure(120)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyse_physio_1")

            p = 0
            for ix in range(2):
                for iy in range(2):
                    axs[ix, iy].plot(dt_array, cc_2d[:, p], '-', linewidth=1)
                    axs[ix, iy].plot([best_dt_per_par[p], best_dt_per_par[p]], [cc_2d[:, p].min(), cc_2d[:, p].max()], 'r-', linewidth=1)
                    axs[ix, iy].plot([best_dt, best_dt], [cc_2d[:, p].min(), cc_2d[:, p].max()], 'r--', linewidth=1)
                    axs[ix, iy].set_xlabel('time shift (ms)')
                    axs[ix, iy].grid('on')
                    p = p + 1

            axs[0, 0].set_ylabel('R amplitude vs. resp.')
            axs[0, 1].set_ylabel('R linewidth. vs. resp.')
            axs[1, 0].set_ylabel('R frequency vs. resp.')
            axs[1, 1].set_ylabel('R phase vs. resp.')
            fig.tight_layout()

            # display time signals
            fig = plt.figure(121)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyse_physio_2")

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
            fig.tight_layout()

            # display correlation plots
            fig = plt.figure(122)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyse_physio_3")

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
            fig.tight_layout()

            # display FFT plots
            fig = plt.figure(123)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyse_physio_4")

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
            fig.tight_layout()

        # done

    def correct_analyse_and_reject(self, peak_range, moving_Naverages, params_min, params_max, relativeToMean=False, auto_adjust_bounds=False, display=True):
        """
        Analyse peak in each average in terms intensity, linewidth, chemical shift and phase and reject data if one of these parameters goes out of the min / max bounds. Usefull to understand what the hell went wrong during your acquisition when you have the raw data (TWIX) and to try to improve things a little...

        Parameters
        ----------
        peak_range : array [2]
            Range in ppm used to analyse peak phase when no reference signal is specified
        moving_Naverages : int
            Number of averages to perform when using moving average, need to be an odd number
        params_min : array [4]
            Minimum values allowed for peak analysis parameters in order to keep the data: amplitude (%), linewidth (Hz), chemical shift (ppm) and phase (rd)
        params_max : array [4]
            Maximum values allowed for peak analysis parameters in order to keep the data: amplitude (%), linewidth (Hz), chemical shift (ppm) and phase (rd)
        relativeToMean : boolean
            The variation calculation needs a reference. Use the average over the whole data (True) or the first shot (False)
        auto_adjust_bounds : boolean
            Try to adjust the peak parameter thresholds automatically (LW only for now)
        display : boolean
            Display correction process (True) or not (False)

        Returns
        -------
        s_cor : MRSData2 numpy array [averages,timepoints]
           Data remaining after data rejection stored in a MRSData2 object.
        """
        # hello
        cprint(">> mrs.reco.MRSData2.correct_analyse_and_reject:", 'green')

        # init
        s = self.copy()
        if(s.shape[0] == 1):
            print(" > single-shot signal, nothing to analyse!")
            return(s)

        ppm = s.frequency_axis_ppm()

        # build moving averaged data
        s_ma = self._build_moving_average_data(moving_Naverages)

        # perform peak analysis
        params_trace, params_trace_rel = s_ma._analyse_peak(peak_range, relativeToMean)

        # choose if absolute or relative will be checked
        params_trace_check = params_trace_rel * 0.0
        # amplitude: relative in %
        params_trace_check[:, 0] = params_trace_rel[:, 0]
        # linewidth: absolute in Hz
        params_trace_check[:, 1] = params_trace[:, 1]
        # frequency: relative in ppm
        params_trace_check[:, 2] = params_trace_rel[:, 2]
        # phase: absolute in rad
        params_trace_check[:, 3] = params_trace[:, 3]

        # choose if absolute or relative will be displayed
        params_trace_disp = params_trace_rel * 0.0
        # amplitude: relative in %
        params_trace_disp[:, 0] = params_trace_rel[:, 0]
        # linewidth: absolute in Hz
        params_trace_disp[:, 1] = params_trace[:, 1]
        # frequency: absolute in ppm
        params_trace_disp[:, 2] = params_trace[:, 2]
        # phase: absolute in rad
        params_trace_disp[:, 3] = params_trace[:, 3]

        # our time scale
        t_ma = np.linspace(0, self.tr * s.shape[0], s_ma.shape[0]) / 1000.0  # s

        # stats
        print(">> peak analysis -> means + / - std. deviations")
        print(" > Rel. peak amplitude = %.2f + / - %.2f %%" % (params_trace_disp[:, 0].mean(), params_trace_disp[:, 0].std()))
        print(" > Abs. linewidth = %.1f + / - %.1f Hz (%.3f + / - %.3f ppm)" % (params_trace_disp[:, 1].mean(), params_trace_disp[:, 1].std(), (params_trace_disp[:, 1] / s_ma.f0).mean(), (params_trace_disp[:, 1] / s_ma.f0).std()))
        print(" > Abs. frequency = %.2f + / - %.2f ppm ( +  / - %.1f Hz)" % (params_trace_disp[:, 2].mean(), params_trace_disp[:, 2].std(), (params_trace_disp[:, 2] * s_ma.f0).std()))
        print(" > Abs. phase = %.2f + / - %.2f rad" % (params_trace_disp[:, 3].mean(), params_trace_disp[:, 3].std()))

        if(auto_adjust_bounds):
            # auto mode

            # copy min / max params
            params_min_auto = params_min.copy()
            params_max_auto = params_max.copy()

            # iterate between max and min for linewidth, and test the resulting data
            print(">> adjusting linewidth bounds ... ", end="", flush=True)
            p = 1
            max_lw_range = np.arange(int(params_min_auto[p]), int(params_max_auto[p]))
            if(len(max_lw_range) == 0):
                raise Exception(" > You chose to automatically adjust the linewidth threshold but you did not leave any range for it %.0f-%.0fHz!" % (int(params_min_auto[p]), int(params_max_auto[p])))

            self._print_progress_bar(0, max_lw_range.shape[0])
            test_snr_list = np.zeros(max_lw_range.shape)
            test_lw_list = np.zeros(max_lw_range.shape)
            test_nrej_list = np.zeros(max_lw_range.shape)
            for (iLW, this_max_lw) in enumerate(max_lw_range):
                # for each max lw, reject data and see
                this_mask_reject_data = np.full([s_ma.shape[0], 4], False)
                params_max_auto[1] = this_max_lw
                for a in range(0, s_ma.shape[0]):
                    for p in range(4):
                        if(params_trace_check[a, p] < params_min_auto[p]):
                            this_mask_reject_data[a, p] = True
                        if(params_trace_check[a, p] > params_max_auto[p]):
                            this_mask_reject_data[a, p] = True

                # reject data now
                this_mask_reject_data_sumup = (this_mask_reject_data.sum(axis=1) > 0)
                this_s_cor = s[(this_mask_reject_data_sumup == False), :]

                # analyse snr / lw and number of rejections
                if(this_mask_reject_data_sumup.sum() < s_ma.shape[0]):
                    test_snr_list[iLW] = this_s_cor._correct_realign()._correct_average()._correct_apodization().analyse_snr(peak_range, [-1, 0], '', False, False, False, [1, 6], True)
                    test_lw_list[iLW] = this_s_cor._correct_realign()._correct_average()._correct_apodization().analyse_linewidth(peak_range, '', False, False, [1, 6], True)
                    test_nrej_list[iLW] = this_mask_reject_data_sumup.sum()

                # progression
                self._print_progress_bar(iLW)

            print(" done.")

            # plot SNR / LW combinaisons
            if(display):
                fig = plt.figure(130)
                fig.clf()
                axs = fig.subplots(1, 2, sharex='all')
                ax2 = axs[0].twinx()
                ax3 = axs[1].twinx()
                fig.canvas.set_window_title(
                    "mrs.reco.MRSData2.correct_analyse_and_reject (auto)")

                axs[0].plot(max_lw_range, test_snr_list, 'rx-', label='SNR')
                axs[0].set_xlabel('Max allowed linewidth (Hz)')
                axs[0].set_ylabel('Estimated SNR (u.a)')
                axs[0].grid('on')
                axs[0].legend(loc='lower left')

                ax2.plot(max_lw_range, test_lw_list, 'bx-', label='Linewidth')
                ax2.set_ylabel('Estimated linewidth (Hz)')
                ax2.legend(loc='lower right')

                axs[1].plot(max_lw_range, test_nrej_list, 'ko-',
                            label='Number of scans rejected')
                axs[1].set_xlabel('Max allowed linewidth (Hz)')
                axs[1].set_ylabel('Number of scans rejected')
                axs[1].grid('on')

                ax3.plot(max_lw_range, test_nrej_list / s_ma.shape[0] * 100, 'ko-', label='Total percentage of scans rejected')
                ax3.set_ylabel('Estimated linewidth (Hz)')
                ax3.set_ylabel('Rejection percentage (%)')

                fig.tight_layout()

        # for each average, check if peak parameters are in the min / max bounds
        print(">> rejecting data ... ", end="", flush=True)
        mask_reject_data = np.full([s_ma.shape[0], 4], False)
        self._print_progress_bar(0, s_ma.shape[0])
        for a in range(0, s_ma.shape[0]):
            for p in range(4):
                if(params_trace_check[a, p] < params_min[p]):
                    mask_reject_data[a, p] = True
                if(params_trace_check[a, p] > params_max[p]):
                    mask_reject_data[a, p] = True

            self._print_progress_bar(a)

        print(" done.")

        # stats regarding data rejection, how many, for what reasons, overall percentage
        print(">> data rejection -> summary")
        print(">> number of averages rejected because of...")
        print(" > amplitude = %d" % mask_reject_data[:, 0].sum())
        print(" > linewidth = %d" % mask_reject_data[:, 1].sum())
        print(" > frequency = %d" % mask_reject_data[:, 2].sum())
        print(" > phase = %d" % mask_reject_data[:, 3].sum())

        # actually reject data now
        mask_reject_data_sumup = (mask_reject_data.sum(axis=1) > 0)
        s_cor = s[(mask_reject_data_sumup == False), :]

        print(" > TOTAL data rejection = %d / %d (%.0f%%)" % (mask_reject_data_sumup.sum(), s_ma.shape[0], (mask_reject_data_sumup.sum() / s_ma.shape[0] * 100)))

        # final display
        if(display):

            fig = plt.figure(131)
            fig.clf()
            axs = fig.subplots(2, 3, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_analyse_and_reject")

            k = 0
            for ix in range(2):
                for iy in range(2):
                    # original data
                    axs[ix, iy].plot(t_ma, params_trace_check[:, k], 'k-x', linewidth=1)
                    # rejected data
                    t_ma_rej = t_ma[mask_reject_data[:, k]]
                    this_params_trace_rej = params_trace_check[mask_reject_data[:, k], k]
                    axs[ix, iy].plot(t_ma_rej, this_params_trace_rej, 'ro', linewidth=1)
                    axs[ix, iy].plot(t_ma, params_min[k] * np.ones(t_ma.shape), '-r', linewidth=1)
                    axs[ix, iy].plot(t_ma, params_max[k] * np.ones(t_ma.shape), '-r', linewidth=1)
                    k = k + 1

            axs[0, 0].set_ylabel('Rel. amplitude change (%)')
            axs[0, 0].grid('on')
            axs[0, 0].set_title("Rel. amplitude = %.2f + / - %.2f %%" % (params_trace_disp[:, 0].mean(), params_trace_disp[:, 0].std()))

            axs[0, 1].set_ylabel('Abs. linewidth (Hz)')
            axs[0, 1].grid('on')
            axs[0, 1].set_title("Abs. linewidth = %.1f + / - %.1f Hz (%.3f + / - %.3f ppm)" % (params_trace_disp[:, 1].mean(), params_trace_disp[:, 1].std(), (params_trace_disp[:, 1] / s_ma.f0).mean(), (params_trace_disp[:, 1] / s_ma.f0).std()))

            axs[1, 0].set_xlabel('Acq. time (s)')
            axs[1, 0].set_ylabel('Abs. frequency (ppm)')
            axs[1, 0].grid('on')
            axs[1, 0].set_title("Abs. frequency = %.2f + / - %.2f ppm ( +  / - %.1f Hz)" % (params_trace_disp[:, 2].mean(), params_trace_disp[:, 2].std(), (params_trace_disp[:, 2] * s_ma.f0).std()))

            axs[1, 1].set_xlabel('Acq. time (s)')
            axs[1, 1].set_ylabel('Abs. phase shift (rd)')
            axs[1, 1].grid('on')
            axs[1, 1].set_title("Abs. phase = %.2f + / - %.2f rad" % (params_trace_disp[:, 3].mean(), params_trace_disp[:, 3].std()))

            # nice plot showing all raw data
            plt.subplot(1, 3, 3)
            ppm = self.frequency_axis_ppm()
            ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
            ystep = np.max(np.mean(s_ma.spectrum().real, axis=0))
            ystep = np.power(10, 1 + (np.floor(np.log10(ystep))))

            ampfactor = 4
            for k in range(s_ma.shape[0]):
                if(mask_reject_data_sumup[k]):
                    plt.plot(ppm, s_ma[k, :].spectrum().real * ampfactor + ystep * k, 'r-', linewidth=1)
                else:
                    plt.plot(ppm, s_ma[k, :].spectrum().real * ampfactor + ystep * k, 'g-', linewidth=1)

                sf_masked = np.real(s_ma[k, :].spectrum())
                sf_masked[ippm_peak_range] = 0
                ippm_peak = np.argmax(sf_masked)
                if(ippm_peak == 0):
                    raise Exception(" > no peak found in specified ppm range or badly phased data!")
                amp_peak = sf_masked[ippm_peak]
                ippm_half_peak = np.where(sf_masked > amp_peak / 2)
                plt.plot(ppm[ippm_half_peak], sf_masked[ippm_half_peak].spectrum().real * ampfactor + ystep * k, 'k-', linewidth=1)

            plt.xlim(peak_range[1], peak_range[0])
            plt.xlabel('chemical shift (ppm)')
            plt.ylabel('individual spectra')
            plt.grid('on')

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fig.tight_layout()

            # plt.pause(0.1)

        # wait, are we removing all data ???
        if(mask_reject_data_sumup.sum() == s.shape[0]):
            raise Exception(" > All data is rejected! You need to readjust your rejection bounds...")

        return(s_cor)

    def correct_realign(self, peak_range, moving_Naverages, display=True, display_range=[1, 6]):
        """Kind of wrapper method for the method just below."""
        return(self._correct_realign(peak_range, moving_Naverages, display, display_range, False))

    def _correct_realign(self, peak_range=[4, 6], moving_Naverages=1, display=False, display_range=[1, 6], beSilent=True):
        """
        Realign each signal of interest in frequency by taking as a reference the first spectra in absolute mode.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to analyse peak phase
        moving_Naverages : int
            Number of averages to perform when using moving average, need to be an odd number
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display
        beSilent : boolean
            No outpur in console (True)

        Returns
        -------
        s_realigned : MRSData2 numpy array [averages,timepoints]
            Resulting frequency realigned data stored in a MRSData2 object
        """
        # hello
        if(not beSilent):
            cprint(">> mrs.reco.MRSData2.correct_realign:", 'green')

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()
        s_realigned = s.copy()

        if(s.shape[0] == 1):
            if(not beSilent):
                print(" > single-shot signal, cannot realign!")
        else:
            # build moving averaged data
            s_ma = self._build_moving_average_data(moving_Naverages, beSilent)

            # find peak in average spectrum absolute mode
            s_avg = np.mean(s, axis=0)
            sf_avg_abs = np.abs(s_avg.spectrum())
            ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
            sf_masked = sf_avg_abs.copy()
            sf_masked[ippm_peak_range] = 0
            ippm_peak_avg = np.argmax(sf_masked, axis=0)
            if(ippm_peak_avg == 0):
                raise Exception(" > no peak found in specified ppm range or badly phased data!")
            ppm_peak_avg = ppm[ippm_peak_avg]
            if(not beSilent):
                print(" > measuring peak properties at %0.2fppm!" % ppm_peak_avg)

            # for each average in moving averaged data
            s_realigned_ma = s_ma.copy()
            df_trace = np.zeros(s_ma.shape[0])
            if(not beSilent):
                print(" > realigning... ", end="", flush=True)
            if(not beSilent):
                self._print_progress_bar(0, s_ma.shape[0])
            for a in range(0, s_ma.shape[0]):

                # measure shift on moving average data
                sf_masked = np.abs(s_ma[a, :].spectrum())
                sf_masked[ippm_peak_range] = 0
                ippm_peak = np.argmax(sf_masked)
                if(ippm_peak == 0):
                    raise Exception(" > no peak found in specified ppm range or badly phased data!")
                ppm_peak = ppm[ippm_peak]

                # estimate frequency shift in Hz compared to average spectrum
                dppm = -(ppm_peak_avg - ppm_peak)
                df_trace[a] = dppm * s_ma.f0

                # correct moving averaged data
                s_realigned_ma[a, :] = s_ma[a, :].adjust_frequency(df_trace[a])

                # correct original data
                s_realigned[a, :] = s[a, :].adjust_frequency(df_trace[a])

                if(not beSilent):
                    self._print_progress_bar(a)
            if(not beSilent):
                print(" done.")

            # final display
            if(display):

                fig = plt.figure(140)
                fig.clf()
                axs = fig.subplots(2, 3, sharex='all', sharey='all')
                fig.canvas.set_window_title("mrs.reco.MRSData2.correct_realign")

                # display original averaged spectrum
                axs[0, 0].plot(ppm, sf_avg_abs, 'k-', linewidth=1)
                axs[0, 0].set_xlim(display_range[1], display_range[0])
                axs[0, 0].set_ylabel('averaged')
                axs[0, 0].grid('on')
                # add peak position
                axs[0, 0].plot(ppm_peak_avg, sf_avg_abs[ippm_peak_avg], 'ro')

                # display original data
                axs[0, 1].plot(ppm, np.abs(s_ma.spectrum().transpose()), 'k-', linewidth=1)
                axs[0, 1].set_xlim(display_range[1], display_range[0])
                axs[0, 1].set_ylabel('original')
                axs[0, 1].grid('on')

                # display corrected spectrum
                axs[1, 0].plot(s_realigned_ma.frequency_axis_ppm(), np.abs(s_realigned_ma.spectrum().transpose()), 'k-', linewidth=1)
                axs[1, 0].set_xlim(display_range[1], display_range[0])
                axs[1, 0].set_xlabel('chemical shift (ppm)')
                axs[1, 0].set_ylabel('corrected')
                axs[1, 0].grid('on')

                # display corrected averaged spectrum
                axs[1, 1].plot(s_realigned_ma.frequency_axis_ppm(), np.abs(np.mean(s_realigned_ma, axis=0).spectrum().transpose()), 'k-', linewidth=1)
                axs[1, 1].set_xlim(display_range[1], display_range[0])
                axs[1, 1].set_xlabel('chemical shift (ppm)')
                axs[1, 1].set_ylabel('averaged & corrected')
                axs[1, 1].grid('on')

                plt.subplot(1, 3, 3)
                plt.plot(df_trace, np.arange(s_ma.shape[0]), 'k-x', linewidth=1)
                plt.xlabel('estimated frequency shift (Hz)')
                plt.ylabel('average index')
                plt.grid('on')

                fig.tight_layout()
                # plt.pause(0.1)

        return(self.inherit(s_realigned))

    def correct_average(self, na=None, display=True, display_range=[1, 6]):
        """Kind of wrapper method for the method just below."""
        return(self._correct_average(na, display, display_range, False))

    def _correct_average(self, na=None, display=False, display_range=[1, 6], beSilent=True):
        """
        Average all averages data into one 1D MRS signal.

        Parameters
        ----------
        na : int
            Number of signals to average
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display
        beSilent : boolean
            No output in console (True)

        Returns
        -------
        s_mean : MRSData2 numpy array [timepoints]
            Resulting frequency realigned data stored in a MRSData2 object
        """
        # hello
        if(not beSilent):
            cprint(">> mrs.reco.MRSData2.correct_average", 'green')

        # init
        s = self.copy()

        if(s.shape[0] == 1):
            if(not beSilent):
                print(" > single-shot signal, nothing to average!")
            s_mean = np.mean(s, axis=0)
            if(not beSilent):
                print(" > reshaped to a " + str(s_mean.shape) + " vector")
        else:
            if(not beSilent):
                print(" > averaging data...", end="", flush=True)

            # check parameters
            if not display:
                display = 0
            if not any(display_range):
                display_range = np.array([0, 5])

            if(na is not None):
                if(na == 1):
                    s_mean = s[0, :]
                else:
                    s_mean = np.mean(s[0:(na - 1), :], axis=0)

                if(not beSilent):
                    print("only " + str(na) + "...", end="", flush=True)
            else:
                s_mean = np.mean(s, axis=0)

            # keep average dimension and set to 1
            # s_mean=s_mean.reshape((1,) + s_mean.shape)

            if(display):
                ppm = s.frequency_axis_ppm()
                ppm_mean = s_mean.frequency_axis_ppm()

                fig = plt.figure(150)
                fig.clf()
                axs = fig.subplots(2, 1, sharex='all', sharey='all')
                fig.canvas.set_window_title("mrs.reco.MRSData2.correct_average")

                axs[0].plot(ppm, s.spectrum().real.transpose(), 'k-', linewidth=1)
                axs[0].set_xlim(display_range[1], display_range[0])
                axs[0].set_xlabel('chemical shift (ppm)')
                axs[0].set_ylabel('all spectra')
                axs[0].grid('on')

                axs[1].plot(ppm_mean, s_mean.spectrum().real.transpose(), 'k-', linewidth=1)
                axs[1].set_xlim(display_range[1], display_range[0])
                axs[1].set_xlabel('chemical shift (ppm)')
                axs[1].set_ylabel('averaged spectrum')
                axs[1].grid('on')

                fig.tight_layout()
                # plt.pause(0.1)

            if(not beSilent):
                print("done.")

        return(self.inherit(s_mean))

    def correct_apodization(self, apo_factor, nPoints_final=4096, display=True, display_range=[1, 6]):
        """Kind of wrapper method for the method just below."""
        return(self._correct_apodization(apo_factor, nPoints_final, display, display_range, False))

    def _correct_apodization(self, apo_factor=1.0, nPoints_final=4096, display=False, display_range=[1, 6], beSilent=True):
        """
        Apodize signal using an exponential window adjusted by a linewidth parameter in Hz. Works only for a 1D MRSData2 object.

        Parameters
        ----------
        apo_factor : float
            Apodization factor in Hz
        nPoints_final : int
            Final number of points (crop)
        display : boolean
            Display correction process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display
        beSilent : boolean
            No output in console (True)

        Returns
        -------
        s_apo : MRSData2 numpy array [timepoints]
            Resulting apodized data stored in a MRSData2 object
        """
        # hello
        if(not beSilent):
            cprint(">> mrs.reco.MRSData2.correct_apodization:", 'green')

        # init
        if(not beSilent):
            print(" > apodizing data", end="", flush=True)
        s = self.copy()
        t = s.time_axis()
        w_apo = np.exp(-apo_factor * t)
        s_apo = s * w_apo
        # crop
        if(nPoints_final < s_apo.shape[0]):
            if(not beSilent):
                print(", cropping data, updating sequence npts ", end="", flush=True)
            s_apo = s_apo[0:nPoints_final]
            self.sequence.npts = nPoints_final
            self.sequence._ready = False

        if(display):
            t_apo = s_apo.time_axis()
            ppm = s.frequency_axis_ppm()
            ppm_apo = s_apo.frequency_axis_ppm()

            fig = plt.figure(160)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='row', sharey='row')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_apodization")

            axs[0, 0].plot(t, np.abs(s), 'k-', linewidth=1, label='fid')
            axs[0, 0].plot(t, w_apo * np.abs(s.max()), 'r-', linewidth=1, label='apodization window')
            axs[0, 0].set_xlabel('time (s)')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(t_apo, np.abs(s_apo), 'k-', linewidth=1)
            axs[0, 1].set_xlabel('time (s)')
            axs[0, 1].set_ylabel('apodized')
            axs[0, 1].grid('on')

            axs[1, 0].plot(ppm, s.spectrum().real, 'k-', linewidth=1)
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('original spectrum')
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].grid('on')

            axs[1, 1].plot(ppm_apo, s_apo.spectrum().real, 'k-', linewidth=1)
            axs[1, 1].set_xlabel("chemical shift (ppm)")
            axs[1, 1].set_ylabel('apodized spectrum')
            axs[1, 1].set_xlim(display_range[1], display_range[0])
            axs[1, 1].grid('on')

            fig.tight_layout()
            # plt.pause(0.1)

        if(not beSilent):
            print("...done.")

        return(self.inherit(s_apo))

    def correct_water_removal(self, hsvd_nComponents, hsvd_range, display=True, display_range=[1, 6]):
        """
        Remove any water residual peak(s) within a ppm range using HSVD.

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
        # hello
        cprint(">> mrs.reco.MRSData2.correct_water_removal:", 'green')

        # init
        print(" > removing residual water peak with HSVD... ", end="", flush=True)
        s = self.copy()
        ppm = s.frequency_axis_ppm()
        self._print_progress_bar(0, 5)

        # estimate HSVD components
        components = suspect.processing.water_suppression.hsvd(s, hsvd_nComponents)
        self._print_progress_bar(1)

        # filter them by keeping the ones contributing to the residual water peak and its sidebands
        water_components = [component for component in components if ((4.7 - component["frequency"] / self.f0) > hsvd_range[0] and (4.7 - component["frequency"] / self.f0) < hsvd_range[1])]
        self._print_progress_bar(2)

        # reconstruct the estimated water peak
        hsvd_fid = suspect.processing.water_suppression.construct_fid(water_components, s.time_axis())
        self._print_progress_bar(3)

        # rebuild object
        hsvd_fid = s.inherit(hsvd_fid)
        self._print_progress_bar(4)

        # and substract it from the fid
        s_water_removed = s - hsvd_fid
        self._print_progress_bar(5)

        # display this over the data
        if(display):
            fig = plt.figure(170)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_water_removal")

            # original spectrum
            axs[0].plot(ppm, s.spectrum().real, 'k-', linewidth=1, label='original data (real part)')
            axs[0].plot(ppm, hsvd_fid.spectrum().real, 'r-', linewidth=1, label='estimated HSVD water peak')
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('original spectrum')
            axs[0].grid('on')
            axs[0].legend()

            # water removed spectrum
            axs[1].plot(ppm, s_water_removed.spectrum().real, 'k-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('water removed spectrum')
            axs[1].grid('on')
            axs[1].legend()

            fig.tight_layout()
            # plt.pause(0.1)

        print(" done.")

        return(self.inherit(s_water_removed))

    def correct_freqshift(self, peak_range=[4, 5], peak_real_ppm=4.7, display=True, display_range=[1, 6]):
        """Kind of wrapper method for the method just below."""
        return(self._correct_freqshift(peak_range, peak_real_ppm, display, display_range, False))

    def _correct_freqshift(self, peak_range=[4, 5], peak_real_ppm=4.7, display=False, display_range=[1, 6], beSilent=True):
        """
        Shift the spectrum in frequency in order to get the right peaks at the right chemical shifts. Peak-picking is done in magnitude mode to reduce sensitivity to phase shit.

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
        # hello
        if(not beSilent):
            cprint(">> mrs.reco.MRSData2.correct_freqshift:", 'green')

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()

        # find maximum peak in range and its chemical shift
        ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
        sf_abs_masked = np.abs(s.spectrum())
        sf_abs_masked[ippm_peak_range] = 0
        ippm_peak = np.argmax(sf_abs_masked)
        if(ippm_peak == 0):
            raise Exception(" > no peak found in specified ppm range or badly phased data!")
        ppm_peak = ppm[ippm_peak]
        if(not beSilent):
            print(" > peak detected at %0.2fppm -> %0.2fppm!" % (ppm_peak, peak_real_ppm))

        if(not beSilent):
            print(" > frequency shifting data...", end="", flush=True)
        # estimate frequency shift in Hz
        dppm = (peak_real_ppm - ppm_peak)
        df = dppm * s.f0
        s_shifted = s.adjust_frequency(-df)

        if(display):
            fig = plt.figure(180)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_freqshift")

            axs[0].plot(ppm, s.spectrum().real, 'k-', linewidth=1)
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('original')
            axs[0].grid('on')

            axs[1].plot(ppm, s_shifted.spectrum().real, 'k-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('shifted')
            axs[1].grid('on')

            fig.tight_layout()
            # plt.pause(0.1)

        if(not beSilent):
            print("done.")

        return(self.inherit(s_shifted))

    def analyse_snr(self, peak_range, noise_range, lbl, area_integrate=False, magnitude_mode=False, display=True, display_range=[1, 6], beSilent=False):
        """
        Estimate the SNR of a peak in the spectrum ; chemical shift ranges for the peak and the noise regions are specified by the user. !Works only for a 1D MRSData2 object.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to find a peak of interest
        noise_range : list [2]
            Range in ppm used to estimate noise
        lbl : string
            Plot label to specify
        area_integrate : boolean
            Integrate spectrum within ppm range instead of peak-picking
        magnitude_mode : boolean
            Analyse signal in magnitude mode (True) or the real part (False)
        display : boolean
            Display process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display
        beSilent : boolean
            No outpur in console (True)

        Returns
        -------
        snr : float
            Resulting SNR value
        """
        # hello
        if(not beSilent):
            cprint(">> mrs.reco.MRSData2.analyse_snr:", 'green')

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()

        # find maximum peak in range and its chemical shift
        sf = s.spectrum()
        if(magnitude_mode):
            sf_analyse = np.abs(sf)
            if(not beSilent):
                print(" > going to analyse the MAGNITUDE spectrum ", end="", flush=True)
        else:
            sf_analyse = np.real(sf)
            if(not beSilent):
                print(" > going to analyse the REAL spectrum ", end="", flush=True)

        if(area_integrate):
            # integrate area underneath the peak
            ippm_peak_range = (peak_range[0] < ppm) & (ppm < peak_range[1])
            if(not beSilent):
                print("by integrating the area within [%0.2f-%0.2f] ppm!" % (peak_range[0], peak_range[1]))
            sf_analyse2 = sf_analyse[ippm_peak_range]
            ppm2 = ppm[ippm_peak_range]
            sf_analyse2 = sf_analyse2[::-1]
            ppm2 = ppm2[::-1]
            snr_signal = np.trapz(sf_analyse2, ppm2)
        else:
            ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
            sf_analyse2 = sf_analyse.copy()
            sf_analyse2[ippm_peak_range] = 0
            ippm_peak = np.argmax(sf_analyse2)
            if(ippm_peak == 0):
                Warning(" > no peak found in specified ppm range or badly phased data!")
                return(np.nan)

            ppm_peak = ppm[ippm_peak]
            if(not beSilent):
                print("by measuring the intensity at %0.2fppm!" % ppm_peak)

            snr_signal = sf_analyse[ippm_peak]

        # estimate noise in user specified spectral region
        if(not beSilent):
            print(" > estimating noise from %0.2f to %0.2fppm region!" % (noise_range[0], noise_range[1]))
        ippm_noise_range = (noise_range[0] < ppm) & (ppm < noise_range[1])
        snr_noise = np.std(sf_analyse[ippm_noise_range])

        # that's it
        snr = snr_signal / snr_noise
        if(not beSilent):
            print(" > results for [" + lbl + "] coming...")
        if(not beSilent):
            print(" > S = %f, N = %f, SNR = %0.2f!" % (snr_signal, snr_noise, snr))

        if(display):
            fig = plt.figure(190)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyse_snr")

            axs[0].plot(ppm, np.real(s.spectrum()), 'k-', linewidth=1)
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('real part')
            axs[0].grid('on')

            axs[1].plot(ppm, np.abs(s.spectrum()), 'k-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('magnitude mode')
            axs[1].grid('on')

            if(magnitude_mode):
                ax = axs[1]
            else:
                ax = axs[0]

            # show peak of interest
            if(area_integrate):
                ax.fill_between(ppm[ippm_peak_range], 0, sf_analyse[ippm_peak_range], facecolor='red')
            else:
                ax.plot(ppm[ippm_peak], sf_analyse[ippm_peak], 'ro')

            # show noise region
            ax.plot(ppm[ippm_noise_range], sf_analyse[ippm_noise_range], 'bo')

            fig.tight_layout()
            # plt.pause(0.1)

        return(snr)

    def analyse_linewidth(self, peak_range, lbl, magnitude_mode=False, display=True, display_range=[1, 6], beSilent=False):
        """
        Estimate the linewidth of a peak in the spectrum ; chemical shift ranges for the peak and the noise regions are specified by the user. !Works only for a 1D MRSData2 object.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to find a peak of interest
        lbl : string
            Plot label to specify
        magnitude_mode : boolean
            Analyse signal in magnitude mode (True) or the real part (False)
        display : boolean
            Display process (True) or not (False)
        display_range : list [2]
            Range in ppm used for display
        beSilent : boolean
            No outpur in console (True)

        Returns
        -------
        lw : float
            Linewidth in Hz
        """
        # hello
        if(not beSilent):
            cprint(">> mrs.reco.MRSData2.analyse_linewidth:", 'green')

        # init
        s = self.copy()
        ppm = s.frequency_axis_ppm()

        # find maximum peak in range and its chemical shift
        sf = s.spectrum()
        if(magnitude_mode):
            sf_analyse = np.abs(sf)
            if(not beSilent):
                print(" > going to analyse the MAGNITUDE spectrum ", end="", flush=True)
        else:
            sf_analyse = np.real(sf)
            if(not beSilent):
                print(" > going to analyse the REAL spectrum ", end="", flush=True)

        # find peak in range
        ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
        sf_analyse[ippm_peak_range] = 0
        ippm_peak = np.argmax(sf_analyse)
        if(ippm_peak == 0):
            raise Exception(" > no peak found in specified ppm range or badly phased data!")
        ppm_peak = ppm[ippm_peak]
        if(not beSilent):
            print("and estimate the linewidth of the peak at %0.2fppm!" % ppm_peak)

        # estimate linewidth in Hz
        amp_peak = sf_analyse[ippm_peak]
        ippm_half_peak = np.where(sf_analyse > amp_peak / 2)
        ippm_min = np.min(ippm_half_peak)
        ippm_max = np.max(ippm_half_peak)
        dppm = np.abs(ppm[ippm_max] - ppm[ippm_min])
        lw = dppm * s.f0
        if(not beSilent):
            print(" > results for [" + lbl + "] coming...")
        if(not beSilent):
            print(" > LW = %0.2f Hz!" % lw)

        if(display):
            fig = plt.figure(200)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyse_linewidth")

            axs[0].plot(ppm, np.real(s.spectrum()), 'k-', linewidth=1)
            axs[0].set_xlim(display_range[1], display_range[0])
            axs[0].set_xlabel('chemical shift (ppm)')
            axs[0].set_ylabel('real part')
            axs[0].grid('on')

            axs[1].plot(ppm, np.abs(s.spectrum()), 'k-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('magnitude mode')
            axs[1].grid('on')

            if(magnitude_mode):
                ax = axs[1]
            else:
                ax = axs[0]

            # show noise region
            ax.plot(ppm[ippm_half_peak], sf_analyse[ippm_half_peak], 'r-')

            fig.tight_layout()
            # plt.pause(0.1)

        return(lw)

    def display_spectrum(self, ifig=1, lbl="", display_range=[1, 6], amp_factor=1.0, yoffset=0.0, magnitude_mode=False):
        """
        Display spectrum in figure 'ifig', overlaying if needed.

        Parameters
        ----------
        ifig: int
            The figure index that shoud host the plot
        s : MRSData2 numpy array [timepoints]
            MRS data to display
        lbl : string
            Plot label to specify
        display_range : list [2]
            Range in ppm used for display
        amp_factor : float
            Spectrum intensity amplification factor
        yoffset : float
            Spectrum intensity offset
        magnitude_mode : boolean
            Displays in magnitude mode (True) or the real part (False)

        Returns
        -------
        fig : matplotlib.figure
            Resulting matplotlib figure
        """
        # hello
        cprint(">> mrs.reco.MRSData2.display_spectrum:", 'green')

        # init
        s = self.copy()
        print(" > displaying stuff!")

        plt.figure(ifig).canvas.set_window_title("mrs.reco.MRSData2.display_spectrum")
        if(magnitude_mode):
            plt.plot(s.frequency_axis_ppm(), np.abs(s.spectrum()) * amp_factor + yoffset, linewidth=1, label=lbl)
        else:
            plt.plot(s.frequency_axis_ppm(), s.spectrum().real * amp_factor + yoffset, linewidth=1, label=lbl)

        # add ytick if offset
        if(yoffset != 0):
            yt = plt.yticks()
            yt2 = np.hstack((yt[0], yoffset))
            yt3 = np.sort(yt2)
            plt.yticks(yt3)

        if any(display_range):
            plt.xlim(display_range[1], display_range[0])

        plt.xlabel('chemical shift (ppm)')
        plt.ylabel('spectrum')
        plt.grid('on')
        plt.legend()

        plt.tight_layout()

        return(plt)

    def save2mat(self, mat_filepath):
        """
        Save the numpy array content to a MATLAB mat file.

        Parameters
        ----------
        mat_filepath: string
            Full absolute file path pointing to mat file
        """
        # hello
        cprint(">> mrs.MRSData2.save2mat:", 'green')

        print(" > saving MRS signal to " + mat_filepath + "...")
        sio.savemat(mat_filepath, {'MRSdata': self})


class pipeline:
    """The pipeline class is used to store all the reconstruction parameters needed for a specific bunch of acquired signals. Once the parameters are set, the pipeline can be run using one of the methods."""

    def __init__(self):

        # --- water-suppressed data ---
        # list of file paths pointaing to datasets we wish to process
        # this can be either one big string with filepaths separated by line breaks (that will be parsed) OR simply a list
        self.data_filepaths = ""
        # parsed list  of dataset filepaths
        self._data_filepaths_list = []
        # list of data signals
        self._data_list = []

        # --- non water-suppressed data ---
        # list of reference (non water suppressed) file paths pointaing to datasets we wish to process
        # this can be either one big string with filepaths separated by line breaks (that will be parsed) OR simply a list
        self.data_ref_filepaths = ""
        # parsed list  of dataset filepaths
        self._data_ref_filepaths_list = []
        # list of reference data signals
        self._data_ref_list = []

        # to avoid any data reshaping misunderstanding, the number of coil channels
        self.data_coil_nChannels = 8

        # --- physio data ---
        # list of filepaths pointing to physio files
        self.data_physio_filepaths = ""
        # parsed list of physio filepaths
        self._data_physio_list = []

        # set reference ppm
        self.ppm0 = 4.7

        # option to process only a set of datasets: lsit of indexes
        self.data_process_only_this_data_index = []
        # by default, process each data/ref.data separatly. If this option is set to true, all data will be concatenated
        self.data_concatenate = False

        # --- FID modulus process ---
        self.fid_modulus = False

        # --- automatic rephasing ---
        self.phase_enable = True
        # weak water suppression was used, meaning we can use the 1st point of the FID (mainly water)
        self.phase_weak_ws_mode = False
        # super high SNR? cool, we can rephase each individual spectrum from each channel separatly (by default, spectra are averaged by channel)
        self.phase_high_snr_mode = False
        # if peak phasing, order of rephasing: 0 means 0th order, 1 means 0th + 1st order (in time) phasing
        self.phase_order = 0
        # add an additional 0th order phase (rd)
        self.phase_offset = 0
        # ppm range to look fo peak used to estimate phase
        self.phase_POI_range_ppm = [1.5, 2.5]
        # display all this process to check what the hell is going on
        self.phase_display = False

        # --- channel recombination ---
        self.recombine_enable = True
        # use non water-suppressed data to recombine and rephase channels
        self.recombine_use_data_ref = True
        # should we rephase (0th order) data while recombining?
        self.recombine_phasing = True
        # boolean mask to switch on/off some Rx channels
        self.recombine_weights = [True]

        # --- zero-filling ---
        self.zerofill_enable = True
        # number of signal points after zf
        self.zerofill_npts = 8192 * 2
        # display all this process to check what the hell is going on
        self.zerofill_display = True

        # --- analyse physio signal ---
        self.analyse_physio_enable = False
        # ppm range to look for a peak to analyse
        self.analyse_physio_POI_range_ppm = [4, 5]
        # time range in (ms) to look around timestamp for the best correlation physio/MRS
        self.analyse_physio_delta_time_ms = 1000.0
        # display all this process to check what the hell is going on
        self.analyse_physio_display = True

        # --- automatic data rejection based on criterias ---
        self.analyse_and_reject_enable = False
        # ppm range to look for a peak to analyse
        self.analyse_and_reject_POI_range_ppm = [4.5, 5]
        # size of moving average window
        self.analyse_and_reject_moving_averages = 1
        # lower parameter bounds for amplitude changes (%), linewidth (Hz), chemical shift changes (ppm) and phase (rd)
        self.analyse_and_reject_min = [-100, 0, -0.5, -3.14]
        # upper parameter bounds for amplitude changes (%), linewidth (Hz), chemical shift changes (ppm) and phase (rd)
        self.analyse_and_reject_max = [+100, 200, +0.5, +3.14]
        # for relative parameters such as amplitude and chemical shift changes, relative to average over whole data (True) or first point (False)?
        self.analyse_and_reject_relative_mean = True
        # automatic linewidth rejection criteria based on SNR/linewidth evolution
        self.analyse_and_reject_auto = False
        # display all this process to check what the hell is going on
        self.analyse_and_reject_display = True

        # --- automatic data frequency realignment ---
        self.realign_enable = True
        # ppm range to look for a peak to analyse
        self.realign_POI_range_ppm = [4.5, 5]
        # size of moving average window
        self.realign_moving_averages = 1
        # display all this process to check what the hell is going on
        self.realign_display = True

        # --- data averaging ---
        self.average_enable = True
        # number of averages to mean
        self.average_na = None
        # display all this process to check what the hell is going on
        self.average_display = True

        # --- data apodization an crop ---
        self.apodize_enable = True
        # exponential damping factor for apodization (Hz)
        self.apodize_damping_hz = 5
        # final number of signal points after crop
        self.apodize_npts = 4096
        # display all this process to check what the hell is going on
        self.apodize_display = True

        # --- water post-acquisition removal ---
        self.remove_water_enable = False
        # number of components when running HSVD
        self.remove_water_hsvd_components = 5
        # ppm range where all components will be remove
        self.remove_water_hsvd_range = [4.6, 4.8]
        # display all this process to check what the hell is going on
        self.remove_water_display = True

        # --- spectrum chemical shift calibration ---
        self.calibrate_enable = True
        # ppm range to look for the peak of interest (NAA by default)
        self.calibrate_POI_range_ppm = [1.8, 2.3]
        # real ppm value for this peak
        self.calibrate_POI_true_ppm = 2.008  # ppm
        # display all this process to check what the hell is going on
        self.calibrate_display = True

        # --- spectrum final display ---
        self.display_enable = True
        # figure index
        self.display_fig_index = 1
        # ppm range used for display
        self.display_range_ppm = [1, 6]  # ppm
        # spectrum intensity scaling factor
        self.display_amp_factor_list = []
        # spectrum y offset parameter
        self.display_offset = 0
        # display spectrum in magnitude mode?
        self.display_magnitude_mode = False
        # display legend captions
        # this can be either one big string with captions separated by line breaks (that will be parsed) OR simply a list
        self.display_legends = ""
        # internal list of legend captions
        self._display_legends_list = []

        # --- SNR analysis ---
        self.analyse_snr_enable = True
        # ppm range to look for a peak to analyse
        self.analyse_snr_s_range_ppm = [1.8, 2.1]  # ppm
        # should we integrate the peak (True) or just get maximum real intensity (False)?
        self.analyse_snr_area_integrate = False
        # ppm range to look for pute noise
        self.analyse_snr_n_range_ppm = [-1, 0]  # ppm
        # should we look at the magnitude (True) or real (False) spectrum for signal estimation?
        self.analyse_snr_magnitude_mode = False
        # display all this process to check what the hell is going on
        self.analyse_snr_display = True

        # --- SNR analysis during signal processing ---
        # SNR will be estimated by max of real peak intensity / standard deviation of noise
        self.analyse_snr_evol_enable = True
        # list of measured SNR
        self.analyse_snr_evol_list = []
        # list of captions
        self.analyse_snr_evol_list_labels = []
        # list of SNR measured per processed signals
        self.analyse_snr_final_list = []

        self.analyse_linewidth_enable = True
        # ppm range to look for a peak to analyse
        self.analyse_linewidth_range_ppm = [4.5, 5]  # ppm
        # should we look at the magnitude (True) or real (False) spectrum for linewidth estimation?
        self.analyse_linewidth_magnitude_mode = False
        # display all this process to check what the hell is going on
        self.analyse_linewidth_display = True

        # --- Linewidth analysis during signal processing ---
        self.analyse_linewidth_evol_enable = True
        # list of measured linewidths
        self.analyse_linewidth_evol_list = []
        # list of captions
        self.analyse_linewidth_evol_list_labels = []
        # list of linewidhs measured per processed signals
        self.analyse_linewidth_final_list = []

    def get_te_list(self):
        """
        Return the TEs for all the data signals in this pipeline.

        Returns
        -------
        te_list : list
            TEs for all signals stored in here
        """
        te_list = []
        for s in self._data_list:
            te_list.append(s.te)

        return(te_list)

    def run_pipeline_std(self):
        """
        Run the standard pipeline for MRS reco. It includes.

            * reading the data
            * phasing
            * combining channels
            * zero-filling
            * peak analysis & data rejection
            * realigning
            * averaging
            * apodizing
            * removing water reidue with HSVD
            * shifting the spectrum
            * analyzing the SNR
            * analyzing the linewidth
            * displaying the spectrum

        Returns
        -------
        s : MRSData2 numpy array [timepoints]
            Final MRS data signal obtained from reconstruction pipeline stored in a MRSData2 object
        """
        print("")
        print("------------------------------------------------------------------")
        print("mrs.reco.pipeline.run_pipeline_std:")
        print("------------------------------------------------------------------")
        print("> phase_enable = ", end="", flush=True)
        cprint("%r" % self.phase_enable, ('green' if self.phase_enable else 'red'), attrs=['bold'])
        print("> fid_modulus = ", end="", flush=True)
        cprint("%r" % self.fid_modulus, ('green' if self.fid_modulus else 'red'), attrs=['bold'])
        print("> recombine_enable = ", end="", flush=True)
        cprint("%r" % self.recombine_enable, ('green' if self.recombine_enable else 'red'), attrs=['bold'])
        print("> data_concatenate = ", end="", flush=True)
        cprint("%r" % self.data_concatenate, ('green' if self.data_concatenate else 'red'), attrs=['bold'])
        print("> zerofill_enable = ", end="", flush=True)
        cprint("%r" % self.zerofill_enable, ('green' if self.zerofill_enable else 'red'), attrs=['bold'])
        print("> analyse_physio_enable = ", end="", flush=True)
        cprint("%r" % self.analyse_physio_enable, ('green' if self.analyse_physio_enable else 'red'), attrs=['bold'])
        print("> analyse_and_reject_enable = ", end="", flush=True)
        cprint("%r" % self.analyse_and_reject_enable, ('green' if self.analyse_and_reject_enable else 'red'), attrs=['bold'])
        print("> realign_enable = ", end="", flush=True)
        cprint("%r" % self.realign_enable, ('green' if self.realign_enable else 'red'), attrs=['bold'])
        print("> average_enable = ", end="", flush=True)
        cprint("%r" % self.average_enable, ('green' if self.average_enable else 'red'), attrs=['bold'])
        print("> apodize_enable = ", end="", flush=True)
        cprint("%r" % self.apodize_enable, ('green' if self.apodize_enable else 'red'), attrs=['bold'])
        print("> remove_water_enable = ", end="", flush=True)
        cprint("%r" % self.remove_water_enable, ('green' if self.remove_water_enable else 'red'), attrs=['bold'])
        print("> calibrate_enable = ", end="", flush=True)
        cprint("%r" % self.calibrate_enable, ('green' if self.calibrate_enable else 'red'), attrs=['bold'])
        print("> display_enable = ", end="", flush=True)
        cprint("%r" % self.display_enable, ('green' if self.display_enable else 'red'), attrs=['bold'])
        print("> analyse_snr_enable = ", end="", flush=True)
        cprint("%r" % self.analyse_snr_enable, ('green' if self.analyse_snr_enable else 'red'), attrs=['bold'])
        print("> analyse_linewidth_enable = ", end="", flush=True)
        cprint("%r" % self.analyse_linewidth_enable, ('green' if self.analyse_linewidth_enable else 'red'), attrs=['bold'])
        print("> analyse_snr_evol_enable = ", end="", flush=True)
        cprint("%r" % self.analyse_snr_evol_enable, ('green' if self.analyse_snr_evol_enable else 'red'), attrs=['bold'])
        print("> analyse_linewidth_evol_enable = ", end="", flush=True)
        cprint("%r" % self.analyse_linewidth_evol_enable, ('green' if self.analyse_linewidth_evol_enable else 'red'), attrs=['bold'])

        print("")
        print("------------------------------------------------------------------")
        print("> checking some stuff...")
        print("------------------------------------------------------------------")

        # data files: parse
        print("> parsing data file paths...")
        if(self.data_filepaths == []):
            parsed_list = None
        elif(type(self.data_filepaths) == list):
            parsed_list = self.data_filepaths
        else:
            parsed_list = _parse_string_into_list(self.data_filepaths)

        if(parsed_list is not None):
            self._data_filepaths_list = parsed_list
        else:
            raise Exception("> no data scan files specified!")

        # ref data files: parse
        print("> parsing ref. data file paths...")

        if(self.data_ref_filepaths == []):
            parsed_list = None
        elif(type(self.data_ref_filepaths) == list):
            parsed_list = self.data_ref_filepaths
        else:
            parsed_list = _parse_string_into_list(self.data_ref_filepaths)

        if(parsed_list is not None):
            self._data_ref_filepaths_list = parsed_list
        else:
            self._data_ref_filepaths_list = [None] * len(self._data_filepaths_list)

        # physio files: parse
        print("> parsing physio data file paths...")
        if(self.data_physio_filepaths == []):
            parsed_list = None
        elif(type(self.data_physio_filepaths) == list):
            parsed_list = self.data_physio_filepaths
        else:
            parsed_list = _parse_string_into_list(self.data_physio_filepaths)

        if(parsed_list is not None):
            self._data_physio_list = parsed_list
        else:
            self._data_physio_list = [None] * len(self._data_filepaths_list)

        # legend captions: parse
        print("> parsing legend captions...")
        if(self.display_legends == []):
            parsed_list = None
        elif(type(self.display_legends) == list):
            parsed_list = self.display_legends
        else:
            parsed_list = _parse_string_into_list(self.display_legends)

        if(parsed_list is not None):
            self._display_legends_list = parsed_list
        else:
            self._display_legends_list = [""] * len(self._data_filepaths_list)

        # amplitude factors: "parse"
        if(len(self.display_amp_factor_list) == 0):
            self.display_amp_factor_list = np.ones([len(self._data_filepaths_list), ])

        # before reading data and stuff, let's check the list dimensions are consistent
        if(len(self._data_ref_filepaths_list) == 0):
            Warning("> no data ref. scans specified, that's ok, it is optional...")
        if(len(self._data_physio_list) == 0):
            Warning("> no physio log files specified, that's ok, it is optional...")

        n = len(self._data_filepaths_list)
        if(len(self._data_ref_filepaths_list) != n and len(self._data_ref_filepaths_list) > 0):
            raise Exception("> weird, not the same number of data files (" + str(n) + ") and data ref. scans (" + str(len(self._data_ref_filepaths_list)) + ")?!")
        elif(len(self._data_physio_list) != n and len(self._data_physio_list) > 0):
            raise Exception("> weird, not the same number of data files (" + str(n) + ") and physio log files (" + str(len(self._data_physio_list)) + ")?!")
        elif(len(self._display_legends_list) != n):
            raise Exception("> weird, not the same number of data files (" + str(n) + ") and legend captions (" + str(len(self._display_legends_list)) + ")?!")
        elif(len(self.display_amp_factor_list) != n):
            raise Exception("> weird, not the same number of data files (" + str(n) + ") and display scale factors (" + str(len(self.display_amp_factor_list)) + ")?!")

        # check if there are some blanks
        if(any(d is None for d in self._data_filepaths_list)):
            raise Exception("> hey, you missed a line in the data scan file path list?!")
        if(any(d is None for d in self._display_legends_list)):
            raise Exception("> hey, you missed a line in the legend caption list?!")

        # oh, we want to process only one dataset in the list
        if(len(self.data_process_only_this_data_index) > 0):
            if(len(self._data_filepaths_list) > 0):
                self._data_filepaths_list = [self._data_filepaths_list[i] for i in self.data_process_only_this_data_index]
            if(len(self._data_ref_filepaths_list) > 0):
                self._data_ref_filepaths_list = [self._data_ref_filepaths_list[i] for i in self.data_process_only_this_data_index]
            if(len(self._data_physio_list) > 0):
                self._data_physio_list = [self._data_physio_list[i] for i in self.data_process_only_this_data_index]
            if(len(self._display_legends_list) > 0):
                self._display_legends_list = [self._display_legends_list[i] for i in self.data_process_only_this_data_index]
            if(len(self.display_amp_factor_list) > 0):
                self.display_amp_factor_list = [self.display_amp_factor_list[i] for i in self.data_process_only_this_data_index]

        # now let's read the data files
        print("")
        print("------------------------------------------------------------------")
        print("> reading data files...")
        max_len_patient_name = 0
        for data_fn, physio_fn, data_leg in zip(self._data_filepaths_list, self._data_physio_list, self._display_legends_list):
            print("------------------------------------------------------------------")
            cprint("> reading data [" + data_leg + "]", 'green', attrs=['bold'])
            s = MRSData2(data_fn, self.data_coil_nChannels, physio_fn)
            print("> got a " + str(s.shape) + " vector")
            # store
            self._data_list.append(s)
            # patient name length
            if(s.patient_name is not None and len(s.patient_name) > max_len_patient_name):
                max_len_patient_name = len(s.patient_name)

        print("------------------------------------------------------------------")
        print("> reading ref. data files...")
        for data_ref_fn, data_leg in zip(self._data_ref_filepaths_list, self._display_legends_list):
            if(data_ref_fn is not None):
                print("------------------------------------------------------------------")
                cprint("> reading ref. data [" + data_leg + "]", 'green', attrs=['bold'])
                s = MRSData2(data_ref_fn, self.data_coil_nChannels, None)
                print("> got a " + str(s.shape) + " vector")
                # if several averages, mean now
                s = s.mean(axis=0)
                # add 1 dimension
                s = s.reshape((1,) + s.shape)
                print("> reshaped to a " + str(s.shape) + " vector")
            else:
                s = None
            # store
            self._data_ref_list.append(s)

        # legends: add patient name and format (space pads)
        nchar = 5 + max_len_patient_name + len(max(self._display_legends_list, key=len))
        display_legends_list_formatted = []
        for i, l, d in zip(range(len(self._display_legends_list)), self._display_legends_list, self._data_list):
            if(d.patient_name is not None):
                this_new_legend = "[" + str(i) + "] " + d.patient_name + " - " + l
            else:
                this_new_legend = l

            this_new_legend = this_new_legend.ljust(nchar)
            display_legends_list_formatted.append(this_new_legend)

        self._display_legends_list = display_legends_list_formatted

        # set ppm0 for all signals
        for ind, d in enumerate(self._data_list):
            self._data_list[ind].ppm0 = self.ppm0
        for ind, d in enumerate(self._data_ref_list):
            if d is not None:
                self._data_ref_list[ind].ppm0 = self.ppm0

        print("")

        # rephase
        if(self.phase_enable):
            print("------------------------------------------------------------------")
            print("> rephasing...")
            print("------------------------------------------------------------------")
            for i in range(0, len(self._data_list)):
                cprint("> rephasing [" + self._display_legends_list[i] + "]", 'green', attrs=['bold'])

                s = self._data_list[i]
                s_ref = self._data_ref_list[i]
                s_phased = s.correct_phase(s_ref, self.phase_POI_range_ppm, self.phase_weak_ws_mode, self.phase_high_snr_mode, self.phase_order, self.phase_offset, self.phase_display, self.display_range_ppm)
                # replace / store
                self._data_list[i] = s_phased
                print("------------------------------------------------------------------")

                if(s_ref is not None):
                    cprint("> rephasing ref. data for [" + self._display_legends_list[i] + "]", 'green', attrs=['bold'])

                    s_ref_phased = s_ref.correct_phase(s_ref, self.phase_POI_range_ppm, self.phase_weak_ws_mode, self.phase_high_snr_mode, self.phase_order, self.phase_offset, False, self.display_range_ppm)
                    # replace / store
                    self._data_ref_list[i] = s_ref_phased
                    print("------------------------------------------------------------------")
            print("")

        # FID modulus
        if(self.fid_modulus):
            print("------------------------------------------------------------------")
            print("> FID modulus...")
            print("------------------------------------------------------------------")
            for i in range(0, len(self._data_list)):
                cprint("> FID modulus of [" + self._display_legends_list[i] + "]", 'green', attrs=['bold'])
                s = self._data_list[i]
                s_fidmod = s.correct_fidmodulus()
                # replace / store
                self._data_list[i] = s_fidmod
            print("------------------------------------------------------------------")
            print("")

        # recombine
        if(self.recombine_enable):
            print("------------------------------------------------------------------")
            print("> recombining channels...")
            print("------------------------------------------------------------------")
            # recombine using reference data
            for i in range(0, len(self._data_list)):
                cprint("> recombining [" + self._display_legends_list[i] + "]", 'green', attrs=['bold'])
                # try and use a ref scan for rephasing if available
                # checking if use_data_ref flag is consistent with the data
                if(self.recombine_use_data_ref):
                    if(len(self._data_ref_list) > 0 and self._data_ref_list[i] is not None):
                        s_ref = self._data_ref_list[i]
                    else:
                        Warning("> no ref. data scan available for this data!")
                        s_ref = None
                else:
                    s_ref = None
                s = self._data_list[i]
                s_combined = s.correct_combine_channels(s_ref, self.recombine_phasing, self.recombine_weights)

                # replace / store
                self._data_list[i] = s_combined
                print("------------------------------------------------------------------")

                if(s_ref is not None):
                    cprint("> recombining ref. data for [" + self._display_legends_list[i] + "]", 'green', attrs=['bold'])

                    s_ref_combined = s_ref.correct_combine_channels(s_ref, self.recombine_phasing, self.recombine_weights)

                    # replace / store
                    self._data_ref_list[i] = s_ref_combined
                    print(
                        "------------------------------------------------------------------")
            print("")

        # concatenate data or process / display separatly?
        if(len(self._data_list) > 1 and self.data_concatenate):
            print("------------------------------------------------------------------")
            print("> concatenating data...")
            print("------------------------------------------------------------------")
            s_concatenated = self._data_list[0]
            for i in range(1, len(self._data_list)):
                s = self._data_list[i]
                s_concatenated = np.concatenate((s_concatenated, s))
                s_concatenated = s.inherit(s_concatenated)

            # empty and store the single concatenated signal in the data list
            self._data_list = []
            self._data_list.append(s_concatenated)
            print("------------------------------------------------------------------")
            print("")

        # and go on with the processing
        d_offset = 0
        for (k, s, s_ref, s_legend, d_factor) in zip(range(len(self._data_list)), self._data_list, self._data_ref_list, self._display_legends_list, self.display_amp_factor_list):

            # check initial snr
            if(self.analyse_snr_enable and self.analyse_snr_evol_enable):
                print("------------------------------------------------------------------")
                this_snr = s._correct_realign()._correct_average()._correct_apodization()._correct_freqshift().analyse_snr(self.analyse_snr_s_range_ppm, self.analyse_snr_n_range_ppm, s_legend, self.analyse_snr_area_integrate, self.analyse_snr_magnitude_mode, self.analyse_snr_display, self.display_range_ppm, True)
                print("> initial SNR of [" + s_legend + "] = %.2f" % this_snr)
                self.analyse_snr_evol_list.append([])
                self.analyse_snr_evol_list[k].append(this_snr)
                self.analyse_snr_evol_list_labels.append([])
                self.analyse_snr_evol_list_labels[k].append('initial')
                print("------------------------------------------------------------------")
                print("")

            # check initial lw
            if(self.analyse_linewidth_enable and self.analyse_linewidth_evol_enable):
                print("------------------------------------------------------------------")
                this_lw = s._correct_realign()._correct_average()._correct_apodization()._correct_freqshift().analyse_linewidth(self.analyse_linewidth_range_ppm, s_legend, self.analyse_linewidth_magnitude_mode, self.analyse_linewidth_display, self.display_range_ppm, True)
                print("> initial LW of [" + s_legend + "] = %.2f" % this_lw)
                self.analyse_linewidth_evol_list.append([])
                self.analyse_linewidth_evol_list[k].append(this_lw)
                self.analyse_linewidth_evol_list_labels.append([])
                self.analyse_linewidth_evol_list_labels[k].append('initial')
                print("------------------------------------------------------------------")
                print("")

            # zero-filling
            if(self.zerofill_enable):
                print("------------------------------------------------------------------")
                print("> zero-filling data...")
                print("------------------------------------------------------------------")
                cprint("> zero-filling [" + s_legend + "]", 'green', attrs=['bold'])
                s = s.correct_zerofill(self.zerofill_npts, self.zerofill_display, self.display_range_ppm)

                # check snr
                if(self.analyse_snr_enable and self.analyse_snr_evol_enable):
                    print("------------------------------------------------------------------")
                    this_snr = s._correct_realign()._correct_average()._correct_apodization()._correct_freqshift().analyse_snr(self.analyse_snr_s_range_ppm, self.analyse_snr_n_range_ppm, s_legend, self.analyse_snr_area_integrate, self.analyse_snr_magnitude_mode, self.analyse_snr_display, self.display_range_ppm, True)
                    print("> post zero-filling SNR of [" + s_legend + "] = %.2f" % this_snr)
                    self.analyse_snr_evol_list[k].append(this_snr)
                    self.analyse_snr_evol_list_labels[k].append('post zero-filling')

                # check lw
                if(self.analyse_linewidth_enable and self.analyse_linewidth_evol_enable):
                    print("------------------------------------------------------------------")
                    this_lw = s._correct_realign()._correct_average()._correct_apodization()._correct_freqshift().analyse_linewidth(self.analyse_linewidth_range_ppm, s_legend, self.analyse_linewidth_magnitude_mode, self.analyse_linewidth_display, self.display_range_ppm, True)
                    print("> post zero-filling LW of [" + s_legend + "] = %.2f" % this_lw)
                    self.analyse_linewidth_evol_list[k].append(this_lw)
                    self.analyse_linewidth_evol_list_labels[k].append('post zero-filling')

                print("------------------------------------------------------------------")
                print("")

            # analysis of physiological stuff
            if(self.analyse_physio_enable):
                print("------------------------------------------------------------------")
                print("> data and physio analysis...")
                print("------------------------------------------------------------------")
                cprint("> data and physio analysis [" + s_legend + "]", 'green', attrs=['bold'])
                s.analyse_physio(self.analyse_physio_POI_range_ppm, self.analyse_physio_delta_time_ms, self.analyse_physio_display)
                print("------------------------------------------------------------------")
                print("")

            # analysis and data rejection
            if(self.analyse_and_reject_enable):
                print("------------------------------------------------------------------")
                print("> data analysis and rejection...")
                print("------------------------------------------------------------------")
                cprint("> data analysing / rejecting [" + s_legend + "]", 'green', attrs=['bold'])
                s = s.correct_analyse_and_reject(self.analyse_and_reject_POI_range_ppm, self.analyse_and_reject_moving_averages, self.analyse_and_reject_min, self.analyse_and_reject_max, self.analyse_and_reject_relative_mean, self.analyse_and_reject_auto, self.analyse_and_reject_display)

                # check snr
                if(self.analyse_snr_enable and self.analyse_snr_evol_enable):
                    print("------------------------------------------------------------------")
                    this_snr = s._correct_realign()._correct_average()._correct_apodization()._correct_freqshift().analyse_snr(self.analyse_snr_s_range_ppm, self.analyse_snr_n_range_ppm, s_legend, self.analyse_snr_area_integrate, self.analyse_snr_magnitude_mode, self.analyse_snr_display, self.display_range_ppm, True)
                    print("> post analyse / reject SNR of [" + s_legend + "] = %.2f" % this_snr)
                    self.analyse_snr_evol_list[k].append(this_snr)
                    self.analyse_snr_evol_list_labels[k].append('post analyse / reject')

                # check lw
                if(self.analyse_linewidth_enable and self.analyse_linewidth_evol_enable):
                    print("------------------------------------------------------------------")
                    this_lw = s._correct_realign()._correct_average()._correct_apodization()._correct_freqshift().analyse_linewidth(self.analyse_linewidth_range_ppm, s_legend, self.analyse_linewidth_magnitude_mode, self.analyse_linewidth_display, self.display_range_ppm, True)
                    print("> post analyse / reject LW of [" + s_legend + "] = %.2f" % this_lw)
                    self.analyse_linewidth_evol_list[k].append(this_lw)
                    self.analyse_linewidth_evol_list_labels[k].append('post analyse / reject')

                print("------------------------------------------------------------------")
                print("")

            # realignement
            if(self.realign_enable):
                print("------------------------------------------------------------------")
                print("> frequency realignment...")
                print("------------------------------------------------------------------")
                cprint("> realigning [" + s_legend + "]", 'green', attrs=['bold'])
                s = s.correct_realign(self.realign_POI_range_ppm, self.realign_moving_averages, self.realign_display, self.display_range_ppm)

                # check snr
                if(self.analyse_snr_enable and self.analyse_snr_evol_enable):
                    print("------------------------------------------------------------------")
                    this_snr = s._correct_average()._correct_apodization()._correct_freqshift().analyse_snr(self.analyse_snr_s_range_ppm, self.analyse_snr_n_range_ppm, s_legend, self.analyse_snr_area_integrate, self.analyse_snr_magnitude_mode, self.analyse_snr_display, self.display_range_ppm, True)
                    print("> post realignement SNR of [" + s_legend + "] = %.2f" % this_snr)
                    self.analyse_snr_evol_list[k].append(this_snr)
                    self.analyse_snr_evol_list_labels[k].append('post realignement')

                # check lw
                if(self.analyse_linewidth_enable and self.analyse_linewidth_evol_enable):
                    print("------------------------------------------------------------------")
                    this_lw = s._correct_average()._correct_apodization()._correct_freqshift().analyse_linewidth(self.analyse_linewidth_range_ppm, s_legend, self.analyse_linewidth_magnitude_mode, self.analyse_linewidth_display, self.display_range_ppm, True)
                    print("> post realignement LW of [" + s_legend + "] = %.2f" % this_lw)
                    self.analyse_linewidth_evol_list[k].append(this_lw)
                    self.analyse_linewidth_evol_list_labels[k].append('post realignement')

                print("------------------------------------------------------------------")
                print("")

            # average
            if(self.average_enable):
                print("------------------------------------------------------------------")
                print("> averaging...")
                print("------------------------------------------------------------------")
                cprint("> averaging [" + s_legend + "]", 'green', attrs=['bold'])
                s = s.correct_average(self.average_na, self.average_display, self.display_range_ppm)

                # check snr
                if(self.analyse_snr_enable and self.analyse_snr_evol_enable):
                    print("------------------------------------------------------------------")
                    this_snr = s._correct_apodization()._correct_freqshift().analyse_snr(self.analyse_snr_s_range_ppm, self.analyse_snr_n_range_ppm, s_legend, self.analyse_snr_area_integrate, self.analyse_snr_magnitude_mode, self.analyse_snr_display, self.display_range_ppm, True)
                    print("> post averaging SNR of [" + s_legend + "] = %.2f" % this_snr)
                    self.analyse_snr_evol_list[k].append(this_snr)
                    self.analyse_snr_evol_list_labels[k].append('post averaging')

                # check lw
                if(self.analyse_linewidth_enable and self.analyse_linewidth_evol_enable):
                    print("------------------------------------------------------------------")
                    this_lw = s._correct_apodization()._correct_freqshift().analyse_linewidth(self.analyse_linewidth_range_ppm, s_legend, self.analyse_linewidth_magnitude_mode, self.analyse_linewidth_display, self.display_range_ppm, True)
                    print("> post averaging LW of [" + s_legend + "] = %.2f" % this_lw)
                    self.analyse_linewidth_evol_list[k].append(this_lw)
                    self.analyse_linewidth_evol_list_labels[k].append('post averaging')

                print("------------------------------------------------------------------")
                print("")

            # apodize
            if(self.apodize_enable):
                print("------------------------------------------------------------------")
                print("> apodizing...")
                print("------------------------------------------------------------------")
                cprint("> apodizing [" + s_legend + "]", 'green', attrs=['bold'])

                s = s.correct_apodization(self.apodize_damping_hz, self.apodize_npts, self.apodize_display, self.display_range_ppm)

                if(s_ref is not None):
                    print("------------------------------------------------------------------")
                    cprint("> apodizing (cropping really) ref. data for [" + s_legend + "]", 'green', attrs=['bold'])
                    s_ref = s_ref._correct_average().correct_apodization(1.0, self.apodize_npts, False, self.display_range_ppm)

                # check snr
                if(self.analyse_snr_enable and self.analyse_snr_evol_enable):
                    print("------------------------------------------------------------------")
                    this_snr = s._correct_freqshift().analyse_snr(self.analyse_snr_s_range_ppm, self.analyse_snr_n_range_ppm, s_legend, self.analyse_snr_area_integrate, self.analyse_snr_magnitude_mode, self.analyse_snr_display, self.display_range_ppm, True)
                    print("> post apodization SNR of [" + s_legend + "] = %.2f" % this_snr)
                    self.analyse_snr_evol_list[k].append(this_snr)
                    self.analyse_snr_evol_list_labels[k].append('post apodization')

                # check lw
                if(self.analyse_linewidth_enable and self.analyse_linewidth_evol_enable):
                    print("------------------------------------------------------------------")
                    this_lw = s._correct_freqshift().analyse_linewidth(self.analyse_linewidth_range_ppm, s_legend, self.analyse_linewidth_magnitude_mode, self.analyse_linewidth_display, self.display_range_ppm, True)
                    print("> post apodization LW of [" + s_legend + "] = %.2f" % this_lw)
                    self.analyse_linewidth_evol_list[k].append(this_lw)
                    self.analyse_linewidth_evol_list_labels[k].append('post apodization')

                print("------------------------------------------------------------------")
                print("")

            # water residue removal
            if(self.remove_water_enable):
                print("------------------------------------------------------------------")
                print("> removing water...")
                print("------------------------------------------------------------------")
                cprint("> removing water from [" + s_legend + "]", 'green', attrs=['bold'])
                s = s.correct_water_removal(self.remove_water_hsvd_components, self.remove_water_hsvd_range, self.remove_water_display, self.display_range_ppm)
                print("------------------------------------------------------------------")
                print("")

            # frequency shift calibration
            if(self.calibrate_enable):
                print("------------------------------------------------------------------")
                print("> calibrating spectrum...")
                print("------------------------------------------------------------------")
                cprint("> calibrating [" + s_legend + "]", 'green', attrs=['bold'])
                s = s.correct_freqshift(self.calibrate_POI_range_ppm, self.calibrate_POI_true_ppm, self.calibrate_display, self.display_range_ppm)
                print("------------------------------------------------------------------")
                print("")

            # display
            if(self.display_enable):
                print("------------------------------------------------------------------")
                print("> display spectrum...")
                print("------------------------------------------------------------------")
                cprint("> displaying [" + s_legend + "]", 'green', attrs=['bold'])
                s.display_spectrum(self.display_fig_index, s_legend, self.display_range_ppm, d_factor, d_offset, self.display_magnitude_mode)
                d_offset = d_offset + self.display_offset
                print("------------------------------------------------------------------")
                print("")

            # snr estimation
            if(self.analyse_snr_enable):
                print("------------------------------------------------------------------")
                print("> estimating final SNR...")
                print("------------------------------------------------------------------")
                cprint("> estimating final SNR of [" + s_legend + "]", 'green', attrs=['bold'])
                this_snr = s.analyse_snr(self.analyse_snr_s_range_ppm, self.analyse_snr_n_range_ppm, s_legend, self.analyse_snr_area_integrate, self.analyse_snr_magnitude_mode, self.analyse_snr_display, self.display_range_ppm)
                self.analyse_snr_final_list.append(this_snr)
                if(self.analyse_snr_evol_enable):
                    self.analyse_snr_evol_list[k].append(this_snr)
                    self.analyse_snr_evol_list_labels[k].append('final')
                print("------------------------------------------------------------------")
                print("")

            # lw estimation
            if(self.analyse_linewidth_enable):
                print("------------------------------------------------------------------")
                print("> estimating final linewidth...")
                print("------------------------------------------------------------------")
                cprint("> estimating final linewidth of [" + s_legend + "]", 'green', attrs=['bold'])
                this_lw = s.analyse_linewidth(self.analyse_linewidth_range_ppm, s_legend, self.analyse_linewidth_magnitude_mode, self.analyse_linewidth_display, self.display_range_ppm)
                self.analyse_linewidth_final_list.append(this_lw)
                if(self.analyse_linewidth_evol_enable):
                    self.analyse_linewidth_evol_list[k].append(this_lw)
                    self.analyse_linewidth_evol_list_labels[k].append('final')
                print("------------------------------------------------------------------")
                print("")

        # summary display in console
        if(self.analyse_snr_enable):
            print("------------------------------------------------------------------")
            print("---SNR results summary--------------------------------------------")
            print("------------------------------------------------------------------")
            for (this_legend, this_snr) in zip(self._display_legends_list, self.analyse_snr_final_list):
                print("> " + this_legend + "\tSNR=%.2f" % this_snr)

            # plot snr evolution
            if(self.analyse_snr_evol_enable):
                plt.figure(210).canvas.set_window_title("mrs.reco.pipeline.run_pipeline_std.analyse_snr_evol_enable")
                plt.clf()
                plt.bar(np.arange(len(self.analyse_snr_evol_list[k])), self.analyse_snr_evol_list[k])
                plt.plot(np.arange(len(self.analyse_snr_evol_list[k])), self.analyse_snr_evol_list[k], 'kx-')
                plt.xticks(np.arange(len(self.analyse_snr_evol_list[k])), labels=self.analyse_snr_evol_list_labels[k], rotation=45)
                plt.ylabel('Estimated SNR (u.a)')
                plt.grid('on')
                plt.tight_layout()
            print("------------------------------------------------------------------")
            print("")

        if(self.analyse_linewidth_enable):
            print("------------------------------------------------------------------")
            print("---LW results summary---------------------------------------------")
            print("------------------------------------------------------------------")
            for (this_legend, this_lw) in zip(self._display_legends_list, self.analyse_linewidth_final_list):
                print("> " + this_legend + "\tLW=%.2f Hz" % this_lw)

            # plot snr evolution
            if(self.analyse_linewidth_evol_enable):
                plt.figure(211).canvas.set_window_title("mrs.reco.pipeline.run_pipeline_std.analyse_linewidth_evol_enable")
                plt.clf()
                plt.bar(np.arange(len(self.analyse_linewidth_evol_list[k])), self.analyse_linewidth_evol_list[k])
                plt.plot(np.arange(len(self.analyse_linewidth_evol_list[k])), self.analyse_linewidth_evol_list[k], 'kx-')
                plt.xticks(np.arange(len(self.analyse_linewidth_evol_list[k])), labels=self.analyse_linewidth_evol_list_labels[k], rotation=45)
                plt.ylabel('Estimated LW (Hz)')
                plt.grid('on')
                plt.tight_layout()
            print("------------------------------------------------------------------")
            print("")
            print("------------------------------------------------------------------")
            print("")

        print("> Pipeline terminated!")
        print("> Returning last processed signal and its corresponding ref. signal!")
        return(s, s_ref)


class voi_pipeline:
    """The voi_pipeline class is similar to pipeline class above but for VOI trace signals."""

    def __init__(self):
        self.data_filepaths = ""
        self._data_filepaths_list = []
        self._data_list = []
        self._xdata = []
        self._ydata = []
        self._zdata = []
        self._xdata_axis = []
        self._ydata_axis = []
        self._zdata_axis = []

        self.display_fig_index = 1
        self.display_legends = ""
        self._display_legends_list = []

        self.analyze_selectivity_range_list = [[800, 3550], [-10600, -7800], [-3650, 1850]]
        self.analyze_selectivity_list = np.array([])

    def get_te_list(self):
        """
        Return the TEs for all the data signals in this pipeline.

        Returns
        -------
        te_list : list
            TEs for all signals stored in here
        """
        te_list = []
        for s in self._data_list:
            te_list.append(s[0].te)

        return(te_list)

    def run_pipeline_std(self):
        """
        Run the standard pipeline for VOI trace signals. It includes.

            * reading the data
            * displaying the spatial selection profiles
            * analyzing the selectivity
        """
        print("")
        print("------------------------------------------------------------------")
        print("mrs.reco.voi_pipeline.run_pipeline_std:")
        print("------------------------------------------------------------------")
        print("")

        # parse
        print("------------------------------------------------------------------")
        print("> parsing data file paths...")
        print("------------------------------------------------------------------")
        if(self.data_filepaths == []):
            parsed_list = None
        elif(type(self.data_filepaths) == list):
            parsed_list = self.data_filepaths
        else:
            parsed_list = _parse_string_into_list(self.data_filepaths)

        if(parsed_list is None):
            raise Exception("> no data scan files specified!")
        else:
            self._data_filepaths_list = parsed_list
        print("------------------------------------------------------------------")
        print("")

        # legends
        print("------------------------------------------------------------------")
        print("> parsing legend captions...")
        print("------------------------------------------------------------------")
        if(self.display_legends == []):
            parsed_list = None
        elif(type(self.display_legends) == list):
            parsed_list = self.display_legends
        else:
            parsed_list = _parse_string_into_list(self.display_legends)

        if(parsed_list is None):
            # oh no, no legend captions were specified
            # let's create them using filenames
            self._display_legends_list = _build_legend_list_from_filepath_list(
                self._data_filepaths_list)
        else:
            self._display_legends_list = parsed_list
        print("------------------------------------------------------------------")
        print("")

        # read the 3 VOI traces per dataset
        print("------------------------------------------------------------------")
        print("> reading data files...")
        print("------------------------------------------------------------------")
        for f in self._data_filepaths_list:
            # read data
            print("> looking in folder: ")
            print(f)
            print("> reading the 3 dicom files...")
            sx = suspect.io.load_siemens_dicom(f + "/original-primary_e09_0001.dcm")
            sy = suspect.io.load_siemens_dicom(f + "/original-primary_e09_0002.dcm")
            sz = suspect.io.load_siemens_dicom(f + "/original-primary_e09_0003.dcm")
            self._data_list.append([sx, sy, sz])
            # normalize
            sx_spectrum = np.abs(sx.spectrum()) / np.abs(sx.spectrum()).max()
            sy_spectrum = np.abs(sy.spectrum()) / np.abs(sy.spectrum()).max()
            sz_spectrum = np.abs(sz.spectrum()) / np.abs(sz.spectrum()).max()
            # store
            self._xdata.append(sx_spectrum)
            self._ydata.append(sy_spectrum)
            self._zdata.append(sz_spectrum)
            self._xdata_axis.append(sx.frequency_axis())
            self._ydata_axis.append(sy.frequency_axis())
            self._zdata_axis.append(sz.frequency_axis())
        print("------------------------------------------------------------------")
        print("")

        # analysis
        print("------------------------------------------------------------------")
        print("> evaluating selectivity using a " + str(self.analyze_selectivity_range_list) + "ppm ranges...")
        print("------------------------------------------------------------------")
        self.analyze_selectivity_list = np.zeros([len(self._data_filepaths_list), 3, 2])
        k = 0
        for (sx, sy, sz, sx_ax, sy_ax, sz_ax, leg) in zip(self._xdata, self._ydata, self._zdata, self._xdata_axis, self._ydata_axis, self._zdata_axis, self._display_legends_list):

            # find ppm indexes corresponding to ppm values
            sx_ithreshold_l = np.argmin(np.abs(sx_ax - self.analyze_selectivity_range_list[0][0]))
            sx_ithreshold_r = np.argmin(np.abs(sx_ax - self.analyze_selectivity_range_list[0][1]))
            sy_ithreshold_l = np.argmin(np.abs(sy_ax - self.analyze_selectivity_range_list[1][0]))
            sy_ithreshold_r = np.argmin(np.abs(sy_ax - self.analyze_selectivity_range_list[1][1]))
            sz_ithreshold_l = np.argmin(np.abs(sz_ax - self.analyze_selectivity_range_list[2][0]))
            sz_ithreshold_r = np.argmin(np.abs(sz_ax - self.analyze_selectivity_range_list[2][1]))

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

            print("> selectivity results for [" + leg + "]...")
            print("> [X] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sx_in, sx_out))
            print("> [Y] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sy_in, sy_out))
            print("> [Z] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sz_in, sz_out))

            # store
            self.analyze_selectivity_list[k, :, 0] = [sx_in, sy_in, sz_in]
            self.analyze_selectivity_list[k, :, 1] = [sx_out, sy_out, sz_out]
            k = k + 1

            '''
            print("> mrs.reco.voi_pipeline.run_pipeline_std: Selectivity debug info for [" + leg + "] coming...")
            print("> [X] Intersection at %f Hz and %f Hz" % (sx_ax[sx_ithreshold_l],sx_ax[sx_ithreshold_r]))
            print("> [Y] Intersection at %f Hz and %f Hz" % (sy_ax[sy_ithreshold_l],sy_ax[sy_ithreshold_r]))
            print("> [Z] Intersection at %f Hz and %f Hz" % (sz_ax[sz_ithreshold_l],sz_ax[sz_ithreshold_r]))
            '''
        print("------------------------------------------------------------------")
        print("")

        # display
        sx_min_all = 0
        sx_max_all = 0
        sy_min_all = 0
        sy_max_all = 0
        sz_min_all = 0
        sz_max_all = 0
        print("------------------------------------------------------------------")
        print("> displaying the 3 spatial profiles...")
        print("------------------------------------------------------------------")
        for (sx, sy, sz, sx_ax, sy_ax, sz_ax, leg) in zip(self._xdata, self._ydata, self._zdata, self._xdata_axis, self._ydata_axis, self._zdata_axis, self._display_legends_list):

            # record min / max
            if(sx.min() < sx_min_all):
                sx_min_all = sx.min()
            if(sx.max() > sx_max_all):
                sx_max_all = sx.max()
            if(sy.min() < sy_min_all):
                sy_min_all = sy.min()
            if(sy.max() > sy_max_all):
                sy_max_all = sy.max()
            if(sz.min() < sz_min_all):
                sz_min_all = sz.min()
            if(sz.max() > sz_max_all):
                sz_max_all = sz.max()

            # display
            plt.figure(self.display_fig_index).canvas.set_window_title("mrs.reco.read_trace_data_and_display")

            plt.subplot(2, 4, 1)
            plt.plot(sx_ax, sx, linewidth=1, label=leg)
            plt.xlabel('frequency dispersion (Hz)')
            plt.grid('on')

            plt.subplot(2, 4, 5)
            plt.plot(sx_ax, sx, linewidth=1, label=leg)
            plt.xlabel('frequency dispersion (Hz)')
            plt.ylabel('X spatial profile')
            plt.grid('on')

            plt.subplot(2, 4, 2)
            plt.plot(sy_ax, sy, linewidth=1, label=leg)
            plt.xlabel('frequency dispersion (Hz)')
            plt.grid('on')

            plt.subplot(2, 4, 6)
            plt.plot(sy_ax, sy, linewidth=1, label=leg)
            plt.xlabel('frequency dispersion (Hz)')
            plt.ylabel('Y spatial profile')
            plt.grid('on')

            plt.subplot(2, 4, 3)
            plt.plot(sz_ax, sz, linewidth=1, label=leg)
            plt.xlabel('frequency dispersion (Hz)')
            plt.grid('on')

            plt.subplot(2, 4, 7)
            plt.plot(sz_ax, sz, linewidth=1, label=leg)
            plt.xlabel('frequency dispersion (Hz)')
            plt.ylabel('Z spatial profile')
            plt.grid('on')

        # display ppm range lines
        for (sx, sy, sz, sx_ax, sy_ax, sz_ax) in zip(self._xdata, self._ydata, self._zdata, self._xdata_axis, self._ydata_axis, self._zdata_axis):

            plt.subplot(2, 4, 1)
            sx_minmax = np.array([sx_min_all, sx_max_all])
            sx_ax_l = np.ones(2,) * self.analyze_selectivity_range_list[0][0]
            sx_ax_r = np.ones(2,) * self.analyze_selectivity_range_list[0][1]
            plt.plot(sx_ax_l, sx_minmax, '--', linewidth=0.5)
            plt.plot(sx_ax_r, sx_minmax, '--', linewidth=0.5)
            plt.subplot(2, 4, 5)
            plt.plot(sx_ax_l, sx_minmax, '--', linewidth=0.5)
            plt.plot(sx_ax_r, sx_minmax, '--', linewidth=0.5)

            plt.subplot(2, 4, 2)
            sy_minmax = np.array([sy_min_all, sy_max_all])
            sy_ax_l = np.ones(2,) * self.analyze_selectivity_range_list[1][0]
            sy_ax_r = np.ones(2,) * self.analyze_selectivity_range_list[1][1]
            plt.plot(sy_ax_l, sy_minmax, '--', linewidth=0.5)
            plt.plot(sy_ax_r, sy_minmax, '--', linewidth=0.5)
            plt.subplot(2, 4, 6)
            plt.plot(sy_ax_l, sy_minmax, '--', linewidth=0.5)
            plt.plot(sy_ax_r, sy_minmax, '--', linewidth=0.5)

            plt.subplot(2, 4, 3)
            sz_minmax = np.array([sz_min_all, sz_max_all])
            sz_ax_l = np.ones(2,) * self.analyze_selectivity_range_list[2][0]
            sz_ax_r = np.ones(2,) * self.analyze_selectivity_range_list[2][1]
            plt.plot(sz_ax_l, sz_minmax, '--', linewidth=0.5)
            plt.plot(sz_ax_r, sz_minmax, '--', linewidth=0.5)
            plt.subplot(2, 4, 7)
            plt.plot(sz_ax_l, sz_minmax, '--', linewidth=0.5)
            plt.plot(sz_ax_r, sz_minmax, '--', linewidth=0.5)

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, mode='expand', borderaxespad=0)
        plt.tight_layout()
        plt.show()
        print("")


def _parse_string_into_list(big_string):
    """
    Parse a big string that contains lines of strings in to a list of stripped strings. Usefull for list of filepaths. Any empty lines will give a None element in the list.

    Parameters
    ----------
    big_string : string
        A list of strings separated by line breaks.

    Returns
    -------
    string_list : list
        List of strings
    """
    # first check
    if(len(big_string) == 0):
        return(None)

    # parse
    string_list = []
    big_string_splitted = big_string.splitlines()

    # remove first line
    big_string_splitted = big_string_splitted[1:]
    for p in big_string_splitted:
        this_string = p.strip()
        if(len(this_string) > 0):
            string_list.append(this_string)
        else:
            string_list.append(None)

    return(string_list)


def _build_legend_list_from_filepath_list(filepath_list):
    """
    Try to build a list of legend captions based on the file paths. Any empty filepath will give a None element in the list.

    Parameters
    ----------
    filepath_list : list
        List of file paths
    Returns
    -------
    legend_list : list
        List of captions
    """
    # first check
    if(len(filepath_list) == 0):
        return(None)

    legend_list = []
    for f in filepath_list:
        if(f is None):
            legend_list.append("[No Data]")
        else:
            f_filename, f_extension = os.path.splitext(f.lower())
            if(f_extension == ".dcm"):
                # that is a dicom, extract folder name
                f_a, f_b = os.path.split(f.lower())
                f_a, f_b = os.path.split(f_a)
                legend_list.append(f_b)
            elif(f_extension == ".dat"):
                # that is a twix, extract file name
                f_a, f_b = os.path.split(f.lower())
                legend_list.append(f_b)
            else:
                # no idea what we are dealing with
                f_a, f_b = os.path.split(f.lower())
                legend_list.append(f_b)


def load_broken_twix_vb(filename):
    """
    Read broken TWIX files. A modified version of the load_twix / load_twix_vb function from the suspect library. I added some acq_end exception handling to deal with broken twix files that are obtained after interrupting an acquisition.

    Parameters
    ----------
    filename : string
        Full path to the broken TWIX file

    Returns
    -------
    builder.build_mrsdata() : MRSData object
        Resulting constructed MRSData object
    """
    with open(filename, 'rb') as fin:

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

            meas_uid, scan_counter, time_stamp, pmu_time_stamp = struct.unpack(
                "IIII", fin.read(16))

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
