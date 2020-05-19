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
import suspect.io.twix as sit
import numpy as np
import pydicom
from scipy import signal
import scipy.io as sio
import matplotlib.pylab as plt
import os
from shutil import copyfile
from datetime import datetime, timedelta
import struct
import pickle
import mrs.sim as sim
import mrs.log as log
import mrs.paths as default_paths

import pdb


class data_db():
    """A class used to deal with storage of reconstructed signals with their respective reco pipeline in pkl files. The dumped data consist of a dict with all our processed data, sorted by patient name, dataset name. Also deals with conflicts like already existing patient name/dataset name."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("You are trying to dynamically create a new attribute (%s) to this object and that is not cool! I will not let you do that because I believe it is a bad habit and can lead to terrible bugs. A clean way of doing this is to initialize your attribute (%s) in the __init__ method of this class. Bisou, bye :)" % (key, key))
        object.__setattr__(self, key, value)

    def __init__(self, reco_data_db_file=default_paths.DEFAULT_RECO_DATA_DB_FILE):
        """
        Initialize the reconstructed data storage, basically creates an empty PKL file if nothing already exists.

        Parameters
        ----------
        reco_data_db_file: string
            PKL file where all the data is stored
        """
        # if pkl does not exist, it is our very first time :heart:
        # let's write an empty dict
        if(not os.path.isfile(reco_data_db_file)):
            log.info("creating storage file [%s]..." % reco_data_db_file)
            # write pkl file
            with open(reco_data_db_file, 'wb') as f:
                pickle.dump([{}], f)
        else:
            log.info("storage file [%s] already exists!" % reco_data_db_file)

        # save filepath
        self.db_file = reco_data_db_file

    def read(self):
        """
        Return content of PKL file.

        Returns
        -------
        pkl_data_dict : dict
            Big dictionnary in PKL file
        """
        log.debug("reading db file [%s]..." % self.db_file)

        # now open pkl file
        with open(self.db_file, 'rb') as f:
            [pkl_data_dict] = pickle.load(f)

        return(pkl_data_dict)

    def get_latest_dataset(self):
        """
        Return the most recent dataset saved.

        Returns
        -------
        dataset : MRSData2 object
            Found dataset
        reco_pipe : pipeline object
            Corresponding pipeline
        """
        log.info("looking for latest dataset in db file [%s]..." % self.db_file)

        # open pkl file
        pkl_data_dict = self.read()

        # for each patient
        ts_diff = timedelta(days=99999999)
        for pnk in list(pkl_data_dict.keys()):
            # each dataset
            for dnk in list(pkl_data_dict[pnk].keys()):
                this_ts = pkl_data_dict[pnk][dnk]["timestamp"]
                if(ts_diff > (datetime.now() - this_ts)):
                    ts_diff = (datetime.now() - this_ts)
                    found_pnk = pnk
                    found_dnk = dnk

        log.info("found [%s/%s]! :)" % (found_pnk, found_dnk))
        latest_data = pkl_data_dict[found_pnk][found_dnk]["data"]
        latest_pipeline = pkl_data_dict[found_pnk][found_dnk]["pipeline"]

        return(latest_data, latest_pipeline)

    def save(self, d, p=None):
        """
        Save MRSData2 object and its reco pipeline and deal with conflicts.

        Parameters
        ----------
        d: MRSData2 object
            Reconstructed data to save
        p: pipeline
            Reco pipeline used to get this data
        """
        log.debug("saving dataset to file [%s]..." % self.db_file)

        # first open pkl file
        pkl_data_dict = self.read()

        # if we reached here, that means the PKL file is not corrupted
        # let's make a backup of it
        copyfile(self.db_file, self.db_file + ".bak")

        # we already have this patient in the db
        if(d.patient_name in pkl_data_dict):
            log.debug("patient name [%s] already exists!" % d.patient_name)
            nd = len(list(pkl_data_dict[d.patient_name].keys()))
            log.debug("already contains %d dataset(s)!" % nd)
        else:
            # create patient entry
            pkl_data_dict[d.patient_name] = {}

        # add/update with the dataset
        log.debug("adding/updating dataset [%s]/[%s]..." % (d.patient_name, d.display_label))
        ts = datetime.now()
        pkl_data_dict[d.patient_name][d.display_label] = {"data": d, "pipeline": p, "timestamp": ts}

        # write back pkl file
        with open(self.db_file, 'wb') as f:
            pickle.dump([pkl_data_dict], f)




class SIEMENS_data_file_reader():
    """A class used to scrap parameters out of DICOM and TWIX files, sometimes in a very dirty way."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("You are trying to dynamically create a new attribute (%s) to this object and that is not cool! I will not let you do that because I believe it is a bad habit and can lead to terrible bugs. A clean way of doing this is to initialize your attribute (%s) in the __init__ method of this class. Bisou, bye :)" % (key, key))
        object.__setattr__(self, key, value)

    def __init__(self, data_fullfilepath):
        """
        Initialize by reading a DICOM or TWIX file in text mode and extracting the DICOM header with pydicom, if required.

        Parameters
        ----------
        data_fullfilepath : string
            Full path to the DICOM/TWIX file
        """
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
        if(self.is_dicom()):
            log.debug("reading DICOM header...")
            self.dcm_header = pydicom.dcmread(self.fullfilepath)
        else:
            self.dcm_header = None

    def is_dicom(self):
        """
        Return true if data is read from a DICOM file.

        Returns
        -------
        self.file_ext == '.dcm' : boolean
            True if dicom file
        """
        return(self.file_ext == '.dcm')

    def read_data(self):
        """
        Read MRS data from file and return a suspect's MRSData object.

        Returns
        -------
        MRSData_obj : MRSData object
            MRS data read from DICOM/TWIX file
        """
        if(self.is_dicom()):
            log.debug("reading DICOM file...")
            MRSData_obj = suspect.io.load_siemens_dicom(self.fullfilepath)
        elif(self.file_ext == '.dat'):
            # try and read this TWIX file
            try:
                log.debug("reading TWIX file...")
                MRSData_obj = suspect.io.load_twix(self.fullfilepath)
            except:
                # well maybe it is broken, maybe the acquisition was interrupted
                # let's try to read it using this modified verion of suspect.io.load_twix
                log.debug("reading broken TWIX file...")
                MRSData_obj = load_broken_twix_vb(self.fullfilepath)
        else:
            log.error("unknown file extension!")

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
        if(self.is_dicom()):
            patient_name_str = str(self.dcm_header.PatientName)
        elif(self.file_ext == '.dat'):
            # find patient name dirty way
            a = self.file_content_str.find("tPatientName")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find("\"", a)
            b = self.file_content_str.find("\"", a + 1)
            patient_name_str = self.file_content_str[(a + 1):b]
            patient_name_str = patient_name_str.strip()
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
        if(self.is_dicom()):
            patient_birthday_datetime = datetime.strptime(str(self.dcm_header.PatientBirthDate), '%Y%m%d')
        elif(self.file_ext == '.dat'):
            # find patient birthday dirty way
            a = self.file_content_str.find("PatientBirthDay")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find("\"", a)
            birthday_str = self.file_content_str[(a + 1):(a + 9)]
            birthday_str = birthday_str.strip()
            if(birthday_str):
                patient_birthday_datetime = datetime.strptime(birthday_str, '%Y%m%d')
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
        if(self.is_dicom()):
            patient_sex_str = str(self.dcm_header.PatientSex)
        elif(self.file_ext == '.dat'):
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
        if(self.is_dicom()):
            patient_weight_kgs = float(self.dcm_header.PatientWeight)
        elif(self.file_ext == '.dat'):
            # find patient weight dirty way
            a = self.file_content_str.find("flUsedPatientWeight")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find(">", a)
            b = self.file_content_str.find("}", a + 1)
            patient_weight_str = self.file_content_str[(a + 5):(b - 2)]
            patient_weight_kgs = float(patient_weight_str.strip())
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
        if(self.is_dicom()):
            patient_height_m = float(self.dcm_header.PatientSize)
        elif(self.file_ext == '.dat'):
            # find patient height dirty way
            a = self.file_content_str.find("flPatientHeight")
            a = self.file_content_str.find("{", a)
            a = self.file_content_str.find(">", a)
            a = self.file_content_str.find(">", a + 1)
            b = self.file_content_str.find("}", a + 1)
            patient_height_str = self.file_content_str[(a + 6):(b - 3)]
            patient_height_m = float(patient_height_str.strip()) / 1000.0
        else:
            log.error("unknown file extension!")

        return(patient_height_m)

    def get_sequence_name(self, file_index=0):
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

    def get_timestamp(self):
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


class MRSData2(suspect.mrsobjects.MRSData):
    """A class based on suspect's MRSData to store MRS data."""

    def __new__(cls, data_filepath, coil_nChannels=8, physio_log_file=None, obj=None, dt=None, f0=None, te=None, ppm0=None, voxel_dimensions=None, transform=None, metadata=None, data_ref=None, label="", offset_display=0.0, timestamp=None, patient_name="", patient_birthday=None, patient_sex=None, patient_weight=None, patient_height=None, tr=None, vref=None, shims=None, sequence_obj=None, noise_level=None, na_pre_data_rejection=None, na_post_data_rejection=None, is_concatenated=None, is_dicom=None):
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
        data_ref : MRSData2 object
            Reference data acquired for this signal
        label : string
            Label for this signal
        offset_display : float
            Y-axis offset display
        timestamp : float
            timestamp in ms
        patient_name : string
            patient name
        patient_birthday : int
            birthyear of patient
        patient_sex : string
            sex of patient ('M', 'F' or 'O')
        patient_weight : float
            patient weight in kgs
        patient_height : float
            patient high in meters
        tr : float
            TR in ms
        vref : float
            reference voltage (V)
        shims : list of floats
            list of shim voltages in volts
        sequence_obj : sim.mrs_sequence object
            sequence object
        noise_level : float
            noise level measured on real FID
        na_pre_data_rejection : int
            number of averages when reading the data file
        na_post_data_rejection : int
            number of averages after data rejection
        is_concatenated : boolean
            was this signal the result of a concatenation?
        is_dicom : boolean
            was this signal read from a DICOM file?

        Returns
        -------
        obj : MRSData2 numpy array [averages,channels,timepoints]
            Resulting constructed MRSData2 object
        """
        if(data_filepath == [] and coil_nChannels == []):
            # calling the parent class' constructor
            pdb.set_trace()
            obj = super(suspect.mrsobjects.MRSData, cls).__new__(cls, obj, dt, f0, te, ppm0, voxel_dimensions, transform, metadata)
            # adding attributes
            obj.data_ref = data_ref
            obj._display_label = label
            obj._display_offset = offset_display
            obj._tr = tr
            obj._timestamp = timestamp
            obj._patient_name = patient_name
            obj._patient_birthday = patient_birthday
            obj._patient_sex = patient_sex
            obj._patient_weight = patient_weight
            obj._patient_height = patient_height
            obj._vref = vref
            obj._shims = shims
            obj._sequence = sequence_obj
            obj._noise_level = noise_level
            obj._na_pre_data_rejection = na_pre_data_rejection
            obj._na_post_data_rejection = na_post_data_rejection
            obj._is_concatenated = is_concatenated
            obj._is_dicom = is_dicom
            # bye
            return(obj)

        # hello
        log.debug("creating object...")

        # open data file
        log.info("reading data file...")
        log.info(data_filepath)
        # read header
        mfr = SIEMENS_data_file_reader(data_filepath)
        # read data and get a suspect MRSData object
        MRSData_obj = mfr.read_data()

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
        obj = super(suspect.mrsobjects.MRSData, cls).__new__(cls, MRSData_obj, MRSData_obj.dt, MRSData_obj.f0, MRSData_obj.te, MRSData_obj.ppm0, MRSData_obj.voxel_dimensions, MRSData_obj.transform, MRSData_obj.metadata)

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

        # TR
        TR_ms = mfr.read_param_num("alTR[0]") / 1000.0
        if(TR_ms is not None):
            log.debug("extracted TR value (%.0fms)" % TR_ms)

        # reference voltage
        vref_v = mfr.read_param_num("flReferenceAmplitude")
        log.debug("extracted reference voltage (%.2fV)" % vref_v)

        # 1st order shim X
        shim_1st_X = mfr.read_param_num("lOffsetX")
        log.debug("extracted 1st order shim value for X (%.2f)" % shim_1st_X)

        # 1st order shim Y
        shim_1st_Y = mfr.read_param_num("lOffsetY")
        log.debug("extracted 1st order shim value for Y (%.2f)" % shim_1st_Y)

        # 1st order shim Z
        shim_1st_Z = mfr.read_param_num("lOffsetZ")
        log.debug("extracted 1st order shim value for Z (%.2f)" % shim_1st_Z)

        # 2nd / 3rd shims
        shims_2nd_3rd_list = []
        for i_shim in range(5):
            this_shim_val = mfr.read_param_num("alShimCurrent[%d]" % i_shim)
            shims_2nd_3rd_list.append(this_shim_val)
        log.debug("extracted 2nd/3rd shims (" + str(shims_2nd_3rd_list) + ")")

        # merge shims values into one vector
        shims_values = [shim_1st_X, shim_1st_Y, shim_1st_Z] + shims_2nd_3rd_list

        # nucleus (used for sequence object)
        nucleus = mfr.get_nucleus()
        log.debug("extracted nuclei (%s)" % nucleus)

        # number of points
        npts = int(mfr.read_param_num("lVectorSize"))
        log.debug("extracted number of points (%d)" % npts)

        # special timestamp
        ulTimeStamp_ms = mfr.get_timestamp()
        log.debug("resulting in a timestamp (%.0fms)" % ulTimeStamp_ms)

        # --- sequence-specific parameters ---
        sequence_name = mfr.get_sequence_name()
        log.debug("extracted sequence name (%s)" % sequence_name)

        if(sequence_name == "eja_svs_slaser"):
            log.debug("this a semi-LASER acquisition, let's extract some specific parameters!")

            # afp pulse (fake) flip angle
            pulse_laser_rfc_fa = mfr.read_param_num("adFlipAngleDegree[1]")
            log.debug("extracted LASER AFP refocussing pulse flip angle (%.2f)" % pulse_laser_rfc_fa)

            # afp pulse length
            pulse_laser_rfc_length = mfr.read_param_num("sWiPMemBlock.alFree[1]")
            log.debug("extracted LASER AFP refocussing pulse duration (%.2f)" % pulse_laser_rfc_length)

            # afp pulse R
            pulse_laser_rfc_r = mfr.read_param_num("sWiPMemBlock.alFree[49]")
            log.debug("extracted LASER AFP refocussing pulse R (%.2f)" % pulse_laser_rfc_r)

            # afp pulse N
            pulse_laser_rfc_n = mfr.read_param_num("sWiPMemBlock.alFree[48]")
            log.debug("extracted LASER AFP refocussing pulse N (%.2f)" % pulse_laser_rfc_n)

            # afp pulse voltage
            pulse_laser_rfc_voltage = mfr.read_param_num("aRFPULSE[1].flAmplitude")
            log.debug("extracted LASER AFP refocussing pulse voltage (%.2f)" % pulse_laser_rfc_voltage)

            # exc pulse duration
            pulse_laser_exc_length = mfr.read_param_num("sWiPMemBlock.alFree[24]")
            log.debug("extracted LASER exc. pulse length (%.2f)" % pulse_laser_exc_length)

            # exc pulse voltage
            pulse_laser_exc_voltage = mfr.read_param_num("aRFPULSE[0].flAmplitude")
            log.debug("extracted LASER excitation pulse voltage (%.2f)" % pulse_laser_exc_voltage)

            # spoiler length
            spoiler_length = mfr.read_param_num("sWiPMemBlock.alFree[12]")
            log.debug("extracted spoiler length (%.2f)" % spoiler_length)

            # build sequence object
            sequence_obj = sim.mrs_seq_eja_svs_slaser(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0, pulse_laser_exc_length / 1000.0, pulse_laser_exc_voltage, pulse_laser_rfc_length / 1000.0, pulse_laser_rfc_fa, pulse_laser_rfc_r, pulse_laser_rfc_n, pulse_laser_rfc_voltage, vref_v, spoiler_length / 1000.0)

        elif(sequence_name == "eja_svs_press"):
            sequence_obj = sim.mrs_seq_eja_svs_press(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)

        elif(sequence_name == "eja_svs_steam"):
            sequence_obj = sim.mrs_seq_eja_svs_steam(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)

        elif(sequence_name == "fid"):
            sequence_obj = sim.mrs_seq_fid(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)

        elif(sequence_name == "svs_st"):
            sequence_obj = sim.mrs_seq_svs_st(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)

        elif(sequence_name == "svs_st_vapor_643"):
            log.debug("this is CMRR's STEAM sequence, let's extract some specific parameters!")

            # TM value
            TM_ms = mfr.read_param_num("alTD[0]") / 1000.0
            log.debug("extracted TM value (%.0fms)" % TM_ms)

            # build sequence object
            sequence_obj = sim.mrs_seq_svs_st_vapor_643(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0, TM_ms)

        elif(sequence_name == "bow_isis_15"):
            # TODO : create a sequence implementation for ISIS?
            sequence_obj = sim.mrs_seq_fid(obj.te, obj.tr, nucleus, npts, 1.0 / obj.dt, obj.f0, 1.0)

        elif(sequence_name is None):
            sequence_obj = None

        else:
            # unknown!
            log.error("ouch unknown sequence!")

        # --- build MRSData2 object ---
        # adding MRSData2 attributes
        obj.data_ref = data_ref
        obj._patient_name = patient_name_str
        obj._patient_birthday = patient_birthday_datetime
        obj._patient_sex = patient_sex_str
        obj._patient_weight = patient_weight_kgs
        obj._patient_height = patient_height_m
        obj._tr = TR_ms
        obj._vref = vref_v
        obj._shims = shims_values
        obj._sequence = sequence_obj
        obj._noise_level = 0.0
        obj._timestamp = ulTimeStamp_ms
        obj._na_pre_data_rejection = None
        obj._na_post_data_rejection = None
        obj._is_concatenated = False
        obj._is_dicom = mfr.is_dicom()

        # those need to be called now, because they the attributes above
        obj.set_display_label()
        obj.set_display_offset()

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
        log.debug("reading physio data...")

        # check file extension
        resp_log_filename_short, resp_log_file_extension = os.path.splitext(
            physio_filepath)
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

    def __array_finalize__(self, obj):
        """
        Overload of special numpy array method called when playing around with stuff relative to object copy etc...

        Parameters
        ----------
        obj : MRSData2 numpy array [1,channels,timepoints]
        """
        super().__array_finalize__(obj)
        self.data_ref = getattr(obj, 'data_ref', None)
        self._resp_trace = getattr(obj, 'resp_trace', None)
        self._display_label = getattr(obj, 'display_label', None)
        self._display_offset = getattr(obj, 'display_offset', 0.0)
        self._patient_birthday = getattr(obj, 'patient_birthday', None)
        self._patient_sex = getattr(obj, 'patient_sex', None)
        self._patient_name = getattr(obj, 'patient_name', None)
        self._patient_weight = getattr(obj, 'patient_weight', None)
        self._patient_height = getattr(obj, 'patient_height', None)
        self._tr = getattr(obj, 'tr', None)
        self._vref = getattr(obj, 'vref', None)
        self._shims = getattr(obj, 'shims', None)
        self._timestamp = getattr(obj, 'timestamp', None)
        self._sequence = getattr(obj, 'sequence', None)
        self._noise_level = getattr(obj, 'noise_level', None)
        self._na_pre_data_rejection = getattr(obj, 'na_pre_data_rejection', None)
        self._na_post_data_rejection = getattr(obj, 'na_post_data_rejection', None)
        self._is_concatenated = getattr(obj, 'is_concatenated', None)
        self._is_dicom = getattr(obj, 'is_dicom', None)

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
        obj2._tr = getattr(self, 'tr', None)
        obj2._display_label = getattr(self, 'display_label', None)
        obj2._display_offset = getattr(self, 'display_offset', 0.0)
        obj2._patient_birthday = getattr(self, 'patient_birthday', None)
        obj2._patient_sex = getattr(self, 'patient_sex', None)
        obj2._resp_trace = getattr(self, 'resp_trace', None)
        obj2._patient_name = getattr(self, 'patient_name', None)
        obj2._patient_weight = getattr(self, 'patient_weight', None)
        obj2._patient_height = getattr(self, 'patient_height', None)
        obj2._vref = getattr(self, 'vref', None)
        obj2._shims = getattr(self, 'shims', None)
        obj2._timestamp = getattr(self, 'timestamp', None)
        obj2._sequence = getattr(self, 'sequence', None)
        obj2._noise_level = getattr(self, 'noise_level', 0.0)
        obj2._na_pre_data_rejection = getattr(self, 'na_pre_data_rejection', 0.0)
        obj2._na_post_data_rejection = getattr(self, 'na_post_data_rejection', 0.0)
        obj2._is_concatenated = getattr(self, 'is_concatenated', 0.0)
        obj2._is_dicom = getattr(self, 'is_dicom', 0.0)
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
        if((lbl == "" or lbl == []) and self.patient_name is not None and self.sequence.name is not None):
            # create a usefull label based on patient name, sequence and object id
            new_lbl = ""
            if(self.patient_name is not None):
                new_lbl = new_lbl + self.patient_name + " | "
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
    def patient_name(self):
        """
        Property get function for patient_name.

        Returns
        -------
        self._patient_name : string
            Name of patient
        """
        return(self._patient_name)

    @property
    def patient_birthday(self):
        """
        Property get function for patient_birthday.

        Returns
        -------
        self._patient_birthday : int
            Birthday of patient
        """
        return(self._patient_birthday)

    @property
    def patient_sex(self):
        """
        Property get function for patient_sex.

        Returns
        -------
        self._patient_sex : string
            sex of patient ('M', 'F' or 'O')
        """
        return(self._patient_sex)

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
    def vref(self):
        """
        Property get function for vref.

        Returns
        -------
        self._vref : float
            Reference voltage
        """
        return(self._vref)

    @property
    def shims(self):
        """
        Property get function for shims.

        Returns
        -------
        self._shims : list of floats
            Shim voltages
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
    def timestamp(self):
        """
        Property get function for timestamp (ms).

        Returns
        -------
        self._timestamp : float
            Timestamp in (ms)
        """
        return(self._timestamp)

    @property
    def na_pre_data_rejection(self):
        """
        Property get function for na_pre_data_rejection.

        Returns
        -------
        self._na_pre_data_rejection : int
            Original number of averages
        """
        return(self._na_pre_data_rejection)

    @property
    def na_post_data_rejection(self):
        """
        Property get function for na_post_data_rejection.

        Returns
        -------
        self._na_post_data_rejection : int
            Number of averages after data rejection
        """
        return(self._na_post_data_rejection)

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
    def is_dicom(self):
        """
        Property get function for is_dicom.

        Returns
        -------
        self._is_dicom : boolean
            True if this current signal was originally read from a DCM file
        """
        return(self._is_dicom)

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

    def correct_intensity_scaling_nd(self, scaling_factor=1e8):
        """
        Amplify the FID signals. Sounds useless but can actually help during quantification! Yes, it is not a good idea to fit signals which have intensities around 1e-6 or lower because of various fit tolerances and also digital problems (epsilon).

        * Works with multi-dimensional signals.

        Parameters
        ----------
        scaling_factor : float
            Amplification factor

        Returns
        -------
        s : MRSData2 numpy array [whatever dimensions]
            Resulting amplified MRSData2 object
        """
        log.debug("intensity scaling [%s]..." % self.display_label)
        log.debug("multiplying time-domain signals by %E..." % scaling_factor)
        # make this louder
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
        s_zf = self.inherit(np.concatenate((s, np.zeros(s_new_shape)), axis=s.ndim - 1))
        log.debug("%d-pts signal + %d zeros = %d-pts zero-filled signal..." % (s.np, nZeros, nPoints_final))

        if(display):
            s_disp = s.copy()
            s_zf_disp = s_zf.copy()
            while(s_disp.ndim > 1):
                s_disp = np.mean(s_disp, axis=0)
                s_zf_disp = np.mean(s_zf_disp, axis=0)

            fig = plt.figure(110)
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

    def correct_phase_3d(self, use_ref_data=True, peak_range=[4.5, 5], weak_ws=False, high_snr=False, phase_order=0, phase_offset=0.0, display=False, display_range=[1, 6]):
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
            if(weak_ws):
                phase_method = 3
            else:
                phase_method = 5

            # and on SNR
            if(high_snr):
                phase_method -= 1

        if(display):
            # prepare subplots
            fig = plt.figure(100)
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
                sf_avg = s[:, c, :].mean(axis=0).spectrum()

                # find maximum peak in range and its phase
                ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
                sf_avg_abs_masked = np.abs(sf_avg)
                sf_avg_abs_masked[ippm_peak_range] = 0
                ippm_peak = np.argmax(sf_avg_abs_masked)
                ppm_peak = ppm[ippm_peak]
                phase_peak_avg = np.angle(sf_avg[ippm_peak])
                if(ippm_peak == 0):
                    log.error("no peak found in specified ppm range or badly phased data!")
                if(c == 0):
                    log.debug("measuring phase at %0.2fppm on 1st channel!" % ppm_peak)

            # late init progress bar
            if(c == 0):
                pbar = log.progressbar("phasing", s.shape[1] * s.shape[0])

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
                        axs[0, 2].axvline(x=ppm[ippm_peak], color='r', linestyle='--')
                    if(phase_method == 4):
                        axs[0, 2].plot(ppm[ippm_peak], np.real(sf[ippm_peak]), 'ro')
                        axs[0, 2].axvline(x=ppm[ippm_peak], color='r', linestyle='--')
                    axs[0, 2].set_xlim(display_range[1], display_range[0])
                    axs[0, 2].set_xlabel('time (s)')
                    axs[0, 2].set_ylabel('original')
                    axs[0, 2].grid('on')
                    axs[0, 2].legend()
                    axs[0, 2].set_title("Meta channel #" + str(c + 1) + " average #" + str(a + 1))

                    # display corrected meta FID
                    axs[1, 1].cla()
                    axs[1, 1].plot(t, np.real(s_phased[a, c, :]), 'b-', linewidth=1, label='real part')
                    axs[1, 1].set_xlabel('time (s)')
                    axs[1, 1].set_ylabel('corrected')
                    axs[1, 1].grid('on')
                    axs[1, 1].legend()

                    # display corrected meta spectrum
                    axs[1, 2].cla()
                    axs[1, 2].plot(ppm, s_phased[a, c, :].spectrum().real, 'b-', linewidth=1, label='real part')
                    axs[1, 2].set_xlim(display_range[1], display_range[0])
                    axs[1, 2].set_xlabel('frequency (ppm)')
                    axs[1, 2].set_ylabel('corrected')
                    axs[1, 2].grid('on')
                    axs[1, 2].legend()

                    fig.subplots_adjust()
                    fig.show()

                pbar.update(c * s.shape[0] + a + 1)

        pbar.finish("done")

        # convert back to MRSData2
        s_phased = self.inherit(s_phased)

        # if any ref data available, we phase it too (silently)
        if(s_phased.data_ref is not None):
            s_phased.data_ref = s_phased.data_ref.correct_phase_3d(True)

        return(s_phased)

    def correct_combine_channels_3d(self, use_ref_data=True, phasing=True, channels_onoff=[True]):
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
        if(self.ndim != 2 or data.ndim != 2):
            log.error("this method only works for 2D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

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

        # first, find peak of interest in range
        ppm = self.frequency_axis_ppm()
        s_avg = np.mean(self, axis=0)
        sf_avg_abs = np.abs(s_avg.spectrum())
        ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
        sf_masked = sf_avg_abs.copy()
        sf_masked[ippm_peak_range] = 0
        ippm_peak_avg = np.argmax(sf_masked, axis=0)
        if(ippm_peak_avg == 0):
            log.error("no peak found in specified ppm range or badly phased data!")
        ppm_peak_avg = ppm[ippm_peak_avg]
        log.debug("found peak of interest at %0.2fppm!" % ppm_peak_avg)

        # for each average in moving averaged data
        peak_trace = np.zeros([self.shape[0], 4])
        pbar = log.progressbar("analyzing", self.shape[0])
        for a in range(0, self.shape[0]):

            # first, measure shift in ppm
            sf_masked_abs = np.abs(self[a, :].spectrum())
            sf_masked_abs[ippm_peak_range] = 0
            ippm_peak = np.argmax(sf_masked_abs)
            if(ippm_peak == 0):
                log.error("no peak found in specified ppm range or badly phased data!")
            ppm_peak = ppm[ippm_peak]
            peak_trace[a, 2] = ppm_peak

            # estimate amplitude
            sf_masked_real = np.real(self[a, :].spectrum())
            sf_masked_real[ippm_peak_range] = 0
            ippm_peak = np.argmax(sf_masked_real)
            if(ippm_peak == 0):
                log.error("no peak found in specified ppm range or badly phased data!")
            ppm_peak = ppm[ippm_peak]
            amp_peak = sf_masked_real[ippm_peak]
            peak_trace[a, 0] = amp_peak

            # estimate linewidth in Hz
            sf_real = np.real(self[a, :].spectrum())
            ippm_max = ippm_peak + np.argmax(sf_real[ippm_peak:] < amp_peak / 2)
            ippm_min = ippm_peak - np.argmax(sf_real[ippm_peak::-1] < amp_peak / 2)
            dppm = np.abs(ppm[ippm_max] - ppm[ippm_min])
            lw = dppm * self.f0
            peak_trace[a, 1] = lw

            # estimate phase in rad
            sf_cmplx = self[a, :].spectrum()
            peak_trace[a, 3] = np.angle(sf_cmplx[ippm_peak])

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
        if(self.resp_trace is None):
            # no physio signal here, exiting
            return()

        # perform peak analysis
        peak_prop_abs, _, _ = self.correct_zerofill_nd()._analyze_peak_2d(peak_range)

        # physio signal
        resp_t = self.resp_trace[0]
        resp_s = self.resp_trace[3]

        # init
        mri_t = np.linspace(self.timestamp, self.timestamp + self.tr * peak_prop_abs.shape[0], peak_prop_abs.shape[0])
        dt_array = np.arange(-delta_time_range / 2.0, delta_time_range / 2.0, 1.0)
        cc_2d = np.zeros([dt_array.shape[0], 4])

        # shift signal and calculate corr coeff
        pbar = log.progressbar("correlating signals", dt_array.shape[0])
        for idt, dt in enumerate(dt_array):
            # build time scale
            this_resp_t_interp = mri_t.copy() + dt
            this_resp_s_interp = np.interp(this_resp_t_interp, resp_t, resp_s)

            # now crop the signals to have the same length
            final_length = min(
                this_resp_s_interp.shape[0], peak_prop_abs.shape[0])
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
        st_ms = self.timestamp
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("data timestamp=\t" + str(st_ms) + "ms\t" + st_str)
        log.info("best start time for...")

        st_ms = self.timestamp + best_dt_per_par[0]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("amplitude=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_per_par[1]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("linewidth=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_per_par[2]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("frequency=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_per_par[3]
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("phase=\t\t" + str(st_ms) + "ms\t" + st_str)
        st_ms = self.timestamp + best_dt_all
        st_str = datetime.fromtimestamp(st_ms / 1000 - 3600).strftime('%H:%M:%S')
        log.info("total=\t\t" + str(st_ms) + "ms\t" + st_str)

        imaxR = np.argmax(best_dt_per_par)
        best_dt = best_dt_per_par[imaxR]
        st_ms = self.timestamp + best_dt
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

    def correct_analyze_and_reject_2d(self, peak_range=[4.5, 5], moving_Naverages=1, peak_properties_ranges={"amplitude (%)": None, "linewidth (Hz)": 30.0, "chemical shift (ppm)": 0.5, "phase std. factor (%)": 60.0}, peak_properties_rel2mean=True, auto_adjust_lw_bound=False, auto_adjust_allowed_snr_change=-10.0, display=False):
        """
        Analyze peak in each average in terms intensity, linewidth, chemical shift and phase and reject data if one of these parameters goes out of the min / max bounds. Usefull to understand what the hell went wrong during your acquisition when you have the raw data (TWIX) and to try to improve things a little...

        * Works only with a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : list
            Range in ppm used to analyze peak phase when no reference signal is specified
        moving_Naverages : int
            Number of averages to perform when using moving average, need to be an odd number
        peak_properties_ranges : dict
            Dictionnary that contains 4 entries, 4 rejection criterias for
                "amplitude (%)": amplitude relative changes: keep data if within +/-val % range
                "linewidth (Hz)": linewidth changes: keep data is below val Hz
                "chemical shift (ppm)": chemical shift changes: keep data is within +/-val ppm
                "phase std. factor (%)": phase changes: keep data if within +/- val/100 * std(phase) rd
        peak_properties_rel2mean : boolean
            Relative peak properties (amplitude, chemical shift and phase) should be caculated based on the mean value over the whole acquisition (True) or only the first acquired point (False)
        auto_adjust_lw_bound : boolean
            Try to adjust the peak linewidth rejection criteria automatically
        auto_adjust_allowed_snr_change : float
            Allowed change in SNR (%), a positive or negative relative to the initial SNR without data rejection
        display : boolean
            Display correction process (True) or not (False)

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
        s._na_pre_data_rejection = s.shape[0]
        if(s.shape[0] == 1):
            log.warning("single-shot signal, nothing to analyze!")
            return(s)
        ppm = s.frequency_axis_ppm()

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
        peak_prop_analyze[:, 3] = peak_prop_abs[:, 3]

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

        # special for phase: rejection range is a factor of std
        phase_std = peak_prop_analyze[:, 3].std()
        phase_std_reject_range = peak_properties_ranges_list[3] / 100.0 * phase_std

        # prepare rejection ranges
        peak_prop_min = [-peak_properties_ranges_list[0],
                         0.0,
                         -peak_properties_ranges_list[2],
                         -phase_std_reject_range]

        peak_prop_max = [+peak_properties_ranges_list[0],
                         peak_properties_ranges_list[1],
                         +peak_properties_ranges_list[2],
                         +phase_std_reject_range]

        if(auto_adjust_lw_bound):
            # automatic rejection based on lw

            # first find optimal lw sweeping range
            lw_min = (np.floor(peak_prop_analyze[:, 1].min() / 10.0)) * 10.0
            lw_max = (np.floor(peak_prop_analyze[:, 1].max() / 10.0) + 1.0) * 10.0
            lw_range = np.arange(lw_min, lw_max, 1.0)

            # iterate between max and min for linewidth, and test the resulting data
            pbar = log.progressbar("adjusting linewidth threshold", lw_range.shape[0])

            test_snr_list = np.zeros(lw_range.shape)
            test_lw_list = np.zeros(lw_range.shape)
            test_nrej_list = np.zeros(lw_range.shape)
            # test each lw
            for (ilw, this_lw) in enumerate(lw_range):
                # rebuild min/max rejection bounds
                peak_prop_min_auto = peak_prop_min.copy()
                peak_prop_min_auto[1] = 0.0
                peak_prop_max_auto = peak_prop_max.copy()
                peak_prop_max_auto[1] = this_lw

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
                    old_level = log.getLevel()
                    log.setLevel(log.ERROR)
                    test_snr_list[ilw], _, _ = this_s_cor.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_1d().analyze_snr_1d(peak_range)
                    test_lw_list[ilw] = this_s_cor.correct_zerofill_nd().correct_realign_2d().correct_average_2d().correct_apodization_1d().analyze_linewidth_1d(peak_range)
                    log.setLevel(old_level)
                    test_nrej_list[ilw] = this_mask_reject_data_sumup.sum()

                # progression
                pbar.update(ilw)

            pbar.finish("done")

            # analyze SNR curve choose LW threshold
            snr_initial = test_snr_list[-1]
            snr_threshold = snr_initial + snr_initial * auto_adjust_allowed_snr_change / 100.0
            test_snr_list_rel = test_snr_list / snr_initial * 100.0 - 100.0
            test_snr_list_mask = (test_snr_list_rel > auto_adjust_allowed_snr_change)

            if(not test_snr_list_mask.any()):
                # that was a bit ambitious
                log.warning("sorry but your exceptation regarding the SNR was a bit ambitious! You are refusing to go under %.0f%% SNR change. While trying to adjust the linewidth criteria for data rejection, the best we found was a %.0f%% SNR change :(" % (auto_adjust_allowed_snr_change, test_snr_list_rel.max()))
                # set optimal LW to max
                optim_lw = lw_max
            else:
                # we found SNR values which matches our request
                # let's choose the one with the lowest LW
                ind_lowest_lw = np.argmax(test_snr_list_rel > auto_adjust_allowed_snr_change)
                optim_lw = lw_range[ind_lowest_lw]
                log.info("found optimal linewidth for data rejection = %.1f Hz" % optim_lw)
                # adjusting bounds
                peak_prop_max[1] = optim_lw

            # plot SNR / LW combinaisons and optimal choice
            if(display):
                fig = plt.figure(130)
                fig.clf()
                axs = fig.subplots(1, 2, sharex='all')
                ax2 = axs[0].twinx()
                ax3 = axs[1].twinx()
                fig.canvas.set_window_title("mrs.reco.MRSData2.correct_analyze_and_reject_2d (auto)")
                fig.suptitle("adjusting linewidth data rejection criteria for [%s]" % self.display_label)

                axs[0].plot(lw_range, test_snr_list, 'rx-', label='SNR')
                axs[0].plot([optim_lw, optim_lw], [test_snr_list.min(), test_snr_list.max()], 'm--', label='Optimal linewidth')
                axs[0].plot([lw_range.min(), lw_range.max()], [snr_threshold, snr_threshold], 'g--', label='SNR threshold')
                axs[0].set_xlabel('Max allowed linewidth (Hz)')
                axs[0].set_ylabel('Estimated SNR (u.a)')
                axs[0].grid('on')
                axs[0].legend(loc='lower left')

                ax2.plot(lw_range, test_lw_list, 'bx-', label='Linewidth')
                ax2.set_ylabel('Estimated linewidth (Hz)')
                ax2.legend(loc='lower right')

                axs[1].plot(lw_range, test_nrej_list, 'ko-', label='Number of scans rejected')
                axs[1].set_xlabel('Max allowed linewidth (Hz)')
                axs[1].set_ylabel('Number of scans rejected')
                axs[1].grid('on')

                ax3.plot(lw_range, test_nrej_list / s_ma.shape[0] * 100, 'ko-', label='Total percentage of scans rejected')
                ax3.set_ylabel('Estimated linewidth (Hz)')
                ax3.set_ylabel('Rejection percentage (%)')

                fig.subplots_adjust()
                fig.show()

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
        log.info("data rejection: summary")
        log.info("number of averages rejected because of...")
        log.info("amplitude = %d" % mask_reject_data[:, 0].sum())
        log.info("linewidth = %d" % mask_reject_data[:, 1].sum())
        log.info("frequency = %d" % mask_reject_data[:, 2].sum())
        log.info("phase = %d" % mask_reject_data[:, 3].sum())

        # actually reject data now
        mask_reject_data_sumup = (mask_reject_data.sum(axis=1) > 0)
        s_cor = s[(mask_reject_data_sumup == False), :]

        log.info("TOTAL data rejection = %d / %d (%.0f%%)" % (mask_reject_data_sumup.sum(), s_ma.shape[0], (mask_reject_data_sumup.sum() / s_ma.shape[0] * 100)))

        # final display
        if(display):

            fig = plt.figure(131)
            fig.clf()
            axs = fig.subplots(2, 3, sharex='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_analyze_and_reject_2d")
            fig.suptitle("analyzing data and rejecting some for [%s]" % self.display_label)

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

                # build lineshape segment
                sf_analyze = np.real(s_ma[k, :].spectrum())
                sf_analyze_cropped = sf_analyze.copy()
                sf_analyze_cropped[ippm_peak_range] = 0
                ippm_peak = np.argmax(sf_analyze_cropped)
                if(ippm_peak == 0):
                    log.error("no peak found in specified ppm range or badly phased data!")

                # estimate linewidth in Hz
                amp_peak = sf_analyze[ippm_peak]
                ippm_max = ippm_peak + np.argmax(sf_analyze[ippm_peak:] < amp_peak / 2)
                ippm_min = ippm_peak - np.argmax(sf_analyze[ippm_peak::-1] < amp_peak / 2)
                ippm_half_peak = np.arange(ippm_min, ippm_max)
                plt.plot(ppm[ippm_half_peak], sf_analyze[ippm_half_peak] * ampfactor + ystep * k, 'k-', linewidth=1)

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

        # keep the amount of data we reject in mind
        s_cor._na_post_data_rejection = s_cor.shape[0]

        return(s_cor)

    def correct_realign_2d(self, peak_range=[4.5, 5], moving_Naverages=1, display=False, display_range=[1, 6]):
        """
        Realign each signal of interest in frequency by taking as a reference the first spectra in absolute mode.

        * Works only with a 2D [averages,timepoints] signal.

        Parameters
        ----------
        peak_range : list [2]
            Range in ppm used to analyze peak phase
        moving_Naverages : int
            Number of averages to perform when using moving average, need to be an odd number
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
            log.debug("measuring peak properties at %0.2fppm!" % ppm_peak_avg)

            # for each average in moving averaged data
            s_realigned_ma = s_ma.copy()
            df_trace = np.zeros(s_ma.shape[0])
            pbar = log.progressbar("realigning", s_ma.shape[0])
            for a in range(0, s_ma.shape[0]):

                # measure shift on moving average data
                sf_masked = np.abs(s_ma[a, :].spectrum())
                sf_masked[ippm_peak_range] = 0
                ippm_peak = np.argmax(sf_masked)
                if(ippm_peak == 0):
                    log.error("no peak found in specified ppm range or badly phased data!")
                ppm_peak = ppm[ippm_peak]

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
                axs[0, 0].plot(ppm, sf_avg_abs, 'k-', linewidth=1)
                axs[0, 0].set_xlim(display_range[1], display_range[0])
                axs[0, 0].set_ylabel('averaged')
                axs[0, 0].grid('on')
                # add peak position
                axs[0, 0].plot(ppm_peak_avg, sf_avg_abs[ippm_peak_avg], 'ro')
                axs[0, 0].axvline(x=ppm_peak_avg, color='r', linestyle='--')

                # display original data
                axs[0, 1].plot(ppm, np.abs(s_ma.spectrum().transpose()), 'k-', linewidth=1)
                axs[0, 1].set_xlim(display_range[1], display_range[0])
                axs[0, 1].set_ylabel('original')
                axs[0, 1].grid('on')

                # display corrected spectrum
                axs[1, 0].plot(s_realigned_ma.frequency_axis_ppm(), np.abs(s_realigned_ma.spectrum().transpose()), 'b-', linewidth=1)
                axs[1, 0].set_xlim(display_range[1], display_range[0])
                axs[1, 0].set_xlabel('chemical shift (ppm)')
                axs[1, 0].set_ylabel('corrected')
                axs[1, 0].grid('on')

                # display corrected averaged spectrum
                axs[1, 1].plot(s_realigned_ma.frequency_axis_ppm(), np.abs(np.mean(s_realigned_ma, axis=0).spectrum().transpose()), 'b-', linewidth=1)
                axs[1, 1].set_xlim(display_range[1], display_range[0])
                axs[1, 1].set_xlabel('chemical shift (ppm)')
                axs[1, 1].set_ylabel('averaged & corrected')
                axs[1, 1].grid('on')

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

    def analyze_noise_1d(self, n_pts=100):
        """
        Measure noise level in time domain and store it in the "noise_level" attribute. This is usefull to keep track of the original noise level for later use, CRB normalization durnig quantification for example.

        * Works only with a 1D [timepoints] signal.

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
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        log.debug("estimating noise level in FID using last %d points..." % n_pts)
        s = self.copy()
        s_real = np.real(s)
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
            s.data_ref = s.data_ref.analyze_noise_1d(n_pts)

        return(s)

    def correct_apodization_1d(self, apo_factor=1.0, display=False, display_range=[1, 6]):
        """
        Apodize signal using an exponential window adjusted by a linewidth parameter in Hz. Works only for a 1D MRSData2 object.

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
        s_apo : MRSData2 numpy array [timepoints]
            Resulting apodized data stored in a MRSData2 object
        """
        log.debug("apodizing [%s]..." % self.display_label)
        # dimensions check
        if(self.ndim != 1):
            log.error("this method only works for 1D signals! You are feeding it with %d-dimensional data. :s" % self.ndim)

        # init
        s = self.copy()
        t = s.time_axis()
        w_apo = np.exp(-apo_factor * t)
        s_apo = s * w_apo

        if(display):
            ppm = s.frequency_axis_ppm()
            ppm_apo = s_apo.frequency_axis_ppm()

            fig = plt.figure(160)
            fig.clf()
            axs = fig.subplots(2, 2, sharex='row', sharey='row')
            fig.canvas.set_window_title("mrs.reco.MRSData2.correct_apodization")
            fig.suptitle("apodizing [%s]" % self.display_label)

            axs[0, 0].plot(t, np.abs(s), 'k-', linewidth=1, label='fid')
            axs[0, 0].plot(t, w_apo * np.abs(s.max()), 'r-', linewidth=1, label='apodization window')
            axs[0, 0].set_xlabel('time (s)')
            axs[0, 0].set_ylabel('original')
            axs[0, 0].grid('on')

            axs[0, 1].plot(t, np.abs(s_apo), 'b-', linewidth=1)
            axs[0, 1].set_xlabel('time (s)')
            axs[0, 1].set_ylabel('apodized')
            axs[0, 1].grid('on')

            axs[1, 0].plot(ppm, s.spectrum().real, 'k-', linewidth=1)
            axs[1, 0].set_xlabel('chemical shift (ppm)')
            axs[1, 0].set_ylabel('original spectrum')
            axs[1, 0].set_xlim(display_range[1], display_range[0])
            axs[1, 0].grid('on')

            axs[1, 1].plot(ppm_apo, s_apo.spectrum().real, 'b-', linewidth=1)
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
            s_apo.data_ref = s_apo.data_ref.correct_apodization_1d(apo_factor, False)

        return(s_apo)

    def correct_crop_1d(self, nPoints_final=4096, display=False, display_range=[1, 6]):
        """
        Crop signal in time-domain to remove last points.

        * Works only with a 1D [timepoints] signal.

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
            if(self.sequence is not None):
                log.debug("updating sequence.npts...")
                self.sequence.npts = nPoints_final
                self.sequence._ready = False
        else:
            s_crop = self.copy()
            log.debug("no cropping needed, getting bored...")

        if(display):
            t = s.time_axis()
            t_crop = s_crop.time_axis()
            ppm = s.frequency_axis_ppm()
            ppm_crop = s_crop.frequency_axis_ppm()

            fig = plt.figure(170)
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
            fig = plt.figure(180)
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
        ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
        sf_abs_masked = np.real(s.spectrum())
        sf_abs_masked[ippm_peak_range] = 0
        ippm_peak = np.argmax(sf_abs_masked)
        if(ippm_peak == 0):
            log.error("no peak found in specified ppm range or badly phased data!")
        ppm_peak = ppm[ippm_peak]
        log.debug("peak detected at %0.2fppm -> %0.2fppm!" % (ppm_peak, peak_real_ppm))

        # estimate frequency shift in Hz
        log.debug("frequency shifting data...")
        dppm = (peak_real_ppm - ppm_peak)
        df = dppm * s.f0
        s_shifted = s.adjust_frequency(-df)

        if(display):
            fig = plt.figure(190)
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
            axs[0].plot(ppm_peak, s.spectrum()[ippm_peak].real, 'ro')
            axs[0].axvline(x=ppm_peak, color='r', linestyle='--')

            axs[1].plot(ppm, s_shifted.spectrum().real, 'b-', linewidth=1)
            axs[1].set_xlim(display_range[1], display_range[0])
            axs[1].set_xlabel('chemical shift (ppm)')
            axs[1].set_ylabel('shifted')
            axs[1].grid('on')
            # add new peak position
            axs[1].plot(peak_real_ppm, s.spectrum()[ippm_peak].real, 'ro')
            axs[1].axvline(x=peak_real_ppm, color='r', linestyle='--')

            fig.subplots_adjust()
            fig.show()

        # convert back to MRSData2
        s_shifted = self.inherit(s_shifted)

        return(s_shifted)

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
            fig = plt.figure(200)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_snr_1d")
            fig.suptitle("analyzing SNR for [%s]" % self.display_label)

        # find maximum peak in range and its chemical shift
        ppm = s.frequency_axis_ppm()
        sf = s.spectrum()
        if(magnitude_mode):
            sf_analyze = np.abs(sf)
            log.debug("going to analyze the MAGNITUDE spectrum...")
        else:
            sf_analyze = np.real(sf)
            log.debug("going to analyze the REAL spectrum...")

        ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
        sf_analyze2 = sf_analyze.copy()
        sf_analyze2[ippm_peak_range] = 0
        ippm_peak = np.argmax(sf_analyze2)
        if(ippm_peak == 0):
            log.warning("no peak found in specified ppm range or badly phased data!")
            return(np.nan, np.nan, np.nan)

        ppm_peak = ppm[ippm_peak]
        log.debug("by measuring the intensity at %0.2fppm!" % ppm_peak)
        snr_signal = sf_analyze[ippm_peak]

        # estimate noise in user specified spectral region
        log.debug("estimating noise from %0.2f to %0.2fppm region!" % (noise_range[0], noise_range[1]))
        ippm_noise_range = (noise_range[0] < ppm) & (ppm < noise_range[1])
        snr_noise = np.std(sf_analyze[ippm_noise_range])

        if(display):
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
            ax.plot(ppm[ippm_peak], sf_analyze[ippm_peak], 'ro')
            ax.axvline(x=ppm[ippm_peak], color='r', linestyle='--')
            # show noise region
            ax.plot(ppm[ippm_noise_range], sf_analyze[ippm_noise_range], 'bo')

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
        ppm = s.frequency_axis_ppm()

        # find maximum peak in range and its chemical shift
        sf = s.spectrum()
        if(magnitude_mode):
            sf_analyze = np.abs(sf)
            log.debug("going to analyze the MAGNITUDE spectrum...")
        else:
            sf_analyze = np.real(sf)
            log.debug("going to analyze the REAL spectrum...")

        # find peak in range
        ippm_peak_range = (peak_range[0] > ppm) | (ppm > peak_range[1])
        sf_analyze_cropped = sf_analyze.copy()
        sf_analyze_cropped[ippm_peak_range] = 0
        ippm_peak = np.argmax(sf_analyze_cropped)
        if(ippm_peak == 0):
            log.error("no peak found in specified ppm range or badly phased data!")
        ppm_peak = ppm[ippm_peak]
        log.debug("estimating the peak linewidth at %0.2fppm!" % ppm_peak)

        # estimate linewidth in Hz
        amp_peak = sf_analyze[ippm_peak]
        ippm_max = ippm_peak + np.argmax(sf_analyze[ippm_peak:] < amp_peak / 2)
        ippm_min = ippm_peak - np.argmax(sf_analyze[ippm_peak::-1] < amp_peak / 2)
        ippm_half_peak = np.arange(ippm_min, ippm_max)
        dppm = np.abs(ppm[ippm_max] - ppm[ippm_min])
        lw = dppm * s.f0
        log.info("results for [" + s.display_label + "] coming...")
        log.info("LW = %0.2f Hz!" % lw)

        if(display):
            fig = plt.figure(210)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='all', sharey='all')
            fig.canvas.set_window_title("mrs.reco.MRSData2.analyze_linewidth_1d")
            fig.suptitle("analyzing peak linewidth for [%s]" % self.display_label)

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
            ax.plot(ppm[ippm_half_peak], sf_analyze[ippm_half_peak], 'r-')

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
        subj.patientName = self.patient_name
        subj.patientID = self.patient_name
        # bug here...
        # subj.patientBirthdate = self.patient_birthday.strftime('%Y-%m-%d')
        if(self.patient_sex == 0):
            subj.patientGender = "M"
        elif(self.patient_sex == 1):
            subj.patientGender = "F"
        elif(self.patient_sex == 2):
            subj.patientGender = "O"
        else:
            log.error("patient gender unknown!")
        subj.patientWeight_kg = self.patient_weight
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
            raise TypeError("You are trying to dynamically create a new attribute (%s) to this object and that is not cool! I will not let you do that because I believe it is a bad habit and can lead to terrible bugs. A clean way of doing this is to initialize your attribute (%s) in the __init__ method of this class. Bisou, bye :)" % (key, key))
        object.__setattr__(self, key, value)

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

        # option to process only a set of datasets: list of indexes
        self.data_process_only_this_data_index = []

        # --- physio data ---
        # list of filepaths pointing to physio files
        self.data_physio_filepaths = ""
        # parsed list of physio filepaths
        self._data_physio_list = []

        # --- display legend captions ---
        # this can be either one big string with captions separated by line breaks (that will be parsed) OR simply a list
        self.display_legends = ""
        # internal list of legend captions
        self._display_legends_list = []

        # --- display options ---
        self.display_offset = 0.0
        # set reference ppm
        self.ppm0 = 4.7

        # --- available jobs and their parameters ---
        self.jobs = {}
        # start with display because other jobs depend on some display parameters
        # --- job: spectrum final display ---
        self.jobs["displaying"] = {0:
                                   {"func": MRSData2.display_spectrum_1d, "name": "displaying"},
                                   # figure index
                                   "fig_index": 1,
                                   # ppm range used for display
                                   "range_ppm": [1, 6],
                                   # display spectrum in magnitude mode?
                                   "magnitude_mode": False
                                   }

        # --- job: automatic rephasing ---
        self.jobs["phasing"] = {0:
                                {"func": MRSData2.correct_phase_3d, "name": "phasing"},
                                # use reference data is available?
                                "using_ref_data": True,
                                # ppm range to look fo peak used to estimate phase
                                "POI_range_ppm": [1.5, 2.5],
                                # weak water suppression was used, thus we can use the 1st pt of the FID
                                "weak_ws_mode": False,
                                # high SNR? we rephase each individual spectrum from each channel separately
                                "high_snr_mode": False,
                                # order of phasing in time: 0th or 1st order
                                "order": 0,
                                # add an additional 0th order phase (rd)
                                "offset": 0.0,
                                # display all this process to check what the hell is going on
                                "display": False,
                                "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                }

        # --- job: amplification ---
        self.jobs["scaling"] = {0:
                                {"func": MRSData2.correct_intensity_scaling_nd, "name": "scaling intensity"},
                                "scaling_factor": 1e8
                                }

        # --- job: FID modulus ---
        self.jobs["FID modulus"] = {0:
                                    {"func": MRSData2.correct_fidmodulus_nd, "name": "FID modulus"}
                                    }

        # --- job: channel combination ---
        self.jobs["channel-combining"] = {0:
                                          {"func": MRSData2.correct_combine_channels_3d, "name": "channel-combining"},
                                          # use non water-suppressed data to recombine and rephase channels
                                          "using_ref_data": True,
                                          # should we rephase (0th order) data while combining?
                                          "phasing": True,
                                          # boolean mask to switch on/off some Rx channels
                                          "weights": [True]
                                          }

        # --- job: concatenate ---
        self.jobs["concatenate"] = {0:
                                    {"func": MRSData2.concatenate_2d, "name": "concatenate"}
                                    }

        # --- job: zero-filling ---
        self.jobs["zero-filling"] = {0:
                                     {"func": MRSData2.correct_zerofill_nd, "name": "zero-filling"},
                                     # number of signal points after zf
                                     "npts": 8192 * 2,
                                     # display all this process to check what the hell is going on
                                     "display": True,
                                     "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                     }

        # --- job: analyze physio signal ---
        self.jobs["physio-analysis"] = {0:
                                        {"func": MRSData2.analyze_physio_2d, "name": "analyzing physio. signals"},
                                        # ppm range to look for a peak to analyze
                                        "POI_range_ppm": [4.5, 5],
                                        # time range in (ms) to look around timestamp for correlation physio/MRS
                                        "delta_time_ms": 1000.0,
                                        # display all this process to check what the hell is going on
                                        "display": True
                                        }

        # --- job: automatic data rejection based on criterias ---
        self.jobs["data-rejecting"] = {0:
                                       {"func": MRSData2.correct_analyze_and_reject_2d, "name": "data rejecting"},
                                       # ppm range to look for a peak to analyze
                                       "POI_range_ppm": [4.5, 5],
                                       # size of moving average window
                                       "moving_averages": 1,
                                       # rejection criterias for
                                       # amplitude relative changes: keep data if within +/-val % range
                                       # linewidth changes: keep data is below val Hz
                                       # chemical shift changes: keep data is within +/-val ppm
                                       # phase changes: keep data if within +/-val/100 * std(phase) rd
                                       "ranges": {"amplitude (%)": None,
                                                  "linewidth (Hz)": 30.0,
                                                  "chemical shift (ppm)": 0.5,
                                                  "phase std. factor (%)": None},
                                       # for amplitude, chemical shift and phase, the rejection of data is based on ranges of relative changes of those metrics. Relative to what? The mean value over the whole acquisition (True) or the first acquired point (False)
                                       "rel2mean": True,
                                       # automatic adjustement of linewidth criteria
                                       "auto": False,
                                       # minimum allowed SNR change (%) when adjusting the linewidth criteria, this can be positive (we want to increase SNR +10% by rejecting crappy data) or negative (we are ok in decreasing the SNR -10% in order to get better resolved spectra)
                                       "auto_allowed_snr_change": -10.0,
                                       # display all this process to check what the hell is going on
                                       "display": True
                                       }

        # --- job: automatic data frequency realignment ---
        self.jobs["realigning"] = {0:
                                   {"func": MRSData2.correct_realign_2d, "name": "frequency realigning"},
                                   # ppm range to look for a peak to analyze
                                   "POI_range_ppm": [4.5, 5],
                                   # size of moving average window
                                   "moving_averages": 1,
                                   # display all this process to check what the hell is going on
                                   "display": True,
                                   "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                   }

        # --- job: data averaging ---
        self.jobs["averaging"] = {0:
                                  {"func": MRSData2.correct_average_2d, "name": "averaging"},
                                  # number of averages to mean (None = all)
                                  "na": None,
                                  # display all this process to check what the hell is going on
                                  "display": True,
                                  "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                  }

        # --- job: noise level analysis ---
        self.jobs["noise-estimation"] = {0:
                                         {"func": MRSData2.analyze_noise_1d, "name": "estimating noise level"},
                                         # estimate noise std time-domain on the last 100 pts of the FID
                                         "npts": 100
                                         }

        # --- job: data apodization ---
        self.jobs["apodizing"] = {0:
                                  {"func": MRSData2.correct_apodization_1d, "name": "apodizing"},
                                  # exponential damping factor for apodization (Hz)
                                  "damping_hz": 5,
                                  # display all this process to check what the hell is going on
                                  "display": True,
                                  "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                  }

        # --- job: data cropping ---
        self.jobs["cropping"] = {0:
                                 {"func": MRSData2.correct_crop_1d, "name": "cropping"},
                                 # final number of signal points after crop
                                 "final_npts": 4096,
                                 # display all this process to check what the hell is going on
                                 "display": True,
                                 "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                 }

        # --- job: water post-acquisition removal ---
        self.jobs["water-removal"] = {0:
                                      {"func": MRSData2.correct_water_removal_1d, "name": "removing water peak"},
                                      # number of components when running HSVD
                                      "hsvd_components": 5,
                                      # ppm range where all components will be remove
                                      "hsvd_range": [4.6, 4.8],
                                      # display all this process to check what the hell is going on
                                      "display": True,
                                      "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                      }

        # --- job: spectrum chemical shift calibration ---
        self.jobs["calibrating"] = {0:
                                    {"func": MRSData2.correct_freqshift_1d, "name": "frequency shifting"},
                                    # ppm range to look for the peak of interest (NAA by default)
                                    "POI_range_ppm": [1.8, 2.3],
                                    # real ppm value for this peak
                                    "POI_true_ppm": 2.008,
                                    # display all this process to check what the hell is going on
                                    "display": True,
                                    "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                    }

        # --- job: SNR analysis ---
        self.jobs["analyzing-snr"] = {0:
                                      {"func": MRSData2.analyze_snr_1d, "name": "analyzing SNR"},
                                      # ppm range to look for a peak to analyze
                                      "s_range_ppm": [1.8, 2.1],
                                      # ppm range to look for pure noise
                                      "n_range_ppm": [-1, 0],
                                      # should we look at the magnitude or real spectrum?
                                      "magnitude_mode": False,
                                      # display all this process to check what the hell is going on
                                      "display": True,
                                      "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                      }

        # --- job: LW analysis ---
        self.jobs["analyzing-lw"] = {0:
                                     {"func": MRSData2.analyze_linewidth_1d, "name": "analyzing peak-linewidth"},
                                     # ppm range to look for a peak to analyze
                                     "range_ppm": [4.5, 5],
                                     # should we look at the magnitude or real spectrum?
                                     "magnitude_mode": False,
                                     # display all this process to check what the hell is going on
                                     "display": True,
                                     "display_range_ppm": self.jobs["displaying"]["range_ppm"]
                                     }

        # --- job list ---
        # list of data processing to apply to the data
        # beware, you need to know what you are doing here
        # also, be careful with the dimensionality of data 3D, 2D, 1D along the data processing
        # order is important!
        self.job_list = [self.jobs["phasing"],
                         self.jobs["scaling"],
                         self.jobs["FID modulus"],
                         self.jobs["channel-combining"],
                         self.jobs["concatenate"],
                         self.jobs["zero-filling"],
                         self.jobs["physio-analysis"],
                         self.jobs["data-rejecting"],
                         self.jobs["realigning"],
                         self.jobs["averaging"],
                         self.jobs["noise-estimation"],
                         self.jobs["apodizing"],
                         self.jobs["cropping"],
                         self.jobs["water-removal"],
                         self.jobs["calibrating"],
                         self.jobs["displaying"]]

        # --- analyze job list ---
        # SNR/LW analysis job list
        self.analyze_job_list = [self.jobs["channel-combining"],
                                 self.jobs["zero-filling"],
                                 self.jobs["averaging"],
                                 self.jobs["calibrating"]]

        # --- SNR/LW analysis ---
        self.analyze_enable = True
        self.analyze_display = True
        # list of measured SNR/LW
        self._analyze_results_dict = {}

        # freeze the object and prevent the creation of new attributes
        self.__isfrozen = True

    def _run_job(self, job, data, default_args=False):
        """
        Estimate SNR and/or peak linewidth for this dataset. Values are stored. A mini default pipeline is applied before SNR/LW measurements and can be set with self.analyze_job_list.

        Parameters
        ----------
        job : dict entry from self.jobs
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
        data_snr, _, _ = self._run_job(self.jobs["analyzing-snr"], data)
        data_lw = self._run_job(self.jobs["analyzing-lw"], data)

        # allow outputs
        log.setLevel(old_level)

        # first time we report this dataset: init result dict
        if(data.display_label not in list(self._analyze_results_dict.keys())):
            self._analyze_results_dict[data.display_label] = {"snr": {}, "lw": {}}

        # output and store SNR
        job_label = "post-" + current_job[0]["name"]
        log.info(job_label + " SNR of [%s] = %.2f" % (data.display_label, data_snr))
        self._analyze_results_dict[data.display_label]["snr"][job_label] = data_snr

        # output and store lw
        log.info(job_label + " LW of [%s] = %.2f" % (data.display_label, data_lw))
        self._analyze_results_dict[data.display_label]["lw"][job_label] = data_lw

    def _display_analyze_results(self):
        """Print final SNR and peak-linewidth for each dataset. Plot bargraph showing evolution of SNR and linewidth during data processing (to check that a job didd not destroy the data for example!)."""
        log.info("displaying SNR and linewidth final results...")

        # terminal output
        # first find longest data label
        data_labels = list(self._analyze_results_dict.keys())
        data_label_nchar = len(max(data_labels, key=len))

        # for each dataset (line), print final SNR/LW columns
        log.info_line________________________()
        print("dataset".ljust(data_label_nchar) + "\t" + "SNR (u.a)".ljust(10) + "\t" + "LW (Hz)".ljust(10), flush=True)
        print("", flush=True)
        first_data_key = list(self._analyze_results_dict.keys())[0]
        job_labels = list(self._analyze_results_dict[first_data_key]["snr"].keys())
        last_job_key = job_labels[-1]
        for d in data_labels:
            data_label_padded = d.ljust(data_label_nchar)
            d_snr = self._analyze_results_dict[d]["snr"][last_job_key]
            d_snr_str = ("%.2f" % d_snr).ljust(10)
            d_lw = self._analyze_results_dict[d]["lw"][last_job_key]
            d_lw_str = ("%.2f" % d_lw).ljust(10)
            print("%s\t%s\t%s" % (data_label_padded, d_snr_str, d_lw_str), flush=True)
        log.info_line________________________()

        # display SNR/LW evolution for all data and all jobs
        if(self.analyze_display):
            fig = plt.figure(220)
            fig.clf()
            axs = fig.subplots(2, 1, sharex='row')
            fig.canvas.set_window_title("mrs.reco.pipeline._display_analyze_results")
            fig.suptitle("SNR/peak-linewidth results")

            # prepare bars
            nBars = len(data_labels)
            width = max(-0.15 * nBars + 0.6, 0.1)  # bars get thinner if more data
            pos_bars = np.arange(len(job_labels))
            pos_shift = np.linspace(-width * (nBars - 1) / 2.0, +width * (nBars - 1) / 2.0, nBars)
            snr_ylim = [+np.inf, -np.inf]
            lw_ylim = [+np.inf, -np.inf]

            # iterate in lists and plot bars
            for d, pos in zip(data_labels, pos_shift):
                # snr list for this dataset
                this_snr_list = list(self._analyze_results_dict[d]["snr"].values())
                snr_ylim[0] = min(snr_ylim[0], np.min(this_snr_list))
                snr_ylim[1] = max(snr_ylim[1], np.max(this_snr_list))
                axs[0].bar(pos_bars + pos, this_snr_list, width, label=d)
                # lw list for this dataset
                this_lw_list = list(self._analyze_results_dict[d]["lw"].values())
                lw_ylim[0] = min(lw_ylim[0], np.min(this_lw_list))
                lw_ylim[1] = max(lw_ylim[1], np.max(this_lw_list))
                axs[1].bar(pos_bars + pos, this_lw_list, width, label=d)

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
            axs[1].set_xticklabels(job_labels, rotation=45)
            axs[1].set_ylim(np.add(lw_ylim, [-3, +3]))
            axs[1].grid('on')

            fig.subplots_adjust(bottom=0.2)
            fig.show()

    def get_te_list(self):
        """
        Return the TEs for all the data signals in this pipeline.

        Returns
        -------
        [s.te for s in self._data_list] : list
            TEs for all signals stored in here
        """
        return([s.te for s in self._data_list])

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
        # --- init stuff ---
        # init some private attributes
        self._data_list = []
        self._data_ref_filepaths_list = []
        self._data_ref_list = []
        self._data_physio_list = []
        self._display_legends_list = []
        self._analyze_results_dict = {}

        # --- parsing some attributes: filepaths, legends, etc. ---
        log.info_line________________________()
        log.info("checking some stuff...")
        log.info_line________________________()

        # data files: parse
        log.info("parsing data file paths...")
        if(self.data_filepaths == []):
            parsed_list = None
        elif(type(self.data_filepaths) == list):
            parsed_list = self.data_filepaths
            parsed_list = [e.strip() for e in parsed_list]
            parsed_list = [None if e == "" else e for e in parsed_list]
        else:
            parsed_list = _parse_string_into_list(self.data_filepaths)

        if(parsed_list is not None):
            self._data_filepaths_list = parsed_list
        else:
            log.error("no data scan files specified!")

        # ref data files: parse
        log.info("parsing ref. data file paths...")

        if(self.data_ref_filepaths == []):
            parsed_list = None
        elif(type(self.data_ref_filepaths) == list):
            parsed_list = self.data_ref_filepaths
            parsed_list = [e.strip() for e in parsed_list]
            parsed_list = [None if e == "" else e for e in parsed_list]
        else:
            parsed_list = _parse_string_into_list(self.data_ref_filepaths)

        if(parsed_list is not None):
            self._data_ref_filepaths_list = parsed_list
        else:
            self._data_ref_filepaths_list = [None] * len(self._data_filepaths_list)

        # physio files: parse
        log.info("parsing physio data file paths...")
        if(self.data_physio_filepaths == []):
            parsed_list = None
        elif(type(self.data_physio_filepaths) == list):
            parsed_list = self.data_physio_filepaths
            parsed_list = [e.strip() for e in parsed_list]
            parsed_list = [None if e == "" else e for e in parsed_list]
        else:
            parsed_list = _parse_string_into_list(self.data_physio_filepaths)

        if(parsed_list is not None):
            self._data_physio_list = parsed_list
        else:
            self._data_physio_list = [None] * len(self._data_filepaths_list)

        # legend captions: parse
        log.info("parsing legend captions...")
        if(self.display_legends == []):
            parsed_list = None
        elif(type(self.display_legends) == list):
            parsed_list = self.display_legends
            parsed_list = [e.strip() for e in parsed_list]
            parsed_list = [None if e == "" else e for e in parsed_list]
        else:
            parsed_list = _parse_string_into_list(self.display_legends)

        if(parsed_list is not None):
            self._display_legends_list = parsed_list
        else:
            self._display_legends_list = [""] * len(self._data_filepaths_list)

        # --- checking consistency of some pipeline attributes ---
        # before reading data and stuff, let's check the list dimensions are consistent
        log.info("checking parameter list sizes consistency...")
        if(len(self._data_ref_filepaths_list) == 0):
            log.warning("no data ref. scans specified, that's ok, it is optional...")
        if(len(self._data_physio_list) == 0):
            log.warning("no physio log files specified, that's ok, it is optional...")

        n = len(self._data_filepaths_list)
        if(len(self._data_ref_filepaths_list) != n and len(self._data_ref_filepaths_list) > 0):
            log.error("weird, not the same number of data files (" + str(n) + ") and data ref. scans (" + str(len(self._data_ref_filepaths_list)) + ")?!")
        if(len(self._data_physio_list) != n and len(self._data_physio_list) > 0):
            log.error("weird, not the same number of data files (" + str(n) + ") and physio log files (" + str(len(self._data_physio_list)) + ")?!")
        if(len(self._display_legends_list) != n):
            log.error("weird, not the same number of data files (" + str(n) + ") and legend captions (" + str(len(self._display_legends_list)) + ")?!")

        # check if there are some blanks
        if(any(d is None for d in self._data_filepaths_list)):
            log.error("hey, you missed a line in the data scan file path list?!")
        if(any(d is None for d in self._display_legends_list)):
            log.error("hey, you missed a line in the legend caption list?!")

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
        log.info_line________________________()

        # --- reading data ---
        # now let's read the data files
        log.info_line_break()
        log.info_line________________________()
        log.info("reading data files...")
        log.info_line________________________()
        max_len_patient_name = 0
        for data_fn, physio_fn, data_leg in zip(self._data_filepaths_list, self._data_physio_list, self._display_legends_list):
            log.info_line_break()
            log.info_line________________________()
            log.info("reading data [" + data_leg + "]")
            log.info_line________________________()
            s = MRSData2(data_fn, self.data_coil_nChannels, physio_fn)
            log.debug("got a " + str(s.shape) + " vector")
            # store
            self._data_list.append(s)
            # patient name length
            if(s.patient_name is not None and len(s.patient_name) > max_len_patient_name):
                max_len_patient_name = len(s.patient_name)
            log.info_line________________________()

        log.info_line_break()
        log.info_line________________________()
        log.info("reading ref. data files...")
        log.info_line________________________()
        for data_ref_fn, data_leg in zip(self._data_ref_filepaths_list, self._display_legends_list):
            if(data_ref_fn is not None):
                log.info_line_break()
                log.info_line________________________()
                log.info("reading ref. data [" + data_leg + "]")
                log.info_line________________________()
                s = MRSData2(data_ref_fn, self.data_coil_nChannels, None)
                log.info("got a " + str(s.shape) + " vector")
                # if several averages, mean now (that could be a problem!?)
                s = s.mean(axis=0)
                # add 1 dimension
                s = s.reshape((1,) + s.shape)
                log.debug("reshaped to a " + str(s.shape) + " vector")
            else:
                s = None
            # store
            self._data_ref_list.append(s)
            log.info_line________________________()

        # --- setting up legends and offsets ---
        # legends: add patient name and format (space pads)
        log.debug_line_break()
        log.debug_line________________________()
        log.debug("processing displays legends and offsets...")
        display_legends_list_indexed = []
        for i, l in enumerate(self._display_legends_list):
            # add an index
            this_new_legend = ("#%d " % i) + l
            # set new legend
            display_legends_list_indexed.append(this_new_legend)
            self._data_list[i].set_display_label(this_new_legend)
            # new display label for ref. data reference, if any
            if(self._data_ref_list[i] is not None):
                self._data_ref_list[i].set_display_label(this_new_legend + " [REF]")

            # y offset
            self._data_list[i].set_display_offset(self.display_offset * i)

        self._display_legends_list = display_legends_list_indexed

        log.debug_line________________________()
        log.debug("referencing all scans to %.2fppm..." % self.ppm0)
        # set ppm0 for all signals
        for i, d in enumerate(self._data_list):
            self._data_list[i].ppm0 = self.ppm0
        for i, d in enumerate(self._data_ref_list):
            if d is not None:
                self._data_ref_list[i].ppm0 = self.ppm0

        log.debug_line________________________()
        log.debug("linking each dataset to its reference...")
        # --- linking each data set to its reference
        for i, d in enumerate(self._data_ref_list):
            # the ref scan for this dataset
            self._data_list[i].data_ref = d

        # --- reading job list ---
        log.info_line________________________()
        log.info("reading your job list...")
        log.info_line________________________()
        # first, for each job reset the display_range_ppm parameter if any
        # BECAUSE PYTHON CANNOT DO POINTERS -_-
        for j in list(self.jobs.keys()):
            if("display_range_ppm" in self.jobs[j]):
                self.jobs[j]["display_range_ppm"] = self.jobs["displaying"]["range_ppm"]

        # for each job
        for k, j in enumerate(self.job_list):
            this_job_name = j[0]["name"]
            log.info("#%d %s" % (k, this_job_name))
        log.info_line________________________()

        # --- running job list now ---
        # exception here: if any concatenate, we need to run the job list in two parts (before and after concatenation in order to reload the processed data)
        log.info_line_break()
        log.info_line________________________()
        log.info("running job list...")
        log.info_line________________________()
        # init job stacks
        concatenate_loop = int(self.jobs["concatenate"] in self.job_list) + 1
        jobs_stack_init = self.job_list[::-1].copy()
        jobs_stack = jobs_stack_init.copy()
        jobs_done_stack_init = []
        jobs_done_stack = jobs_done_stack_init.copy()
        for k in range(concatenate_loop):
            # resuming job stack after concatenate
            jobs_stack_init = jobs_stack.copy()
            jobs_done_stack_init = jobs_done_stack.copy()
            # for each dataset, the list of data processing functions with the right arguments
            for k, data in enumerate(self._data_list):
                log.info_line_break()
                log.info_line________________________()
                log.info("processing [" + data.display_label + "]")
                log.info_line________________________()

                # run job list for this dataset
                jobs_stack = jobs_stack_init.copy()
                jobs_done_stack = jobs_done_stack_init.copy()
                while(len(jobs_stack) > 0):
                    # job pop
                    job = jobs_stack.pop()
                    # if concatenate, get out of this dataset job loop
                    if(job == self.jobs["concatenate"]):
                        break

                    # run job on this dataset
                    job_result = self._run_job(job, data)

                    # push job in the stack of jobs done
                    jobs_done_stack.append(job)

                    # replace with processed signal
                    if(type(job_result) == MRSData2):
                        data = job_result
                        # measure SNR/LW after this process?
                        if(self.analyze_enable):
                            self._analyze(job_result, job, jobs_done_stack)
                            log.info_line________________________()

                # we finish running all jobs on a dataset, storing
                self._data_list[k] = data

            # if last job was concatenate, that means all datasets were processed and are ready for it
            if(job == self.jobs["concatenate"]):
                job_name = self.jobs["concatenate"][0]["name"]
                log.info("%s..." % job_name)
                s_concatenated = self._data_list[0]
                for i in range(1, len(self._data_list)):
                    s_concatenated = s_concatenated.concatenate_2d(self._data_list[i])

                # empty and store the single concatenated signal in the data list
                s_concatenated.set_display_label(s_concatenated.display_label + " (concatenated)")
                self._data_list = []
                self._data_list.append(s_concatenated)

                # if analyze going on, empty the results dict to avoid conflicts
                self._analyze_results_dict = {}

        # --- summary final linewidths ---
        if(self.analyze_enable):
            self._display_analyze_results()

        log.info("pipeline terminated!")
        log.info("returning processed signals!")
        return(self._data_list)

    def display_results(self):
        """Plot final processed datasets."""
        log.info("displaying final processed data...")

        # get display job
        disp_job = self.jobs["displaying"]

        # run a diplay job for each dataset
        for d in self._data_list:
            # run job on this dataset
            self._run_job(disp_job, d)

    def save(self, rdb):
        """
        Save each data set and the current pipeline to the PKL data storage file.

        Parameters
        ----------
        rdb: data_db object
            Data storage object
        """
        log.info("saving %d dataset(s) to file [%s]..." % (len(self._data_list), rdb.db_file))

        # for each dataset
        for d in self._data_list:
            rdb.save(d, self)


class voi_pipeline:
    """The voi_pipeline class is similar to pipeline class above but for VOI trace signals."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("You are trying to dynamically create a new attribute (%s) to this object and that is not cool! I will not let you do that because I believe it is a bad habit and can lead to terrible bugs. A clean way of doing this is to initialize your attribute (%s) in the __init__ method of this class. Bisou, bye :)" % (key, key))
        object.__setattr__(self, key, value)

    def __init__(self):
        # full paths to data files
        self.data_filepaths = ""

        # private attributes
        self._data_filepaths_list = []
        self._data_list = []
        self._xdata = []
        self._ydata = []
        self._zdata = []
        self._xdata_axis = []
        self._ydata_axis = []
        self._zdata_axis = []

        # display options
        self.display_fig_index = 1
        self.display_legends = ""
        self._display_legends_list = []

        # regions in the spatial profile to analyze (number of points)
        self.analyze_selectivity_range_list = [[800, 3550], [-10600, -7800], [-3650, 1850]]
        self.analyze_selectivity_list = np.array([])

        # freeze objet
        self.__isfrozen = True

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

    def run(self):
        """
        Run the standard pipeline for VOI trace signals. It includes.

            * reading the data
            * displaying the spatial selection profiles
            * analyzing the selectivity
        """
        # parse
        log.info_line________________________()
        log.info("parsing data file paths...")
        log.info_line________________________()
        if(self.data_filepaths == []):
            parsed_list = None
        elif(type(self.data_filepaths) == list):
            parsed_list = self.data_filepaths
        else:
            parsed_list = _parse_string_into_list(self.data_filepaths)

        if(parsed_list is None):
            log.error("no data scan files specified!")
        else:
            self._data_filepaths_list = parsed_list
        log.info_line________________________()

        # legends
        log.info_line________________________()
        log.info("parsing legend captions...")
        log.info_line________________________()
        if(self.display_legends == []):
            parsed_list = None
        elif(type(self.display_legends) == list):
            parsed_list = self.display_legends
        else:
            parsed_list = _parse_string_into_list(self.display_legends)

        if(parsed_list is None):
            # oh no, no legend captions were specified
            # let's create them using filenames
            self._display_legends_list = _build_legend_list_from_filepath_list(self._data_filepaths_list)
        else:
            self._display_legends_list = parsed_list
        log.info_line________________________()

        # read the 3 VOI traces per dataset
        log.info_line________________________()
        log.info("reading data files...")
        log.info_line________________________()
        for f in self._data_filepaths_list:
            # read data
            log.info("looking in folder: ")
            log.info(f)
            log.info("reading the 3 dicom files...")
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
        log.info_line________________________()

        # analysis
        log.info_line________________________()
        log.info("evaluating selectivity using a " + str(self.analyze_selectivity_range_list) + "ppm ranges...")
        log.info_line________________________()
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

            log.info("selectivity results for [" + leg + "]...")
            log.info(" [X] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sx_in, sx_out))
            log.info(" [Y] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sy_in, sy_out))
            log.info(" [Z] IN=%.0fHz-1 OUT=%.0f Hz-1" % (sz_in, sz_out))

            # store
            self.analyze_selectivity_list[k, :, 0] = [sx_in, sy_in, sz_in]
            self.analyze_selectivity_list[k, :, 1] = [sx_out, sy_out, sz_out]
            k = k + 1

            '''
            log.info("selectivity debug info for [" + leg + "] coming...")
            log.info(" [X] Intersection at %f Hz and %f Hz" % (sx_ax[sx_ithreshold_l],sx_ax[sx_ithreshold_r]))
            log.info(" [Y] Intersection at %f Hz and %f Hz" % (sy_ax[sy_ithreshold_l],sy_ax[sy_ithreshold_r]))
            log.info(" [Z] Intersection at %f Hz and %f Hz" % (sz_ax[sz_ithreshold_l],sz_ax[sz_ithreshold_r]))
            '''
        log.info_line________________________()

        # display
        log.info_line________________________()
        log.info("displaying the 3 spatial profiles...")
        log.info_line________________________()
        fig = plt.figure(self.display_fig_index)
        fig.clf()
        axs = fig.subplots(2, 3)
        fig.canvas.set_window_title("mrs.reco.voi_pipeline.run")
        fig.suptitle("spatial selection profiles")
        for (sx, sy, sz, sx_ax, sy_ax, sz_ax, leg) in zip(self._xdata, self._ydata, self._zdata, self._xdata_axis, self._ydata_axis, self._zdata_axis, self._display_legends_list):

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
            axs[0, i].axvline(x=self.analyze_selectivity_range_list[i][0], linestyle='--')
            axs[0, i].axvline(x=self.analyze_selectivity_range_list[i][1], linestyle='--')
            axs[1, i].axvline(x=self.analyze_selectivity_range_list[i][0], linestyle='--')
            axs[1, i].axvline(x=self.analyze_selectivity_range_list[i][1], linestyle='--')

        axs[1, 2].legend()
        fig.subplots_adjust()
        fig.show()


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
