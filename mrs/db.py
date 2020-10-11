#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stuff related to database handling. These classes help to store reconstructed data for later group analysis.

@author: tangir
"""

from mrs import paths as default_paths
from mrs import log
import os
import pickle
from datetime import datetime, timedelta
from shutil import copyfile
import numpy as np

import pdb


def _check_func_test(dataset_entry):
    """
    Return true if you want to select this dataset when getting data from database using the get_datasets or print methods. This is only a dummy example. Please feel free to recode a function with the same prototype to select for example only 3T scans, or short-TE scans, etc.

    Parameters
    ----------
    dataset_entry: dict
        Dict entry in db

    Returns
    -------
    r : boolean
        True if we keep this dataset/pipeline
    """
    r = (dataset_entry["study"].isdigit() and int(dataset_entry["study"]) < 5 and
         dataset_entry["dataset"]["raw"]["data"].te > 0.0 and
         dataset_entry["dataset"]["raw"]["data"].na > 0)  # that is stupid, just an example
    return(r)


class data_db(dict):
    """A class used to deal with storage of reconstructed signals with their respective reco and fit pipelines in pkl files. The database consist of a dict, keys are patient IDs, subkeys are study IDs. A study contains a dataset (see mrs.reco.pipeline), reco and fit pipelines and a timestamp."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, reco_data_db_file=default_paths.DEFAULT_RECO_DATA_DB_FILE):
        """
        Initialize the reconstructed data storage, basically creates an empty PKL file if nothing already exists.

        Parameters
        ----------
        reco_data_db_file: string
            PKL file where all the data is stored
        """
        # init dict
        super(data_db, self).__init__()

        # attach dict to pkl file
        if(reco_data_db_file is None):
            self.db_file = None
        else:
            self.db_file = reco_data_db_file

            # if file does not exist, create it
            if(not os.path.isfile(self.db_file)):
                log.info("creating storage file [%s]..." % self.db_file)
                # write pkl file
                with open(self.db_file, 'wb') as f:
                    pickle.dump([{}], f)

            # first read
            self._read_db_file()

    def _is_linked_to_file(self):
        """
        Return true if this db is mirrored to a file on the disk.

        Returns
        -------
        r : bool
            True if this data_db object is stored on the disk
        """
        r = (self.db_file is not None)
        return(r)

    def _read_db_file(self):
        """Read the content of PKL file and store it."""
        log.debug("reading db file [%s]..." % self.db_file)

        # now open pkl file
        with open(self.db_file, 'rb') as f:
            [pkl_data_dict] = pickle.load(f)

        for d in list(pkl_data_dict.keys()):
            self[d] = pkl_data_dict[d]

    def _save_db_file(self):
        """Save the content of the dict to the PKL file."""
        log.debug("saving to db file [%s]..." % self.db_file)

        # write back pkl file
        with open(self.db_file, 'wb') as f:
            pickle.dump([self], f)

    def select_datasets(self, check_func=_check_func_test):
        """
        Return datasets which passes the check function.

        Parameters
        ----------
        check_func : function
            Function with one argument (database dict subentry) that you should code and which returns True if you want to extract this entry out of the database.

        Returns
        -------
        selected_data_db : data_db
            New data_db object after filtering. This data_db object will not be attached to any file.
        """
        log.info("getting datasets using check function [%s]..." % check_func)

        # init
        selected_data_db = data_db(None)

        # reload db file if needed
        if(self._is_linked_to_file()):
            self._read_db_file()

        # for each dataset
        for patient_key in self:
            for study_key in self[patient_key]:
                study_entry = self[patient_key][study_key]
                if(check_func(study_entry)):
                    # passed the test --> adding this study to the new db

                    # init patient if needed
                    if(patient_key not in list(selected_data_db.keys())):
                       selected_data_db[patient_key] = {}

                    # add study
                    selected_data_db[patient_key][study_entry] = study_entry

        # print extracted datasets
        selected_data_db.print()

        return(selected_data_db)

    def get_latest_dataset(self):
        """
        Return the most recent study saved.

        Returns
        -------
        selected_data_db : data_db
            New data_db object after filtering. This data_db object will not be attached to any file.
        """
        log.info("looking for latest dataset in db file [%s]..." % self.db_file)

        # init
        selected_data_db = data_db(None)
        found_patient_key = None
        found_study_key = None

        # reload db file if needed
        if(self._is_linked_to_file()):
            self._read_db_file()

        # for each dataset
        ts_diff = timedelta(days=99999999)
        for patient_key in self:
            for study_key in self[patient_key]:
                study_entry = self[patient_key][study_key]
                this_ts = study_entry["timestamp"]
                if(ts_diff > (datetime.now() - this_ts)):
                    ts_diff = (datetime.now() - this_ts)
                    found_patient_key = patient_key
                    found_study_key = study_key

        # add found study
        selected_data_db[found_patient_key] = {}
        selected_data_db[found_patient_key][found_study_key] = self[found_patient_key][found_study_key]

        return(selected_data_db)

    def print(self):
        """
        Print a summary of the database content.
        """
        # init
        patient_key_max_len = 0

        # reload db file if needed
        if(self._is_linked_to_file()):
            self._read_db_file()

        # first scan to get longest patient name
        for patient_key in self:
            if(len(str(patient_key)) > patient_key_max_len):
                patient_key_max_len = len(str(patient_key))

        if(patient_key_max_len < 10):
            patient_key_max_len = 10

        # now browse and print
        print("[Patient]".ljust(patient_key_max_len) + "[Study]".ljust(10) + "[RAW]".ljust(10) + "[DCM]".ljust(10) + "[RECO]".ljust(10) + "[FIT]".ljust(10))
        for patient_key in self:
            for study_key in self[patient_key]:
                study_entry = self[patient_key][study_key]
                print(str(patient_key).ljust(patient_key_max_len) +
                      str(study_key).ljust(10) +
                      str(study_entry["dataset"]["raw"] is not None).replace("True", " x ").replace("False", "").ljust(10) +
                      str(study_entry["dataset"]["dcm"] is not None).replace("True", " x ").replace("False", "").ljust(10) +
                      str(study_entry["reco_pipeline"] is not None).replace("True", " x ").replace("False", "").ljust(10) +
                      str(study_entry["fit_pipeline"] is not None).replace("True", " x ").replace("False", "").ljust(10))

    def _extract_patient_study_num(self, d):
        """
        Extract the "numero d'inclusion" and "numero de passage" from patient name, when possible.

        Parameters
        ----------
        d: mrs.reco.pipeline dataset entry (dict)
            Dict containing data and stuff

        Returns
        -------
        patient_num : int
            Found dataset
        study_num : int
            Corresponding pipeline
        """
        # find the first non-digit character
        if(d["raw"]["data"] is not None):
            data_raw_patient_name = d["raw"]["data"].patient["name"]
        elif(d["dcm"]["data"] is not None):
            data_raw_patient_name = d["dcm"]["data"].patient["name"]
        else:
            log.error("cannot find patient info. in the raw or dcm file headers!")

        log.debug("extracting patient and study IDs from [%s]..." % data_raw_patient_name)
        data_raw_patient_name = data_raw_patient_name.replace("-", "_")
        data_raw_patient_name = data_raw_patient_name.replace(".", "_")
        data_raw_patient_name_splitted = data_raw_patient_name.split("_")

        patient_num = data_raw_patient_name_splitted[0]
        if(patient_num.isdigit()):
            patient_num = int(patient_num)
        else:
            patient_num = None

        for s in data_raw_patient_name_splitted:
            study_num = s
            if(study_num[0] == 'P' and study_num[1:].isdigit()):
                study_num = int(study_num[1:])
                break
            else:
                study_num = None

        return(patient_num, study_num)

    def save_dataset(self, d, rp=None):
        """
        Save mrs.reco.pipeline dataset and its reco pipeline and deal with conflicts.

        Parameters
        ----------
        d: mrs.reco.pipeline dataset entry (dict)
            Dict containing data and stuff
        rp: pipeline
            Reco pipeline used to get this data
        """
        log.debug("saving dataset to file [%s]..." % self.db_file)

        # reload db file if needed
        if(self._is_linked_to_file()):
            self._read_db_file()

        # if we reached here, that means the PKL file is not corrupted
        # let's make a backup of it
        copyfile(self.db_file, self.db_file + ".bak")

        # try to extract patient/study nums if any
        patient_id, study_id = self._extract_patient_study_num(d)

        if(patient_id is None or study_id is None):
            log.warning("could not extract any patient/study IDs!")
            patient_id = d["raw"]["data"].patient["name"]
            study_id = 0
            log.warning("using patient name instead: [#%d / P%d] !" % (patient_id, study_id))

        # check if dataset already stored, based in patient/study IDs
        if(patient_id in list(self.keys())):
            if(study_id in list(self[patient_id].keys())):
                log.debug("dataset [#%d / P%d] already exists!" % (patient_id, study_id))
                log.debug("--> updating dataset!")
            else:
                log.debug("dataset patient [#%d] already exists!" % patient_id)
                log.debug("--> adding dataset as a new study [P%d]!" % study_id)
        else:
            log.debug("--> creating new patient [#%d] already exists!" % patient_id)
            log.debug("--> adding dataset as a new study [P%d]!" % study_id)
            # creating key
            self[patient_id] = {}

        # add/update with the dataset
        ts = datetime.now()
        self[patient_id][study_id] = {"patient": patient_id,
                                      "study": study_id,
                                      "dataset": d,
                                      "reco_pipeline": rp,
                                      "fit_pipeline": None,
                                      "timestamp": ts}

        # save db file if needed
        if(self._is_linked_to_file()):
            self._save_db_file()
