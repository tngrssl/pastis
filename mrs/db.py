#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stuff related to database handling. These classes help to store reconstructed data for later group analysis.

@author: tangir
"""

from mrs import paths as default_paths
from mrs import reco
from mrs import sim
from mrs import log
import os
import pickle
from datetime import datetime, timedelta
from shutil import copyfile
import numpy as np
import pandas as pd

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
            # create file on disk etc
            self.link_to_file(reco_data_db_file)

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

    def link_to_file(self, data_db_file):
        """
        Initialize link to file on disk.

        Parameters
        ----------
        data_db_file: string
            PKL file where all the data is stored
        """
        self.db_file = data_db_file

        # if file does not exist, create it
        if(not os.path.isfile(self.db_file)):
            log.info("creating storage file [%s]..." % self.db_file)
            # write pkl file
            with open(self.db_file, 'wb') as f:
                pickle.dump([{}], f)

        # first read
        self._read_db_file()

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
                for hash_key in self[patient_key][study_key]:
                    this_entry = self[patient_key][study_key][hash_key]
                    if(check_func(this_entry)):
                        # passed the test --> adding this study to the new db

                        # init patient if needed
                        if(patient_key not in list(selected_data_db.keys())):
                           selected_data_db[patient_key] = {}
                        if(study_key not in list(selected_data_db[patient_key].keys())):
                           selected_data_db[patient_key][study_key] = {}

                        # add study
                        selected_data_db[patient_key][study_key][hash_key] = this_entry

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
        found_hash_key = None

        # reload db file if needed
        if(self._is_linked_to_file()):
            self._read_db_file()

        # for each dataset
        ts_diff = timedelta(days=99999999)
        for patient_key in self:
            for study_key in self[patient_key]:
                for hash_key in self[patient_key][study_key]:
                    this_entry = self[patient_key][study_key][hash_key]
                    this_ts = this_entry["timestamp"]
                    if(ts_diff > (datetime.now() - this_ts)):
                        ts_diff = (datetime.now() - this_ts)
                        found_patient_key = patient_key
                        found_study_key = study_key
                        found_hash_key = hash_key

        # add found study
        selected_data_db[found_patient_key] = {}
        selected_data_db[found_patient_key][found_study_key] = {}
        selected_data_db[found_patient_key][found_study_key][found_hash_key] = self[found_patient_key][found_study_key][found_hash_key]

        return(selected_data_db)

    def _scrap_data(self, var, prefix_str=None):
        """
        Scrap data in db reading dicts, lists and some objects' attributes and converting it to lists.

        Parameters
        ----------
        var : ?
            Could be anything but likely a dict, list, a mrs.* object or int/float/etc.

        Returns
        -------
        par_name_list : list
            Name of parameters scraped.
        par_val_list : list
            Values of parameters scraped.
        """
        # init lists
        par_name_list = []
        par_val_list = []

        # if scraping a dict
        if(type(var) is dict):
            # browse
            for this_dict_key in var:
                # recursive call
                name_list, val_list = self._scrap_data(var[this_dict_key], str(this_dict_key))

                # add the resulting parameter names
                par_name_list = par_name_list + name_list
                par_val_list = par_val_list + val_list

        # if scraping a list
        elif(type(var) is list):
            # browse
            for ind, this_item in enumerate(var):
                name_list, val_list = self._scrap_data(this_item, "[" + str(ind) + "]")

                # add the resulting parameter names
                par_name_list = par_name_list + name_list
                par_val_list = par_val_list + val_list

        # if scraping an object instance of a class we wrote in the package
        elif(isinstance(var, (reco.MRSData2, reco.pipeline, sim.mrs_sequence))):
            # scrap the dict attribute
            name_list, val_list = self._scrap_data(var.__dict__)

            # add the resulting parameter names
            par_name_list = par_name_list + name_list
            par_val_list = par_val_list + val_list

        # if scraping something else, probably some int or stuff
        else:
            # make it simple
            par_name_list = par_name_list + [""]
            par_val_list = par_val_list + [var]

        # format the field name, usefull later in pandas df
        if(prefix_str is not None):
            if(len(par_name_list) > 1 and type(var) is not list):
                prefix_str += "_"

            par_name_list = [prefix_str + s for s in par_name_list]

        return(par_name_list, par_val_list)

    def to_dataframe(self):
        """
        Convert this huge db structure into a pandas dataframe object.

        Returns
        -------
        par_name_list : dataframe
            Resulted df object.
        """
        log.info("converting db into pandas dataframe...")

        # init a list of df
        df_list = []

        # browse db
        for patient_key in self:
            for study_key in self[patient_key]:
                for hash_key in self[patient_key][study_key]:
                    # scrap this scan entry
                    this_entry = self[patient_key][study_key][hash_key]
                    this_column_list, this_values_list = self._scrap_data(this_entry)
                    # convert to df (that makes one line)
                    df = pd.DataFrame([this_values_list], columns=this_column_list)
                    # store
                    df_list.append(df)

        # creating final df by merging all lines
        df = pd.DataFrame()
        df = df.append(df_list)
        # hash of scan will be our index (should be unique)
        df.set_index('hash')

        return(df)

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
        print("[Patient]".ljust(patient_key_max_len) + "[Study]".ljust(10) + "[Scan (hash)]".ljust(34) + "[RAW]".ljust(6) + "[DCM]".ljust(6) + "[RECO]".ljust(7) + "[FIT]".ljust(6))
        for patient_key in self:
            for study_key in self[patient_key]:
                for hash_key in self[patient_key][study_key]:
                    this_entry = self[patient_key][study_key][hash_key]
                    print(str(patient_key).ljust(patient_key_max_len) +
                          str(study_key).ljust(10) +
                          str(hash_key).ljust(34) +
                          str(this_entry["dataset"]["raw"] is not None).replace("True", " x ").replace("False", "").ljust(6) +
                          str(this_entry["dataset"]["dcm"] is not None).replace("True", " x ").replace("False", "").ljust(6) +
                          str(this_entry["reco_pipeline"] is not None).replace("True", " x ").replace("False", "").ljust(7) +
                          str(this_entry["fit_pipeline"] is not None).replace("True", " x ").replace("False", "").ljust(6))

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

        # find hash
        if(d["raw"]["data"] is None):
            h = d["dcm"]["data"].data_file_hash
        else:
            h = d["raw"]["data"].data_file_hash

        # check if dataset already stored, based in patient/study IDs/hash
        if(patient_id in list(self.keys())):
            if(study_id in list(self[patient_id].keys())):
                log.debug("study [#%d / P%d] already exists!" % (patient_id, study_id))
                log.debug("--> updating study!")
                if(h in list(self[patient_id][study_id].keys())):
                    log.debug("dataset [#%d / P%d / %s] already exists!" % (patient_id, study_id, h))
                    log.debug("--> updating dataset!")
            else:
                log.debug("dataset patient [#%d] already exists!" % patient_id)
                log.debug("--> adding dataset as a new study [P%d]!" % study_id)
                self[patient_id][study_id] = {}
        else:
            log.debug("--> creating new patient [#%d] already exists!" % patient_id)
            log.debug("--> adding dataset as a new study [P%d]!" % study_id)
            # creating key
            self[patient_id] = {}
            self[patient_id][study_id] = {}

        # add/update with the dataset
        ts = datetime.now()

        self[patient_id][study_id][h] = {"patient": patient_id,
                                          "study": study_id,
                                          "hash": h,
                                          "dataset": d,
                                          "reco_pipeline": rp,
                                          "fit_pipeline": None,
                                          "timestamp": ts}

        # save db file if needed
        if(self._is_linked_to_file()):
            self._save_db_file()

    def save_fit(self, data_hash, fit_key, fit_results):
        """
        Save mrs.reco.pipeline dataset and its reco pipeline and deal with conflicts.

        Parameters
        ----------
        data_hash: mrs.reco.pipeline dataset entry (dict)
            Dict containing data and stuff
        fit_key: string
            String used to store the fit results
        fit_results: dict
            Fit results stored as a dict
        """
        log.debug("saving fit to file [%s]..." % self.db_file)

        # reload db file if needed
        if(self._is_linked_to_file()):
            self._read_db_file()

        # if we reached here, that means the PKL file is not corrupted
        # let's make a backup of it
        copyfile(self.db_file, self.db_file + ".bak")

        # find data hash in dict
        for patient_key in self:
            for study_key in self[patient_key]:
                for hash_key in self[patient_key][study_key]:
                    if(hash_key == data_hash):
                        if(self[patient_key][study_key][hash_key]["fit_pipeline"] is None):
                            self[patient_key][study_key][hash_key]["fit_pipeline"] = {}
                        # store fit results
                        self[patient_key][study_key][hash_key]["fit_pipeline"][fit_key] = fit_results

        # save db file if needed
        if(self._is_linked_to_file()):
            self._save_db_file()
