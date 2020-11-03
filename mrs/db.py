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
import scipy.optimize as optimize
import pandas as pd
import os
import pickle
from datetime import datetime, timedelta
from shutil import copyfile
import pandas as pd
import copy

import pdb


class data_db():
    """A class used to deal with storage of reconstructed signals in pkl files.
    The database consists of a dataframe."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, data_db_file=default_paths.DEFAULT_RECO_DATA_DB_FILE):
        """
        Initialize the reconstructed data storage, basically creates an empty PKL file if nothing already exists.

        Parameters
        ----------
        data_db_file: string
            PKL file where all the data is stored
        """
        # init df
        self._df = pd.DataFrame(columns = ['hash' ,
                                         'patient',
                                         'study' ,
                                         'dataset',
                                         'reco_pipeline',
                                         'fit_results',
                                         'timestamp'])

        # hash of scan will be our index (should be unique)
        self._df = self._df.set_index('hash')

        self.db_file = data_db_file

        # if file does not exist, create it
        if(not os.path.isfile(self.db_file)):
            log.info("creating storage file [%s]..." % self.db_file)
            self._save_db_file()

        # first read
        self._read_db_file()

    @property
    def df(self):
        """
        Property get function for df.

        Returns
        -------
        self._df : dataframe
            Database
        """
        return(self._df)

    def _read_db_file(self):
        """Read the content of PKL file and store it."""
        log.debug("reading db file [%s]..." % self.db_file)

        # now open pkl file
        with open(self.db_file, 'rb') as f:
            [pkl_db_file, pkl_df] = pickle.load(f)

        self.db_file = pkl_db_file
        self._df = pkl_df
        log.debug("finished reading!")

    def _save_db_file(self):
        """Save the content of the dict to the PKL file."""
        log.debug("saving to db file [%s]..." % self.db_file)

        # write back pkl file
        with open(self.db_file, 'wb') as f:
            pickle.dump([self.db_file, self._df], f)

        log.debug("finished saving!")

    def save_reco_dataset(self, d, rp=None):
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

        # reload db file
        self._read_db_file()

        # if we reached here, that means the PKL file is not corrupted
        # let's make a backup of it
        copyfile(self.db_file, self.db_file + ".bak")

        # get timestamp
        ts = datetime.now()

        # get hash of data
        if(d["raw"]["data"] is None):
            h = d["dcm"]["data"].data_file_hash
        else:
            h = d["raw"]["data"].data_file_hash

        # get patient/study ids
        patient_id, study_id = self._extract_patient_study_num(d)

        # look for hash in df
        if(h in self._df.index):
            log.debug("scan [#%s] already exists!" % h)
            log.debug("--> updating entry in db!")
        else:
            log.debug("creating new entry for scan [#%s]!" % h)

        # add entry in db
        self._df.loc[h] = [patient_id, study_id, d, rp, None, ts]

        # save db file
        self._save_db_file()

    def save_fit(self, data_hash, fit_key, fit_results):
        """
        Save mrs.reco.pipeline dataset and its reco pipeline and deal with conflicts.

        Parameters
        ----------
        data_hash: string
            Hash of fitted data
        fit_key: string
            String used to store the fit results
        fit_results: dict
            Fit results stored as a dict
        """
        log.info("saving fit to file [%s]..." % self.db_file)

        # reload db file
        self._read_db_file()

        # if we reached here, that means the PKL file is not corrupted
        # let's make a backup of it
        copyfile(self.db_file, self.db_file + ".bak")

        # find data hash in df
        if(data_hash not in self._df.index):
            log.error("you are trying to store some fit results for a scan not present in the db! :(")
        else:
            if(self.df.loc[data_hash]["fit_results"] is None):
                # never stored fit results before?
                log.debug("creating a new [fit_results] entry for scan [%s]!" % data_hash)
                self.df.at[data_hash, "fit_results"] = {fit_key: fit_results}
            else:
                # update
                log.debug("updating [fit_results] for scan [%s]!" % data_hash)
                self.df.at[data_hash, "fit_results"][fit_key] = fit_results

        # save db file
        self._save_db_file()

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

            # add the resulting parameter names and the original object
            par_name_list = par_name_list + ["obj"] + name_list
            par_val_list = par_val_list + [var] + val_list

        # if scraping an optim result
        elif(isinstance(var, optimize.OptimizeResult)):
            # this is actually a dict but the __dict__ is empty for some reason (!?)
            name_list = list(var.keys())
            val_list = list(var.values())

            # add the resulting parameter names and the original object
            par_name_list = par_name_list + ["obj"] + name_list
            par_val_list = par_val_list + [var] + val_list

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

    def create_big_df(self):
        """
        Extract info from stored objects and dicts and create an extended dataframe.

        Returns
        -------
        df : dataframe
            Resulted df object.
        """
        log.info("extending dataframe by scraping data...")

        # reload db file
        self._read_db_file()

        # init a list of df
        df_list = []

        # browse db
        for i, row in self.df.iterrows():
            this_column_list, this_values_list = self._scrap_data(row.to_dict())
            # add hash
            this_column_list = ["hash"] + this_column_list
            this_values_list = [i] + this_values_list
            # convert to df (that makes one line)
            df = pd.DataFrame([this_values_list], columns=this_column_list)
            # store
            df_list.append(df)

        # creating final df by merging all lines
        df = pd.DataFrame()
        df = df.append(df_list)
        # hash of scan will be our index (should be unique)
        df = df.set_index('hash')

        return(df)
