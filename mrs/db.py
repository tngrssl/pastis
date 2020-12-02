#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stuff related to database handling. These classes help to store reconstructed data for later group analysis.

@author: tangir
"""

from mrs import paths as default_paths
from mrs import reco
from mrs import sim
from mrs import fit
from mrs import log
import scipy.optimize as optimize
import pandas as pd
import os
import numbers
import numpy as np
from datetime import datetime
from shutil import copyfile

import pdb


class data_db():
    """A class used to deal with storage of reconstructed signals in pkl files. The database consists of dataframes.
    One stores the reconstructed signals, the other the fits."""

    # frozen stuff: a technique to prevent creating new attributes
    # (https://stackoverflow.com/questions/3603502/prevent-creating-new-attributes-outside-init)
    __isfrozen = False

    def __setattr__(self, key, value):
        """Overload of __setattr__ method to check that we are not creating a new attribute."""
        if self.__isfrozen and not hasattr(self, key):
            log.error_new_attribute(key)
        object.__setattr__(self, key, value)

    def __init__(self, data_db_reco_file=default_paths.DEFAULT_RECO_DATA_DB_FILE):
        """
        Initialize the reconstructed data storage, basically creates empty PKL files if nothing already exists.

        Parameters
        ----------
        data_db_reco_file: string
            PKL file where all the reconstructed data is stored
        """
        data_db_fit_file = data_db_reco_file[:-4] + "_fit.pkl"
        log.info("initializing db files [in progress]")
        log.info(data_db_reco_file)
        log.info(data_db_fit_file)

        # init df reco
        self._df_reco = pd.DataFrame(columns = ['scan_hash',
                                                 'patient',
                                                 'study',
                                                 'dataset',
                                                 'reco_pipeline'])

        # hash of scan will be our index (should be unique)
        self._df_reco = self._df_reco.set_index('scan_hash')

        # if file does not exist, create it
        self.data_db_reco_file = data_db_reco_file
        if(not os.path.isfile(self.data_db_reco_file)):
            log.debug("creating storage file [%s]..." % self.data_db_reco_file)
            self.write_pickle_files(True, False)
        else:
            self.read_pickle_files(True, False)

        # init df fit
        self._df_fit = pd.DataFrame(columns = ['fit_hash',
                                                'scan_hash',
                                                'fit_strategy',
                                                'fit_pipeline',
                                                'fit_results'])

        # hash of scan will be our index (should be unique)
        self._df_fit = self._df_fit.set_index('fit_hash')

        # if file does not exist, create it
        self.data_db_fit_file = data_db_fit_file
        if(not os.path.isfile(self.data_db_fit_file)):
            log.debug("creating storage file [%s]..." % self.data_db_fit_file)
            self.write_pickle_files(False, True)
        else:
            self.read_pickle_files(False, True)

        log.info("initializing db files [done]")

    @property
    def df_reco(self):
        """
        Property get function for df_reco.

        Returns
        -------
        self._df_reco : dataframe
            Database
        """
        return(self._df_reco)

    @property
    def df_fit(self):
        """
        Property get function for df_fit.

        Returns
        -------
        self._df_fit : dataframe
            Database
        """
        return(self._df_fit)

    def read_pickle_files(self, read_reco=True, read_fit=True):
        """Read dataframes from pickle files.

        Parameters
        ----------
        read_reco: boolean
            Read df_reco from pickle file
        read_fit: boolean
            Read df_fit from pickle file
        """
        if(not read_reco and not read_fit):
            return()

        if(read_reco):
            log.debug("reading storage file [%s]..." % self.data_db_reco_file)
            self._df_reco = pd.read_pickle(self.data_db_reco_file)

        if(read_fit):
            log.info("reading storage file [%s]..." % self.data_db_fit_file)
            self._df_fit = pd.read_pickle(self.data_db_fit_file)

        log.info("reading storage file [done]...")

    def write_pickle_files(self, write_reco=True, write_fit=True):
        """Write dataframes to pickle files.

        Parameters
        ----------
        write_reco: boolean
            Read df_reco from pickle file
        write_fit: boolean
            Read df_fit from pickle file
        """
        if(not write_reco and not write_fit):
            return()

        if(write_reco):
            log.debug("writing storage file [%s]..." % self.data_db_reco_file)
            self._df_reco.to_pickle(self.data_db_reco_file)

        if(write_fit):
            log.debug("writing storage file [%s]..." % self.data_db_fit_file)
            self._df_fit.to_pickle(self.data_db_fit_file)

        log.info("writing storage files [done]")

    def save_reco_dataset(self, d, rp=None, write2file_now=True):
        """
        Save mrs.reco.pipeline dataset and its reco pipeline and deal with conflicts.

        Parameters
        ----------
        d: mrs.reco.pipeline dataset entry (dict)
            Dict containing data and stuff
        rp: pipeline
            Reco pipeline used to get this data
        write2file_now: boolean
            Write resulting dataframe to pickle file now (takes time)
        """
        log.info("saving dataset to file [%s]..." % self.data_db_reco_file)

        # read storage file to refresh data in memory
        self.read_pickle_files(write2file_now, False)

        # if we reached here, that means the PKL file is not corrupted
        # let's make a backup of it
        copyfile(self.data_db_reco_file, self.data_db_reco_file + ".bak")

        # get hash of data
        if(d["raw"]["data"] is None):
            h = d["dcm"]["data"].data_file_hash
        else:
            h = d["raw"]["data"].data_file_hash

        # get patient/study ids
        patient_id, study_id = self._extract_patient_study_num(d)

        # look for hash in df
        if(h in self._df_reco.index):
            log.debug("scan [#%s] already exists!" % h)
            log.debug("--> updating entry in db!")
        else:
            log.debug("creating new entry for scan [#%s]!" % h)

        # add/update entry in db
        self._df_reco.loc[h] = [patient_id, study_id, d, rp]

        # write data in memory to storage file
        self.read_pickle_files(write2file_now, False)

    def save_fit_results(self, scan_hash, fr, ft=None, fs=None, write2file_now=True):
        """
        Save fit pipeline, results and deal with conflicts.

        Parameters
        ----------
        scan_hash: string
            Hash of fitted data
        fr: dict
            Fit results stored as a dict
        ft: mrs.fit.fit_tool
            Fit pipeline object used
        fs: mrs.fit.fit_strategy
            Fit strategy object used
        write2file_now: boolean
            Write resulting dataframe to pickle file now (takes time)
        """
        log.info("saving fit to file [%s]..." % self.data_db_fit_file)

        # read storage file to refresh data in memory
        self.read_pickle_files(False, write2file_now)

        # if we reached here, that means the PKL file is not corrupted
        # let's make a backup of it
        copyfile(self.data_db_fit_file, self.data_db_fit_file + ".bak")

        # first, check if the data we fitted in stored in the df_reco
        if(scan_hash not in self._df_reco.index):
            log.error("you are trying to store some fit results from a scan which is not present in the reco_db! Weird. :(")

        # get fit hash: mix data hash and fit stategy
        h = fs.hashit(scan_hash)

        # look for hash in df
        if(h in self._df_fit.index):
            log.debug("fit [#%s] already exists!" % h)
            log.debug("--> updating entry in db!")
        else:
            log.debug("creating new entry for fit [#%s]!" % h)

        # add/update entry in db
        self._df_fit.loc[h] = [scan_hash, fs, ft, fr]

        # write data in memory to storage file
        self.read_pickle_files(False, write2file_now)

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

    def extend_df_reco(self, replace_none_raw_with_dcm=False):
        """
        Extract info from stored objects and dicts and create an extended dataframe.

        Parameters
        ----------
        replace_none_raw_with_dcm : boolean
            Deal with the tricky exception of datasets that only have dcm data and no raw

        Returns
        -------
        edf : dataframe
            Resulted df object
        """
        log.info("extending dataframe by scraping data...")

        # fix non raw data?
        if(replace_none_raw_with_dcm):
            df2scrap = self.df_reco.copy()
            for i, row in self.df_reco.iterrows():
                if(row["dataset"]["raw"]["data"] is None):
                    df2scrap.loc[i, "dataset"]["raw"] = row["dataset"]["dcm"]
        else:
            df2scrap = self.df_reco

        edf = _extend_df(df2scrap)
        return(edf)

    def extend_df_fit(self):
        """
        Extract info from stored objects and dicts and create an extended dataframe.

        Returns
        -------
        edf : dataframe
            Resulted df object
        """
        log.info("extending dataframe by scraping data...")

        edf = _extend_df(self.df_fit)
        return(edf)


def _extend_df(df):
    """
    Scrap data in db reading dicts, lists and some objects' attributes and converting it to lists.

    Parameters
    ----------
    df : dataframe
        Panda dataframe to extend

    Returns
    -------
    edf : dataframe
        Resulting extended dataframe
    """
    # init a list of df
    df_list = []

    # find index col of df
    df_index_name = df.index.name

    # browse db
    for i, row in df.iterrows():
        this_column_list, this_values_list = _scrap_data(row.to_dict())
        # add hash
        this_column_list = [df_index_name] + this_column_list
        this_values_list = [i] + this_values_list
        # convert to df (that makes one line)
        this_df = pd.DataFrame([this_values_list], columns=this_column_list)
        # store
        df_list.append(this_df)

    # creating final df by merging all lines
    edf = pd.DataFrame()
    edf = edf.append(df_list)
    # hash of scan will be our index (should be unique)
    edf = edf.set_index(df_index_name)

    return(edf)


def _scrap_data(var, prefix_str=None):
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
            name_list, val_list = _scrap_data(var[this_dict_key], str(this_dict_key))

            # add the resulting parameter names
            par_name_list = par_name_list + name_list
            par_val_list = par_val_list + val_list

    # if scraping a list
    elif(type(var) is list):
        # browse
        for ind, this_item in enumerate(var):
            name_list, val_list = _scrap_data(this_item, "[" + str(ind) + "]")

            # add the resulting parameter names
            par_name_list = par_name_list + name_list
            par_val_list = par_val_list + val_list

    # if scraping an object instance of a class we wrote in the package
    elif(isinstance(var, (reco.MRSData2, fit.fit_stategy, sim.params))):
        # scrap the dict attribute
        name_list, val_list = _scrap_data(var.__dict__)

        # add the resulting parameter names and the original object
        par_name_list = par_name_list + ["obj"] + name_list
        par_val_list = par_val_list + [var] + val_list

    # if scraping an object instance of a class we wrote in the package
    elif(isinstance(var, (reco.pipeline, sim.mrs_sequence, fit.fit_stategy))):
        # scrap the dict attribute
        name_list, val_list = _scrap_data(var.__dict__)

        # add the resulting parameter names
        par_name_list = par_name_list + name_list
        par_val_list = par_val_list + val_list

    # if scraping an object instance of a class we wrote in the package
    elif(isinstance(var, (sim.params))):

        # add the resulting parameter names and the original object
        par_name_list = par_name_list + ["obj"]
        par_val_list = par_val_list + [var]

    # if scraping an optim result
    elif(isinstance(var, optimize.OptimizeResult)):
        # this is actually a dict but the __dict__ is empty for some reason (!?)
        name_list = list(var.keys())
        val_list = list(var.values())

        # add the resulting parameter names and the original object
        par_name_list = par_name_list + name_list
        par_val_list = par_val_list + val_list

    # if scraping some basis types
    elif(var is None or isinstance(var, (numbers.Number,
                                         str, datetime,
                                         sim.gating_signal_source,
                                         reco.data_rejection_method,
                                         tuple,
                                         np.ndarray,
                                         type))):
        # make it simple
        par_name_list = par_name_list + [""]
        par_val_list = par_val_list + [var]

    # if scraping something else, log it
    else:
        log.debug("not scraping variable %s" % str(type(var)))

    # format the field name, usefull later in pandas df
    if(prefix_str is not None):
        if(len(par_name_list) > 1 and type(var) is not list):
            prefix_str += "_"

        par_name_list = [prefix_str + s for s in par_name_list]

    return(par_name_list, par_val_list)
