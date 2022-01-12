#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A file to list all the default file paths used in this package. :)

@author: Tangi Roussel
"""

import os
pkg_path, _ = os.path.split(__file__)

# used in sim.py for simulation
# full file path to metabolite database file
# it contains all the chemical shifts, J-couplings, etc. for each metabolite
# and the metabolite basis set templates (proton by default)
DEFAULT_META_DB_FILE = os.path.join(pkg_path + "/metabolites_1H.xls")

# default metabolite basis set
# metabolite basis set are configured in the "Template_*" tabs of the above excel file
DEFAULT_META_BASIS_SET_NAME = "Brain_7T"

# pkl file where all the processed data is stored
DEFAULT_RECO_DATA_DB_FILE = "/home/tangir/crmbm/acq_db/db.pkl"
