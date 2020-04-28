#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A file to list all the default file paths used in this package. :)

@author: Tangi Roussel
"""

# used in sim.py for simulation
# full file path to metabolite database file
# it contains all the chemical shifts, J-couplings, etc. for each metabolite
DEFAULT_META_DB_FILE = "./mrs/metabolites_db.xls"

# full file path to metabolite basis set file
# it contains the metabolite names that will be used in simulation and fit
DEFAULT_META_BASIS_SET_FILE = "./metabolite_basis_sets/human_brain_7T.xls"

# pkl file where all the processed data is stored
DEFAULT_RECO_DATA_DB_FILE = "/home/tangir/crmbm/data_reco/db.pkl"
