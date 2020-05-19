#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A user script related to simulation.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.aliases as xxx
import mrs.reco as reco
import mrs.sim as sim
import mrs.log as log
import numpy as np
from datetime import datetime
import time
import os
import pickle
import copy

import pdb

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = reco.data_db()
# %% generate sequence database

# metabolite db
meta_bs = sim.metabolite_basis_set()
meta_bs.basis_set_xls_file = './metabolite_basis_sets/muscle_7T.xls'
meta_bs.non_coupled_only = True
meta_bs.initialize()

# generate smart filename
today = datetime.now()
fn = os.path.basename(meta_bs.basis_set_xls_file)
fn, _ = os.path.splitext(fn)
pickle_file_fullpath = "./" + fn + today.strftime("_%Y%m%d") + ".pkl"

# simulation ranges
seq_list = [sim.mrs_seq_press(20.0), sim.mrs_seq_steam(20.0), sim.mrs_seq_eja_svs_slaser(20.0)]
f0_list = [123.198361, 297.205620]  # MHz
te_list = np.arange(1.0, 400.0, 5.0)  # ms
seq_big_list = []

for this_seq in seq_list:
    for this_f0 in f0_list:
        for this_te in te_list:
            # display
            print("* Sequence = %s | f0 = %.0fMHz | TE = %.0fms" % (this_seq.name, this_f0, this_te))
            time.sleep(1.0)

            # setup sequence
            this_seq.te = this_te
            this_seq.f0 = this_f0
            try:
                # run
                this_seq.initialize(meta_bs)
                # store
                seq_big_list.append(copy.deepcopy(this_seq))
            except:
                log.warning("TE = %.0fms was too short! Skipped this sequence..." % this_te)

with open(pickle_file_fullpath, 'wb') as f:
    pickle.dump([seq_big_list], f)

# %% load the metabolite database and run a sLASER sequence

te = 55.0

meta_bs = sim.metabolite_basis_set()
meta_bs.initialize()

seq = sim.mrs_seq_eja_svs_slaser(te)
seq.pulse_rfc_r = 50.0
seq.pulse_rfc_optim_power_enable = True
seq.initialize(meta_bs)

# %% load the metabolite database and run a PRESS sequence

te = 55.0

meta_bs = sim.metabolite_basis_set()
meta_bs.initialize()

seq = sim.mrs_seq_eja_svs_press(te)
seq.initialize(meta_bs)

# %% load the metabolite database and run a STEAM sequence

te = 55.0

meta_bs = sim.metabolite_basis_set()
meta_bs.initialize()

seq = sim.mrs_seq_eja_svs_steam(te)
seq.initialize(meta_bs)

# %% spectrum for each metabolite
# left/right margins = 0.2/0.8 maximized over both screens for A3 format

p_human_min = sim.params(meta_bs).set_default_min()
p_human_max = sim.params(meta_bs).set_default_max()

p_human = (p_human_min + p_human_max) / 2.0
p_human[xxx.m_Water, 0] = 1.0
p_human[:, xxx.p_dd] = 10.0

ss_mod = seq.simulate_signal(p_human)
p_human.print()

fig = plt.figure(100)
fig.clf()
ax = fig.subplots()
fig.canvas.set_window_title("sim_datasets")
ax.set_title("GAMMA simulation of average normal brain MR spectrum (Cf. de Graaf's book) | " + seq.name + " | TE=" + str(seq.te) + "ms")
sim.disp_fit(ax, ss_mod, p_human, seq, True, True)
fig.subplots_adjust()

ss_mod.display_spectrum_1d(200, [1, 5])