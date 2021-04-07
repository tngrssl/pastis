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
import mrs.sim as sim
import mrs.fit as fit
import mrs.log as log
import numpy as np
import pandas as pd
import time

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

# %% generate sequence database

sim.GAMMA_LIB_LOADED = True

# simulation ranges
basis_set_name = 'Muscle'
basis_set_pkl_file = "20210331_muscle.pkl"
seq_list = [sim.mrs_seq_press, sim.mrs_seq_steam, sim.mrs_seq_eja_svs_slaser]
f0_list = [123.198361, 297.205620]  # MHz
te_list = np.arange(1.0, 400.0, 50.0)  # ms

# metabolite db
meta_bs = sim.metabolite_basis_set(basis_set_name=basis_set_name, non_coupled_only=True, one_proton_mode=True)

df_list = []
for this_seq in seq_list:
    for this_f0 in f0_list:
        for this_te in te_list:
            # setup sequence
            this_seq_obj = this_seq(te=this_te, f0=this_f0)

            # run
            try:
                this_seq_obj.initialize(meta_bs)
            except:
                print("* Sequence = %s | f0 = %.0fMHz | TE = %.0fms" % (this_seq_obj.name, this_f0, this_te))
                log.warning("Skipped this sequence... Maybe TE (%.0fms) was too short?" % this_te)
                time.sleep(10.0)

            # store
            df_list.append(this_seq_obj.to_dataframe(True))
            # display
            print("* Sequence = %s | f0 = %.0fMHz | TE = %.0fms" % (this_seq_obj.name, this_f0, this_te))
            time.sleep(1.0)

# build dataframe and store on disk
df = pd.concat(df_list, axis=0)
df.to_pickle(basis_set_pkl_file)

# %% and test

# pretend I cannot load pyGAMMA
sim.GAMMA_LIB_LOADED = False

te = 55.0
basis_set_name = 'Muscle'
basis_set_pkl_file = "20210331_muscle.pkl"
meta_bs = sim.metabolite_basis_set(basis_set_name=basis_set_name, non_coupled_only=True, one_proton_mode=True)

seq = sim.mrs_seq_press(te)
seq.db_file = basis_set_pkl_file
seq.initialize(meta_bs)

# %% load the metabolite database and run a sLASER sequence

sim.GAMMA_LIB_LOADED = True

te = 55.0

meta_bs = sim.metabolite_basis_set()

seq = sim.mrs_seq_eja_svs_slaser(te)
seq.pulse_rfc_r = 50.0
seq.pulse_rfc_optim_power_enable = True
seq.initialize(meta_bs)

# %% spectrum for each metabolite
# left/right margins = 0.2/0.8 maximized over both screens for A3 format

p_human_min = sim.params(meta_bs).set_default_min()
p_human_max = sim.params(meta_bs).set_default_max()
p_human = (p_human_min + p_human_max) / 2.0

ss_mod = seq.simulate_signal(p_human)
p_human.print()

fig = plt.figure(100)
fig.clf()
ax = fig.subplots()
fig.canvas.set_window_title("sim_datasets")
ax.set_title("pyGAMMA simulation of MR spectrum\nTemplate=" + meta_bs.basis_set_name + " | Sequence=" + seq.name + " | TE=" + str(seq.te) + "ms")
fit.disp_fit(ax, ss_mod, p_human, seq, True, True)
