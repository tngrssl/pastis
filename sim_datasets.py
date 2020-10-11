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
import mrs.fit as fit
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
meta_bs.one_proton_mode = True
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
# meta_bs["Cr_CH3"]["metabolites"]["Cr_CH3"]["ppm"] = [3.027, 3.027,3.027,3.027,3.027,3.027]
# meta_bs["Cr_CH3"]["metabolites"]["Cr_CH3"]["iso"] = [1,1,1,1,1,1]
# meta_bs["Cr_CH3"]["metabolites"]["Cr_CH3"]["J"] = np.zeros([6,6])

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
p_human[:, 0] = 0.0
p_human[:, xxx.p_dd] = 30.0
p_human[xxx.m_Cr_CH3, xxx.p_cm] = 1.0

ss_mod = seq.simulate_signal(p_human)
p_human.print()

pf = fit.prefit_tool(ss_mod, seq)
pf.area_integration_peaks = [xxx.m_Cr_CH3]
pf.initialize()
pf.run()

fig = plt.figure(100)
fig.clf()
ax = fig.subplots()
fig.canvas.set_window_title("sim_datasets")
ax.set_title("GAMMA simulation of average normal brain MR spectrum (Cf. de Graaf's book) | " + seq.name + " | TE=" + str(seq.te) + "ms")
fit.disp_fit(ax, ss_mod, p_human, seq, True, True)
fig.subplots_adjust()

ss_mod.display_spectrum_1d(200, [1, 5])

# %% ethanol tests

from mpl_toolkits.mplot3d import Axes3D

# metabolite db
meta_bs = sim.metabolite_basis_set()
meta_bs.basis_set_xls_file = "./metabolite_basis_sets/ethanol_7T.xls"
meta_bs.initialize()

# sequence
seq_effTE = sim.mrs_seq_press(88)
seq = sim.mrs_seq_eja_svs_slaser(90)
seq.allow_evolution_during_hard_pulses = False
seq.pulse_rfc_real_shape_enable = False
seq.pulse_rfc_r = 40
seq.pulse_rfc_optim_power_enable = True

# params
p = sim.params(meta_bs)
p[:] = 0.0
p[xxx.m_Eth, xxx.p_cm] = 1.0
p[xxx.m_Eth, xxx.p_dd] = 20

# test
seq.initialize(meta_bs)
seq_effTE.initialize(meta_bs)
plt.figure(1)
plt.clf()
s = seq_effTE.simulate_signal(p, 0.0, 1, "Original data")
s.display_spectrum_1d(1)
s = seq.simulate_signal(p, 0.0, 1, seq.name + " TE=" + str(seq.te) + "ms")
s.display_spectrum_1d(1)

# %% play with TE

te_list = np.arange(50, 250, 1)
s_list = []
for te in te_list:
    seq.te = te
    s = seq.simulate_signal(p)
    s_list.append(s.spectrum())

# wireframe plot
te2d, ppm2d = np.meshgrid(s.frequency_axis_ppm(), te_list)
s_list_np = np.array(s_list)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(te2d, ppm2d, s_list_np)
ax.set_ylim(te_list[0], te_list[-1])
ax.set_xlim(1, 2)
plt.show()

# look for the middle multiplet peak at 1.3ppm
ppm = s.frequency_axis_ppm()
peak_index, _, peak_val, _, _, _ = s._analyze_peak_1d([1, 2])
# build J evolution curve
peak_te_evol = []
for s in s_list:
    peak_val = np.real(s.spectrum()[peak_index])
    peak_te_evol.append(peak_val)

peak_te_evol = np.array(peak_te_evol)

# evaluate J-modulation amplitude relative to first TE
peak_te_evol_rel = peak_te_evol / peak_te_evol[0] * 100.0

fig = plt.figure(2)
plt.plot(te_list, peak_te_evol_rel)




