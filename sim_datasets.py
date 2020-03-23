#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A user script related to simulation.

@author: Tangi Roussel
"""

from __future__ import division
import matplotlib.pylab as plt
import mrs.metabase as xxx
import mrs.sim as sim

from IPython import get_ipython
import warnings
warnings.filterwarnings("ignore", ".*GUI is implemented*")
get_ipython().magic('clear')

# %% load the metabolite database and run a sLASER sequence

te = 55.0

meta_db = sim.metabolite_db()
meta_db.initialize()

seq = sim.mrs_seq_eja_svs_slaser(te)
seq.pulse_rfc_r = 50.0
seq.pulse_rfc_optim_power_enable = True
seq.initialize(meta_db)

# %% load the metabolite database and run a PRESS sequence

te = 55.0

meta_db = sim.metabolite_db()
meta_db.initialize()

seq = sim.mrs_seq_eja_svs_press(te)
seq.initialize(meta_db)

# %% load the metabolite database and run a STEAM sequence

te = 55.0

meta_db = sim.metabolite_db()
meta_db.initialize()

seq = sim.mrs_seq_eja_svs_steam(te)
seq.initialize(meta_db)

# %% spectrum for each metabolite
# left/right margins = 0.2/0.8 maximized over both screens for A3 format

p_human_min = sim.params(meta_db).set_default_human_brain_min()
p_human_max = sim.params(meta_db).set_default_human_brain_max()

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
plt.tight_layout()

ss_mod.display_spectrum(200, seq.name, [1, 5], 1.0, 0.0, False)
