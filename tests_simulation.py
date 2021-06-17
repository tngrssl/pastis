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
log.setLevel(log.INFO)

# %% generate sequence database

sim.GAMMA_LIB_LOADED = True

# simulation ranges
basis_set_name = 'Muscle'
basis_set_pkl_file = "20210331_muscle.pkl"
seq_list = [sim.mrs_seq_press, sim.mrs_seq_steam, sim.mrs_seq_eja_svs_slaser]
f0_list = [123.198361, 297.205620]  # MHz
te_list = np.arange(1.0, 400.0, 50.0)  # ms

# metabolite db
meta_bs = sim.metabolite_basis_set(basis_set_name=basis_set_name)

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
meta_bs = sim.metabolite_basis_set(basis_set_name=basis_set_name)

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

# %% testing simuation ppm range

sim.GAMMA_LIB_LOADED = True

te = 30.0

meta_bs = sim.metabolite_basis_set()

seq = sim.mrs_seq_press(te)
seq.bandpass_filter_range_ppm =[-2, 5]
seq.initialize(meta_bs)

p_test = ( sim.params(meta_bs).set_default_min() + sim.params(meta_bs).set_default_max() )/ 2
p_test[:, xxx.p_dd] = 10

ss = seq.simulate_signal(p_test)
p_test.print()

ss.display_spectrum_1d(display_range=[-5, 15])

# %% lipid spectrum: simple model with individual lipid components

# acquisition parameters
te = 12.0
scaling_factor = 1
noise_level = 0

# lipid parameters
CL = 20
ndb = 2
nmidb = 1

# load lipids metabolite basis set template
meta_bs = sim.metabolite_basis_set("Lipids")

# create a steam sequence
seq = sim.mrs_seq_steam(te)
seq.initialize(meta_bs)

# create a simulation parameter set
p = sim.params(meta_bs)
# null all concentrations
p[:, xxx.p_cm] = 0.0
# linewidth damping factor (Hz) for all
p[:, xxx.p_dd] = 10
# frequency shift (Hz) for all
p[:, xxx.p_df] = 0.0
# phase shift (rd) for all
p[:, xxx.p_dp] = 0.0

# add a strong water peak
p[xxx.m_Water, xxx.p_cm] = 50

# add the lipid amplitudes
p[xxx.m_LipA, xxx.p_cm] = 9
p[xxx.m_LipB, xxx.p_cm] = 6 * (CL - 4) - 8 * ndb + 2 * nmidb
p[xxx.m_LipC, xxx.p_cm] = 6
p[xxx.m_LipD, xxx.p_cm] = 4 * (ndb - nmidb)
p[xxx.m_LipE, xxx.p_cm] = 6
p[xxx.m_LipF, xxx.p_cm] = 2 * nmidb
p[xxx.m_LipG, xxx.p_cm] = 2
p[xxx.m_LipH, xxx.p_cm] = 2
p[xxx.m_LipI, xxx.p_cm] = 1
p[xxx.m_LipJ, xxx.p_cm] = 2 * ndb

# scale everything if needed
p[:, xxx.p_cm] = p[:, xxx.p_cm] * scaling_factor

# show final simulation parameter set
p.print()

# run simulation
s = seq.simulate_signal(p, noise_level)

# display and break down model into individual lipids
#fig = plt.figure(100)
#fig.clf()
#ax = fig.subplots()
#fit.disp_fit(ax, s, p, seq, display_range=[0, 6])

s.display_spectrum_1d(display_range=[0, 6], allowed_apodization=0)

# %% lipid spectrum: amplitude-linked model

# acquisition parameters
te = 12.0
scaling_factor = 1.0
noise_level = 0

# lipid parameters
CL = 17.5
ndb = 2.7
nmidb = 0.7

# load lipids metabolite basis set template
meta_bs = sim.metabolite_basis_set("Lipids")

# create a steam sequence
seq = sim.mrs_seq_steam(te)
seq.initialize(meta_bs)

# create a simulation parameter set
p = sim.params(meta_bs)
# null all concentrations
p[:, xxx.p_cm] = 0.0
# linewidth damping factor (Hz) for all
p[:, xxx.p_dd] = 100.0
# frequency shift (Hz) for all
p[:, xxx.p_df] = 0.0
# phase shift (rd) for all
p[:, xxx.p_dp] = 0.0

# add a strong water peak
p[xxx.m_Water, xxx.p_cm] = 50

# add the lipid amplitudes
p[xxx.m_LipA1, xxx.p_cm] = 1

p[xxx.m_LipB1, xxx.p_cm] = CL - 4
p[xxx.m_LipB2, xxx.p_cm] = 2 * ndb
p[xxx.m_LipB2, xxx.p_dp] = np.pi
p[xxx.m_LipB3, xxx.p_cm] = 2 * nmidb

p[xxx.m_LipC1, xxx.p_cm] = 1

p[xxx.m_LipD1, xxx.p_cm] = 2 * ndb
p[xxx.m_LipD2, xxx.p_cm] = 2 * nmidb
p[xxx.m_LipD2, xxx.p_dp] = np.pi

p[xxx.m_LipE1, xxx.p_cm] = 1
p[xxx.m_LipF1, xxx.p_cm] = 2 * nmidb
p[xxx.m_LipG1, xxx.p_cm] = 1
p[xxx.m_LipH1, xxx.p_cm] = 1
p[xxx.m_LipI1, xxx.p_cm] = 1
p[xxx.m_LipJ1, xxx.p_cm] = 2 * ndb

# scale everything if needed
p[:, xxx.p_cm] = p[:, xxx.p_cm] * scaling_factor

# forbid phase changes
p._linklock[:, xxx.p_dp] = 1

# same frequency shift for all
p._linklock[:, xxx.p_df] = 10
p._linklock[xxx.m_LipA1, xxx.p_df] = -10

# same damping for all
p._linklock[:, xxx.p_dd] = 20
p._linklock[xxx.m_LipA1, xxx.p_dd] = -20

# link the 2ndb amplitudes
p._linklock[[xxx.m_LipD1, xxx.m_LipJ1], xxx.p_cm] = 30
p._linklock[xxx.m_LipB2, xxx.p_cm] = -30

# link the 2nmib amplitudes
p._linklock[[xxx.m_LipD2, xxx.m_LipF1], xxx.p_cm] = 40
p._linklock[xxx.m_LipB3, xxx.p_cm] = -40

# link the other amplitudes
p._linklock[[xxx.m_LipC1, xxx.m_LipE1, xxx.m_LipG1, xxx.m_LipH1, xxx.m_LipI1], xxx.p_cm] = 50
p._linklock[xxx.m_LipA1, xxx.p_cm] = -50

# lock the simple lipids
p._linklock[[xxx.m_LipA, xxx.m_LipB, xxx.m_LipC, xxx.m_LipD, xxx.m_LipE, xxx.m_LipF, xxx.m_LipG, xxx.m_LipH, xxx.m_LipI, xxx.m_LipJ], xxx.p_cm] = 1

# check if the linklock looks good
p.check()

# unlink/link everything
p = p.toFullParams(p.toFreeParams())

# show final simulation parameter set
p.print()

# run simulation
s = seq.simulate_signal(p, noise_level)

# display and break down model into individual lipids
#fig = plt.figure(101)
#fig.clf()
#ax = fig.subplots()
#fit.disp_fit(ax, s, p, seq, display_range=[0, 6])

s.display_spectrum_1d(display_range=[0, 6])
