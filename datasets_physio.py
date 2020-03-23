#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2019
@author: Tangi Roussel
"""

from __future__ import division
import matplotlib.pylab as plt
import numpy as np
import mrs.reco as reco
from IPython import get_ipython
import warnings
warnings.filterwarnings("ignore", ".*GUI is implemented*")


# %% 27/08/2019 - 308-rs-p1-moelle - Ocha - short TR resp test

get_ipython().magic('clear')
get_ipython().magic('reset -f')

plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 8
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID177_steam_shortTE_SNR+_FID38922.dat
"""

p.display_legends = """
sLASER 20/1 SC short TR
"""

p.phase_enable = True
p.phase_display = False
p.recombine_phasing = True
p.realign_enable = False

p.analyse_and_reject_enable = True
p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-100, 0, -1, -3.14]
p.analyse_and_reject_max = [100, 50, 1, 3.14]

p.apodize_enable = False
p.apodize_damping_hz = 15
p.calibrate_enable = True

p.analyse_linewidth_enable = True
p.analyse_linewidth_magnitude_mode = False

p.analyse_snr_enable = True
p.analyse_snr_magnitude_mode = False

p.remove_water_enable = False
p.data_process_only_this_data_index = [0]

s = p.run_pipeline_std()


# %% 05/09/2019 - 311-sl-p1-moelle - Simon, test respiration

get_ipython().magic('clear')
get_ipython().magic('reset -f')


p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID90_steam_shortTE_SNR+_FID39702.dat
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID91_steam_shortTE_SNR+_FID39703.dat
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID115_steam_shortTE_SNR+_FID39727.dat
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID116_steam_shortTE_SNR+_FID39728.dat
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID117_steam_shortTE_SNR+_FID39729.dat
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID118_steam_shortTE_SNR+_FID39730.dat
"""

p.display_legends = """
STEAM noWS no trig
STEAM noWS resp trig 20%
STEAM noWS no trig
STEAM noWS resp trig 20%
STEAM noWS resp trig 10%
STEAM noWS resp trig 30%
"""

p.zerofill_npts = 4096*8

p.realign_enable = False
p.realign_POI_range_ppm = [4.5, 4.8]
p.realign_moving_averages = 1

p.analyse_and_reject_enable = False
p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-100, 0, -1, -3.14]
p.analyse_and_reject_max = [100, 20, 1, 3.14]
p.analyse_and_reject_POI_range_ppm = [4.5, 4.8]

p.apodize_enable = True
p.apodize_damping_hz = 5

p.calibrate_enable = False

p.analyse_linewidth_enable = True
p.analyse_linewidth_range_ppm = [4.5, 4.8]

p.analyse_snr_enable = True
p.analyse_snr_s_range_ppm = [4.5, 4.8]

p.remove_water_enable = False

p.data_process_only_this_data_index = [1, 2]
s = p.run_pipeline_std()


lbls = ['No trig.', 'No trig. post-corrected',
        'Respiration-triggered', 'Respiration-triggered post-corrected']
xpos = [0, 0.5, 1, 1.5]
plt.subplot(2, 1, 1)
plt.bar(xpos, [2101, 2163, 3227, 4309], width=0.3)
plt.xticks(xpos, lbls)
plt.ylabel('SNR (u.a)')
plt.grid('on')

plt.subplot(2, 1, 2)
plt.bar(xpos, [22.97, 17.09, 17.09, 9.77], width=0.3)
plt.xticks(xpos, lbls)
plt.ylabel('Water peak linewidth (Hz)')
plt.grid('on')
plt.legend(lbls)

plt.tight_layout()

# %% 23/09/2019 - 313-ft-p1-moelle - Fransiska, test respiration

get_ipython().magic('clear')
get_ipython().magic('reset -f')

plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 8
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/313-ft-p1-moelle/meas_MID60_steam_shortTE_SNR+_FID41492.dat
"""

p.data_physio_filepaths = """
/home/tangir/crmbm/acq_physio/313_FT_P1_MOELLE_1.resp
"""

p.display_legends = """
STEAM TR~1070ms NA=90
"""

p.phase_enable = False
p.phase_display = False
p.recombine_phasing = True

p.zerofill_enable = True
p.zerofill_npts = 4096*2

p.realign_enable = False
p.analyse_physio_enable = True
p.analyse_physio_POI_range_ppm = [4, 5]
p.analyse_physio_delta_time_ms = 40000.0

p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-50, 0, -0.1, -0.5]
p.analyse_and_reject_max = [50, 50, 0.1, 0.5]

p.apodize_enable = False
p.calibrate_enable = False
p.analyse_linewidth_enable = False
p.analyse_snr_enable = False
p.remove_water_enable = False

s = p.run_pipeline_std()


# %% 25/09/2019 - 314-yt-p1-moelle - Yolanda, test physio

get_ipython().magic('clear')
get_ipython().magic('reset -f')

plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 8
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID78_steam_shortTE_SNR+_FID41676.dat
"""

p.data_physio_filepaths = """
/home/tangir/crmbm/acq_physio/314_YT_P1_MOELLE_1.resp
"""

p.display_legends = """
sLASER 20/1 SC
"""

p.phase_enable = False
p.phase_display = False
p.recombine_phasing = True

p.zerofill_enable = True
p.zerofill_npts = 4096*4

p.realign_enable = False
p.realign_POI_range_ppm = [4.5, 4.8]
p.realign_moving_averages = 9

p.analyse_physio_enable = True
p.analyse_physio_POI_range_ppm = [4, 5]
p.analyse_physio_delta_time_ms = 40000.0

p.analyse_and_reject_enable = True
p.analyse_and_reject_moving_averages = 1
p.analyse_and_reject_min = [-50, 0, -0.1, -0.5]
p.analyse_and_reject_max = [50, 50, 0.1, 0.5]

p.apodize_enable = False
p.apodize_damping_hz = 10

p.calibrate_enable = False

p.analyse_linewidth_enable = False
p.analyse_linewidth_magnitude_mode = False

p.analyse_snr_enable = False
p.analyse_snr_magnitude_mode = False

p.remove_water_enable = False

s = p.run_pipeline_std()


# %% 03/10/2019 - 316-ap-p1-moelle - Anissa, test resp

get_ipython().magic('clear')
get_ipython().magic('reset -f')

plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 8
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID44_steam_shortTE_SNR+_FID42203.dat
"""

p.data_physio_filepaths = """
/home/tangir/crmbm/acq_physio/316_AP_P1_MOELLE.resp
"""

p.display_legends = """
STEAM TR=1010ms NA=256
"""

p.phase_enable = False
p.phase_display = False
p.recombine_phasing = True

p.zerofill_enable = True
p.zerofill_npts = 4096*2

p.realign_enable = False

p.analyse_and_reject_enable = True

p.analyse_physio_enable = True
p.analyse_physio_POI_range_ppm = [4, 5]
p.analyse_physio_delta_time_ms = 10000.0

p.apodize_enable = False
p.apodize_damping_hz = 10

p.calibrate_enable = False

p.analyse_linewidth_enable = False
p.analyse_linewidth_magnitude_mode = False

p.analyse_snr_enable = False
p.analyse_snr_magnitude_mode = False

p.remove_water_enable = False

s = p.run_pipeline_std()

# %% 05/11/2019 - 328-af-p1-moelle - Anne, test apnea/breath hold

get_ipython().magic('clear')
get_ipython().magic('reset -f')

plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 8
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/328-af-p1-moelle/meas_MID69_steam_shortTE_SNR+_FID45776.dat
"""

p.data_physio_filepaths = """
/home/tangir/crmbm/acq_physio/328_AF_P1_MOELLE.resp
"""

p.display_legends = """
STEAM TR=1010ms NA=300
"""

p.phase_enable = False
p.phase_display = False
p.recombine_phasing = True

p.zerofill_enable = True
p.zerofill_npts = 4096*2

p.realign_enable = False

p.analyse_and_reject_enable = True

p.analyse_physio_enable = True
p.analyse_physio_POI_range_ppm = [4, 5]
p.analyse_physio_delta_time_ms = 10000.0

p.apodize_enable = False
p.apodize_damping_hz = 10

p.calibrate_enable = False

p.analyse_linewidth_enable = False
p.analyse_linewidth_magnitude_mode = False

p.analyse_snr_enable = False
p.analyse_snr_magnitude_mode = False

p.remove_water_enable = False

s = p.run_pipeline_std()
