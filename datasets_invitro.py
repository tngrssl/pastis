#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vitro data.

@author: Tangi Roussel
"""
# %% init
from __future__ import division
import matplotlib.pylab as plt
import mrs.metabase as xxx
import numpy as np
import mrs.reco as reco

from IPython import get_ipython
import warnings
warnings.filterwarnings("ignore", ".*GUI is implemented*")
plt.close('all')
get_ipython().magic('clear')

# %% 27/02/2019 - spinal cord phantom - nice WS STEAM, testing TWIX reco

get_ipython().magic('clear')
plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/braino_009/meas_MID49_svs_st_vapor_643_FID28469.dat
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0024_svs-st-vapor-643
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/braino_009/meas_MID50_svs_st_vapor_643_FID28470.dat
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0025_svs-st-vapor-643
"""

p.display_legends = """
TWIX
DICOM
"""

p.display_amp_factor_list = [57412 / 0.000692, 1]
p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 27/02/2019 - spinal cord phantom - OVS optim

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0026_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0027_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0028_svs-st-vapor-643
"""

p.display_legends = """
original - OVS duration 3000µs
OVS duration 1000µs
OVS duration 5000µs
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 27/02/2019 - spinal cord phantom - spoiler optim

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0029_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0030_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0031_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0032_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0033_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0034_svs-st-vapor-643
"""

p.display_legends = """
original - spoiler 20mT/m for 500µs
spoiler 20mT/m for 250µs
spoiler 20mT/m for 100µs
spoiler 10mT/m for 500µs
spoiler 5mT/m for 500µs
spoiler 0mT/m for 500µs
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 27/02/2019 - spinal cord phantom - various optim

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0035_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0036_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0037_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0038_svs-st-vapor-643
"""

p.display_legends = """
original
asym->sym pulses
2048->4096 pts
0->200µs acq win shift
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 27/02/2019 - spinal cord phantom - TM optim

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0039_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0040_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0041_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0042_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0043_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0044_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0045_svs-st-vapor-643
"""

p.display_legends = """
original - TM=43ms
TM=50ms
TM=60ms
TM=70ms
TM=80ms
TM=90ms
TM=100ms
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 27/02/2019 - spinal cord phantom - compare svs-st-vapor-643, eja-svs-steam and eja_svs_slaser

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0039_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0051_eja-svs-steam
/home/tangir/crmbm/acq/braino_009/braino-phantom/20190227/01_0079_eja-svs-slaser
"""

p.display_legends = """
svs-st-vapor-643
eja-svs-steam
eja-svs-slaser
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 12/03/2019 - spinal cord phantom - TR optim

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0026_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0027_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0028_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0029_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0030_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0031_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0032_svs-st-vapor-643
/home/tangir/crmbm/acq/braino_010/braino-phantom/20190312/01_0033_svs-st-vapor-643
"""

p.display_legends = """
TR=1500ms
TR=2000ms
TR=2500ms
TR=3000ms
TR=3500ms
TR=4000ms
TR=4500ms
TR=5000ms
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0016_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0019_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0022_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0025_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0028_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0031_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0034_svs-st-vapor-643-optim-trig
"""

p.display_legends = """
initial STEAM
sym. pulses
increased settling time
increasing ramp time
short strong spoilers
mini spoiler
OVS on
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - VOI traces

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0018_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0021_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0024_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0027_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0030_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0033_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0036_svs-st-vapor-643-optim-trig
"""

p.display_legends = """
initial STEAM
sym. pulses
increased settling time
increasing ramp time
short strong spoilers
mini spoiler
OVS on
"""

p.run_pipeline_std()


# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - sLASER

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0040_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0044_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0047_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0051_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0054_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0057_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0060_eja-svs-slaser-optim-trig
"""

p.display_legends = """
initial sLASER
OVS on
short strong spoilers
reducing ramp time
increasing settling time
reducing R 20->10
increasing R 20->25
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM versus sLASER

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0016_svs-st-vapor-643-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0040_eja-svs-slaser-optim-trig
"""

p.display_legends = """
initial STEAM
initial sLASER
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 20/03/2019 - fatty_braino_002 - sLASER optim invitro in invivo B1/B0 conditions - VOI 1/2 profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0043_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0050_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0053_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0056_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0059_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0062_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0063_eja-svs-slaser-optim-trig
"""

p.display_legends = """
initial sLASER
OVS on
short strong spoilers
reducing ramp time
increasing settling time
reducing R 20->10
increasing R 20->25
"""

p.run_pipeline_std()


# %% 20/03/2019 - fatty_braino_002 - sLASER optim invitro in invivo B1/B0 conditions, R factor - VOI 1/2 profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0063_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0064_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0065_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0066_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0067_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0068_eja-svs-slaser-optim-trig
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0069_eja-svs-slaser-optim-trig
"""

p.display_legends = """
increasing R=5
increasing R=10
increasing R=15
increasing R=20
increasing R=25
increasing R=30
increasing R=40
"""

p.run_pipeline_std()


# %%20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - testing TWIX reco (1/2)

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/fatty_braino_002/meas_MID54_svs_st_vapor_643_optim_trig_FID29713.dat
"""

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/fatty_braino_002/meas_MID53_svs_st_vapor_643_optim_trig_FID29712.dat
"""

p.display_legends = """
twix - ref 0th order phasing
"""

p.display_amp_factor_list = [1 / 0.00005 * 1 / 1.2]
p.phase_display = False
p.phase_high_snr_mode = False
p.phase_weak_ws_mode = False
p.phase_order = 0
p.remove_water_enable = False
p.run_pipeline_std()


# %%20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - compare previous TWIX with DICOM (2/2)

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0017_svs-st-vapor-643-optim-trig
"""

p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_002/fatty-braino/20190320/01_0016_svs-st-vapor-643-optim-trig
"""

p.display_legends = """
DICOM - 0th order phasing with REF
"""

p.display_amp_factor_list = [1 / 6650]
p.phase_display = False
p.phase_high_snr_mode = False
p.phase_weak_ws_mode = False
p.phase_order = 0
p.remove_water_enable = False
p.run_pipeline_std()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - sym vs asym STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0026_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0024_svs-st-vapor-643-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0027_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0025_svs-st-vapor-643-optim
"""

p.display_legends = """
asym STEAM
sym STEAM
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - sym vs asym STEAM - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0021_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0022_svs-st-vapor-643-optim
"""

p.display_legends = """
sym STEAM
asym sSTEAM
"""

p.analyze_selectivity_range_list = [[-2000, 1600], [-12500, -9000], [-3800, 2600]]
p.run_pipeline_std()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - settling time in STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0028_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0030_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0032_svs-st-vapor-643-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0029_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0031_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0033_svs-st-vapor-643-optim
"""

p.display_legends = """
STEAM 500us settling
STEAM 300us settling
STEAM 10us settling"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - settling time in STEAM - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0034_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0035_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_003/fatty-braino/20190503/01_0036_svs-st-vapor-643-optim
"""

p.display_legends = """
STEAM 10us settling
STEAM 300us settling
STEAM 500us settling
"""

p.run_pipeline_std()


# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - OVS in STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0006_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0008_svs-st-vapor-643-optim/
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0007_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0009_svs-st-vapor-643-optim
"""

p.display_legends = """
STEAM no OVS
STEAM OVS 45deg (SAR limit)
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - OVS in STEAM - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0012_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0013_svs-st-vapor-643-optim/
"""

p.display_legends = """
no OVS
OVS 45deg
"""

p.run_pipeline_std()


# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - ramp time in STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0014_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0016_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0018_svs-st-vapor-643-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0015_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0017_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0019_svs-st-vapor-643-optim
"""

p.display_legends = """
200us ramp
100us ramp
300us ramp
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - OVS in STEAM - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0020_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0021_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0022_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0023_svs-st-vapor-643-optim
"""

p.display_legends = """
200 us ramp
100 us ramp
300 us ramp
100 us ramp FOV 500
"""

p.run_pipeline_std()


# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - spoilers in STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0024_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0027_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0028_svs-st-vapor-643-optim/
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0025_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0026_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0029_svs-st-vapor-643-optim
"""

p.display_legends = """
20mT/m - 500ms
35mT/m - 100ms
20mT/m - 100ms
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - spoilers in STEAM - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0030_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0031_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0032_svs-st-vapor-643-optim
"""

p.display_legends = """
20mT/m - 500ms
35mT/m - 100ms
20mT/m - 100ms
"""

p.run_pipeline_std()


# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - VAPOR BW in STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0033_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0036_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0038_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0043_svs-st-vapor-643-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0034_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0037_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0039_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0044_svs-st-vapor-643-optim
"""

p.display_legends = """
135Hz
200Hz
250Hz
optim
"""

p.phase_display = False
p.remove_water_enable = False
p.run_pipeline_std()


# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - testing 30/30/60 OVS in STEAM - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0045_svs-st-vapor-643-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0046_svs-st-vapor-643-optim
"""

p.display_legends = """
30/30/60 OVS
no OVS
"""

p.run_pipeline_std()


# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - R/pulse length in sLASER - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0051_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0052_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0053_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0054_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0055_eja-svs-slaser-optim
"""

p.display_legends = """
R=20
R=10
OVS 90
+ 30mT/m spoiler"
+ 40mT/m spoiler
"""

p.run_pipeline_std()


# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - R/N pulses in sLASER - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0056_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0057_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0058_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0059_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0060_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0061_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0062_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0063_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0064_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0065_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0066_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0067_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0068_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_004/fatty-braino/20190507/01_0069_eja-svs-slaser-optim
"""

p.display_legends = """
                N=1,R=5
                N=1,R=10
                N=1,R=15
                N=1,R=20
                N=1,R=25
                N=1,R=30
                N=2,R=5
                N=2,R=10
                N=2,R=15
                N=3,R=5
                N=3,R=10
                N=4,R=5
                N=4,R=10
                N=5,R=5
                """

p.analyze_selectivity_range_list = [
    [800, 3550], [-10600, -7800], [-3650, 1450]]

p.run_pipeline_std()

# more info
n_r = np.array([[1, 5],
                [1, 10],
                [1, 15],
                [1, 20],
                [1, 25],
                [1, 30],
                [2, 5],
                [2, 10],
                [2, 15],
                [3, 5],
                [3, 10],
                [4, 5],
                [4, 10],
                [5, 5]])


min_TE = np.array([[24],
                   [32],
                   [36],
                   [44],
                   [44],
                   [48],
                   [18],
                   [24],
                   [28],
                   [16],
                   [22],
                   [16],
                   [22],
                   [16]])

# this is shaped [14,3,2]: 14 scans, 3 axes, 2 conditions IN/OUT
sel_results = p.analyze_selectivity_list

# reshaping to [5,7,3,2]: 5 possible Ns, 7 possible Rs, 3 axes, 2 conditions IN/OUT
sel_results_2d = np.full([5, 7, 3, 2], np.nan)
min_TE2d = np.full([5, 7], np.nan)
for k in range(0, sel_results.shape[0]):
    this_n = n_r[k, 0]
    this_r = n_r[k, 1]
    this_n_ind = int(this_n - 1)
    this_r_ind = int(this_r / 5 - 1)
    sel_results_2d[this_n_ind, this_r_ind, :, :] = sel_results[k, :, :]
    min_TE2d[this_n_ind, this_r_ind] = min_TE[k]

# and plotting
fig, ax = plt.subplots(2, 3)

ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 0]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 0]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 0]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 4)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 1]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 5)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 1]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 6)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 1]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

# ratio
sel_results_2d_ratio = np.full([5, 7, 3], np.nan)
sel_results_2d_ratio = sel_results_2d[:, :, :, 0] / sel_results_2d[:, :, :, 1]

fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 0]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 1]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 2]))
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

# and the min TE heat map
fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(min_TE2d)
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('Mininum TE')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('TE (ms)')
plt.show()


# %% 07/05/2019 - fatty_braino_005, bad B1/B0 conditions - R/N pulses in sLASER TE fixed at 50ms

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0005_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0007_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0009_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0011_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0013_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0015_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0017_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0019_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0021_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0023_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0025_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0027_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0029_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0031_eja-svs-slaser-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0006_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0008_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0010_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0012_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0014_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0016_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0018_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0020_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0022_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0024_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0026_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0028_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0030_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_005/fatty-braino/20190509/01_0032_eja-svs-slaser-optim
"""

p.display_legends = """
                N=1,R=5
                N=1,R=10
                N=1,R=15
                N=1,R=20
                N=1,R=25
                N=1,R=30
                N=2,R=5
                N=2,R=10
                N=2,R=15
                N=3,R=5
                N=3,R=10
                N=4,R=5
                N=4,R=10
                N=5,R=5
                """

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = True

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.run_pipeline_std()

# more info
n_r = np.array([[1, 5],
                [1, 10],
                [1, 15],
                [1, 20],
                [1, 25],
                [1, 30],
                [2, 5],
                [2, 10],
                [2, 15],
                [3, 5],
                [3, 10],
                [4, 5],
                [4, 10],
                [5, 5]])

# this is shaped [14]: 14 scans
snr_results = p.analyze_snr_list

# reshaping to [5,7]: 5 possible Ns, 7 possible Rs
snr_results_2d = np.full([5, 7], np.nan)
for k in range(0, len(snr_results)):
    this_n = n_r[k, 0]
    this_r = n_r[k, 1]
    this_n_ind = int(this_n - 1)
    this_r_ind = int(this_r / 5 - 1)
    snr_results_2d[this_n_ind, this_r_ind] = snr_results[k]

# and plotting
fig, ax = plt.subplots()

ax = plt.subplot(2, 2, 1)
im = ax.imshow(snr_results_2d)
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('SNR on NAA peak')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('SNR')
plt.show()


# %% 16/05/2019 - fatty_braino_007, bad B1/B0 conditions - R/N pulses in sLASER TE minimal

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0006_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0008_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0010_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0012_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0014_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0016_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0018_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0020_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0022_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0024_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0026_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0028_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0030_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0032_eja-svs-slaser-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0007_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0009_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0011_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0013_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0015_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0017_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0019_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0021_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0023_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0025_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0027_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0029_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0031_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0033_eja-svs-slaser-optim
"""

p.display_legends = """
                N=1,R=5
                N=1,R=10
                N=1,R=15
                N=1,R=20
                N=1,R=25
                N=1,R=30
                N=2,R=5
                N=2,R=10
                N=2,R=15
                N=3,R=5
                N=3,R=10
                N=4,R=5
                N=4,R=10
                N=5,R=5
                """

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = True

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.run_pipeline_std()

# more info
n_r = np.array([[1, 5],
                [1, 10],
                [1, 15],
                [1, 20],
                [1, 25],
                [1, 30],
                [2, 5],
                [2, 10],
                [2, 15],
                [3, 5],
                [3, 10],
                [4, 5],
                [4, 10],
                [5, 5]])

# this is shaped [14]: 14 scans
snr_results = p.analyze_snr_list

# reshaping to [5,7]: 5 possible Ns, 7 possible Rs
snr_results_2d = np.full([5, 7], np.nan)
for k in range(0, len(snr_results)):
    this_n = n_r[k, 0]
    this_r = n_r[k, 1]
    this_n_ind = int(this_n - 1)
    this_r_ind = int(this_r / 5 - 1)
    snr_results_2d[this_n_ind, this_r_ind] = snr_results[k]

# and plotting
fig, ax = plt.subplots()

ax = plt.subplot(2, 2, 1)
im = ax.imshow(snr_results_2d)
ax.set_xticks(range(6))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35])
ax.set_xlabel('R')
ax.set_yticks(range(5))
ax.set_yticklabels([1, 2, 3, 4, 5])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('SNR on NAA peak')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('SNR')
plt.show()


# %%16/05/2019 - fatty_braino_007, bad B1/B0 conditions - STEAM - testing TWIX reco

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID31_eja_svs_slaser_optim_FID32220.dat
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0033_eja-svs-slaser-optim
"""

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID30_eja_svs_slaser_optim_FID32219.dat
/home/tangir/crmbm/acq/fatty_braino_007/fatty-braino/20190516/01_0032_eja-svs-slaser-optim
"""

p.display_legends = """
twix 0th
dcm 0th
"""

p.phase_display = True
p.phase_enable = True
p.phase_ref_signal_too = True
p.phase_high_snr_mode = False
p.phase_weak_ws_mode = False
p.phase_order = 0

p.realign_enable = False
p.analyze_snr_enable = True
p.display_amp_factor_list = [80000000, 1]

p.remove_water_enable = False
p.run_pipeline_std()


# %% 21/05/2019 - fatty_braino_007, bad B1/B0 conditions - R/N pulses in sLASER TE minimal

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0008_svs-st-vapor-643
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0015_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0018_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0020_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0022_eja-svs-slaser-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0009_svs-st-vapor-643
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0017_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0019_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0021_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0023_eja-svs-slaser-optim
"""

p.display_legends = """
                steam
                N=1,R=5
                N=1,R=20
                N=2,R=10
                N=5,R=5
                """

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = True

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.run_pipeline_std()


# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions - R/N pulses in sLASER TE minimal

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0026_svs-st-vapor-643
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0030_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0032_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0034_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0036_eja-svs-slaser-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0027_svs-st-vapor-643
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0031_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0033_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0035_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0037_eja-svs-slaser-optim
"""

p.display_legends = """
                steam
                N=1,R=5
                N=1,R=20
                N=2,R=10
                N=5,R=5
                """

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = True

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.run_pipeline_std()


# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions, stronger spoilers - R/N pulses in sLASER TE minimal

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0026_svs-st-vapor-643
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0041_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0043_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0045_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0047_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0049_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0051_eja-svs-slaser-optim
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0027_svs-st-vapor-643
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0042_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0044_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0046_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0048_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0050_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0052_eja-svs-slaser-optim
"""

p.display_legends = """
steam
N=1,R=5
N=1,R=20
N=2,R=10
N=5,R=5
N=1,R=20, big spoil
N=5,R=5, big spoil
"""

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = True

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.run_pipeline_std()


# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions, even more stronger spoilers - R/N pulses in sLASER TE minimal

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0053_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0055_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0057_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0059_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0061_svs-st-vapor-643
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0054_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0056_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0058_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0060_eja-svs-slaser-optim
/home/tangir/crmbm/acq/fatty_braino_008/fatty-braino/20190521/01_0062_svs-st-vapor-643
"""

p.display_legends = """
                N=1,R=5
                N=1,R=20
                N=2,R=10
                N=5,R=5
                steam
                """

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = True

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.run_pipeline_std()


# %% 03/06/2019 - fatty_braino_head_coil_002, test non-WS steam

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0010_steam-shortte-snr
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0012_steam-shortte-snr
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0011_steam-shortte-snr
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0012_steam-shortte-snr
"""

p.display_legends = """
                steam
                steam nonWS
                """

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = True

p.display_range_ppm = [1, 6]

p.remove_water_enable = True
p.run_pipeline_std()


# %% 05/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - VOI profiles

get_ipython().magic('clear')


plt.close("all")

p = reco.voi_pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0040_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0037_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0034_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0031_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0027_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0024_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0021_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0018_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0051_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0068_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0077_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0052_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0071_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0053_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0074_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0056_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0059_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0062_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0065_slaser-r-n
"""

p.display_legends = """
                N=1,R=5
                N=1,R=10
                N=1,R=15
                N=1,R=20
                N=1,R=25
                N=1,R=30
                N=1,R=35
                N=1,R=40
                N=2,R=5
                N=2,R=10
                N=2,R=15
                N=3,R=5
                N=3,R=10
                N=4,R=5
                N=4,R=10
                N=5,R=5
                N=6,R=5
                N=7,R=5
                N=8,R=5
                """

p.analyze_selectivity_range_list = [
    [-6000, -2000], [-9000, -5000], [-5000, -2000]]

p.run_pipeline_std()

# more info
n_r = np.array([[1, 5],
                [1, 10],
                [1, 15],
                [1, 20],
                [1, 25],
                [1, 30],
                [1, 35],
                [1, 40],
                [2, 5],
                [2, 10],
                [2, 15],
                [3, 5],
                [3, 10],
                [4, 5],
                [4, 10],
                [5, 5],
                [6, 5],
                [7, 5],
                [8, 5]])


min_TE = p.get_te_list()

# this is shaped [20,3,2]: 20 scans, 3 axes, 2 conditions IN/OUT
sel_results = p.analyze_selectivity_list

# reshaping to [5,8,3,2]: 8 possible Ns, 8 possible Rs, 3 axes, 2 conditions IN/OUT
sel_results_2d = np.full([8, 8, 3, 2], np.nan)
min_TE2d = np.full([8, 8], np.nan)
for k in range(sel_results.shape[0]):
    this_n = n_r[k, 0]
    this_r = n_r[k, 1]
    this_n_ind = int(this_n - 1)
    this_r_ind = int(this_r / 5 - 1)
    sel_results_2d[this_n_ind, this_r_ind, :, :] = sel_results[k, :, :]
    min_TE2d[this_n_ind, this_r_ind] = min_TE[k]

# and plotting
fig, ax = plt.subplots(2, 3)

ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 0]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 0]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 0]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 4)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 0, 1]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 5)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 1, 1]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

ax = plt.subplot(2, 3, 6)
im = ax.imshow(np.squeeze(sel_results_2d[:, :, 2, 1]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('profile area ($Hz^{-1}$)')
plt.show()

# ratio
sel_results_2d_ratio = np.full([8, 8, 3], np.nan)
sel_results_2d_ratio = sel_results_2d[:, :, :, 0] / sel_results_2d[:, :, :, 1]

fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 0]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[X] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 2)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 1]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Y] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

ax = plt.subplot(2, 3, 3)
im = ax.imshow(np.squeeze(sel_results_2d_ratio[:, :, 2]))
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('[Z] IN/OUT')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('Ratio IN/OUT (u.a)')
plt.show()

# and the min TE heat map
fig, ax = plt.subplots(2, 3)
ax = plt.subplot(2, 3, 1)
im = ax.imshow(min_TE2d)
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('Mininum TE')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('TE (ms)')
plt.show()


# %% 05/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - sLASER spectra

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0038_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0035_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0032_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0028_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0025_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0022_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0019_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0016_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0041_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0066_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0075_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0045_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0069_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0048_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0072_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0054_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0057_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0060_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0063_slaser-r-n
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0039_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0036_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0033_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0029_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0026_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0023_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0020_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0017_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0042_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0067_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0076_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0046_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0070_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0049_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0073_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0055_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0058_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0061_slaser-r-n
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0064_slaser-r-n
"""

p.display_legends = """
                N=1,R=5
                N=1,R=10
                N=1,R=15
                N=1,R=20
                N=1,R=25
                N=1,R=30
                N=1,R=35
                N=1,R=40
                N=2,R=5
                N=2,R=10
                N=2,R=15
                N=3,R=5
                N=3,R=10
                N=4,R=5
                N=4,R=10
                N=5,R=5
                N=5,R=6
                N=5,R=7
                N=5,R=8
                """

p.display_offset = 30000
p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = False

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.run_pipeline_std()

# more info
n_r = np.array([[1, 5],
                [1, 10],
                [1, 15],
                [1, 20],
                [1, 25],
                [1, 30],
                [1, 35],
                [1, 40],
                [2, 5],
                [2, 10],
                [2, 15],
                [3, 5],
                [3, 10],
                [4, 5],
                [4, 10],
                [5, 5],
                [6, 5],
                [7, 5],
                [8, 5]])

# this is shaped [14]: 14 scans
snr_results = p.analyze_snr_list

# reshaping to [8,8]: 8 possible Ns, 8 possible Rs
snr_results_2d = np.full([8, 8], np.nan)
for k in range(0, len(snr_results)):
    this_n = n_r[k, 0]
    this_r = n_r[k, 1]
    this_n_ind = int(this_n - 1)
    this_r_ind = int(this_r / 5 - 1)
    snr_results_2d[this_n_ind, this_r_ind] = snr_results[k]

# and plotting
fig, ax = plt.subplots()

ax = plt.subplot(2, 2, 1)
im = ax.imshow(snr_results_2d)
ax.set_xticks(range(8))
ax.set_xticklabels([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_xlabel('R')
ax.set_yticks(range(8))
ax.set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8])
ax.set_ylabel('N')
ax.grid('on')
ax.set_title('SNR on NAA peak')
cb = fig.colorbar(ax=ax, mappable=im, orientation='vertical')
cb.set_label('SNR')
plt.show()


# %% 05/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - STEAM WS spectra

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0010_steam-shortte-snr
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/fatty_braino_head_coil_002/fatty-braino/20190603/01_0012_steam-shortte-snr
"""

p.display_legends = """
                steam WS
                """

p.phase_enable = True
p.phase_display = False
p.phase_order = 0

p.analyze_snr_s_range_ppm = [1.9, 2.1]
p.analyze_snr_area_integrate = False

p.display_range_ppm = [1, 6]

p.remove_water_enable = False
p.remove_water_hsvd_components = 1
p.remove_water_hsvd_pars = [70, 100]
p.analyze_snr_enable = False
p.calibrate_enable = False
p.run_pipeline_std()


# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing non WS STEAM spectra

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0012_steam-shortte-snr
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr
"""

p.data_ref_filepaths = """
"""

p.display_legends = """
                steam ws 90
                steam ws 30deg
                steam ws 1deg
                """

p.phase_enable = False

p.analyze_snr_enable = True
p.analyze_snr_s_range_ppm = [1.5, 2]
p.analyze_snr_area_integrate = False

p.display_range_ppm = [1, 6]

p.remove_water_enable = True
p.remove_water_hsvd_components = 3
p.remove_water_hsvd_range = [4.6, 4.8]

p.calibrate_enable = False
p.run_pipeline_std()


# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing DCM reco with right ref scans...

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr
"""

p.display_legends = """
                steam dcm not phased
                """

p.phase_enable = False
p.phase_order = 0

p.display_range_ppm = [1, 6]
p.remove_water_enable = False

p.run_pipeline_std()


# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing TWIX reco with right ref scans...

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr
/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID90_steam_shortTE_SNR+_FID33706.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr
/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID92_steam_shortTE_SNR+_FID33708.dat
"""

p.display_legends = """
dcm phased
twix phased
"""

p.display_amp_factor_list = [1, 1e8 * 50800 / 26300]

p.phase_enable = True
p.phase_display = True

p.apodize_enable = True
p.apodize_damping_hz = 20

p.run_pipeline_std()


# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing TWIX reco and non WS STEAM

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID91_steam_shortTE_SNR+_FID33707.dat
/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID92_steam_shortTE_SNR+_FID33708.dat
"""

p.data_ref_filepaths = """
"""

p.display_legends = """
                VAPOR FA=30deg
                VAPOR FA=1deg
                """
p.display_amp_factor_list = [1, 1]

p.phase_enable = True
p.phase_display = False
p.phase_weak_ws_mode = False
p.phase_POI_range_ppm = [4, 5]
p.phase_order = 1

p.recombine_phasing = True
p.apodize_enable = True
p.apodize_damping_hz = 20

p.display_range_ppm = [1, 6]
p.remove_water_enable = True

p.run_pipeline_std()


# %%16/05/2019 - fatty_braino_007, bad B1/B0 conditions - STEAM - testing TWIX reco

get_ipython().magic('clear')


plt.close("all")

for chan_to_turn_off in range(0, 9):

    print("chan_to_turn_off=" + str(chan_to_turn_off))
    p = reco.pipeline()
    p.data_ref_filepaths = """
    /home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID31_eja_svs_slaser_optim_FID32220.dat
    """

    p.data_filepaths = """
    /home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID30_eja_svs_slaser_optim_FID32219.dat
    """

    p.display_legends = "off channel = " + str(chan_to_turn_off)
    p.display_amp_factor_list = [1]

    p.phase_enable = True
    p.phase_display = False
    p.phase_order = 1

    p.recombine_phasing = False
    p.apodize_enable = False

    p.display_range_ppm = [1, 6]
    p.remove_water_enable = False

    p.recombine_weights = np.full([8, ], True)
    if(chan_to_turn_off < 8):
        p.recombine_weights[chan_to_turn_off] = False

    p.run_pipeline_std()


# %%27/11/2019 - phantom_test_inv, sLASER - inv test

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 8
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID145_slaser_R_N=20+_1_longTE_SNR++++_FID47558.dat
/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID149_slaser_R_N=20+_1_longTE_SNR++++_FID47562.dat
/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID151_slaser_R_N=20+_1_longTE_SNR++++_FID47564.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID147_slaser_R_N=20+_1_longTE_SNR++++_FID47560.dat
/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID147_slaser_R_N=20+_1_longTE_SNR++++_FID47560.dat
/home/tangir/crmbm/acq_twix/phantom_test_inv/meas_MID147_slaser_R_N=20+_1_longTE_SNR++++_FID47560.dat
"""

p.display_legends = """
WS
WS TI=770ms
WS double inv TI=770ms
"""

p.run_pipeline_std()

# %%20/12/2019 - test_IR_sLASER - inv tests

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 8
p.data_filepaths = """
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0002.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0003.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0004.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0005.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0006.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0007.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0008.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0009.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0010.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0011.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0012.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0013.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0014.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0015.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0016.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0017.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0018.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0019.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0020.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0021.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0022.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0023.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0024.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0025.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0026.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0027.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0028.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0029.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0030.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0031.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0032.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0033.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0034.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0035.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0036.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0037.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0038.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0039.dcm
/home/tangir/crmbm/acq/test-ir-slaser/20191220/01_0012_slaser-r-n/original-primary_e09_0040.dcm
"""

p.data_ref_filepaths = """"""

p.display_legends = """
15
25
35
45
55
65
75
85
95
105
115
125
135
145
155
165
175
185
195
205
215
225
235
245
255
265
275
285
295
305
315
325
335
345
355
365
375
385
395
405
"""

p.phase_enable = False
p.analyse_snr_enable = False
p.analyse_linewidth_enable = False
p.calibrate_enable = False

p.run_pipeline_std()

# %%20/12/2019 - test_IR_sLASER - head coil - WS/noWS/FID modulus tests

get_ipython().magic('clear')

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID68_slaser_R_N=20+_1_longTE_SNR++++_FID48889.dat
/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID69_slaser_R_N=20+_1_longTE_SNR++++_FID48890.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID70_slaser_R_N=20+_1_longTE_SNR++++_FID48891.dat
/home/tangir/crmbm/acq_twix/test-ir-slaser/meas_MID70_slaser_R_N=20+_1_longTE_SNR++++_FID48891.dat
"""

p.display_legends = """
WS
noWS
"""

p.remove_water_enable = False
# p.remove_water_hsvd_components
p.fid_modulus = False
p.data_process_only_this_data_index = [0]

p.run_pipeline_std()

# %% 22/01/2020 - test_IR_sLASER - head coil - WS tests

get_ipython().magic('clear')
plt.close('all')

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = ["/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID72_slaser_R_N=20+_1_longTE_SNR++++_FID49978.dat",
                    "/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID49979.dat",
                    "/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID74_slaser_R_N=20+_1_longTE_SNR++++_FID49980.dat",
                    "/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID75_slaser_R_N=20+_1_longTE_SNR++++_FID49981.dat",
                    "/home/tangir/crmbm/acq_twix/test-ir-slaser2/meas_MID76_slaser_R_N=20+_1_longTE_SNR++++_FID49982.dat"]

p.data_ref_filepaths = []

p.display_legends = ["110",
                     "45",
                     "30",
                     "10",
                     "0"]

p.analyse_snr_s_range_ppm = [2.8, 3.2]
p.analyse_linewidth_range_ppm = [2.8, 3.2]

p.calibrate_POI_range_ppm = [4, 5]
p.calibrate_POI_true_ppm = 4.7

p.remove_water_enable = True
p.remove_water_hsvd_range = [4.5, 5]

# p.remove_water_hsvd_components
p.fid_modulus = True
p.data_process_only_this_data_index = []

p.run_pipeline_std()

# %% 03/03/2020 - Jojo's data about with several concentration of creatine

get_ipython().magic('clear')
plt.close('all')

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = ["/home/tangir/desktop/200228_phantoms/meas_MID30_slaser_[10]_WSAT_FID53582.dat",
                    "/home/tangir/desktop/200228_phantoms/meas_MID35_steam_[10]_WSAT_FID53587.dat",
                    "/home/tangir/desktop/200228_phantoms/meas_MID54_steam_[25]_WSAT_FID53606.dat",
                    "/home/tangir/desktop/200228_phantoms/meas_MID59_slaser_[25]_WSAT_FID53611.dat"]

p.data_ref_filepaths = ["/home/tangir/desktop/200228_phantoms/meas_MID31_slaser_[10]_NOWSAT_FID53583.dat",
                        "/home/tangir/desktop/200228_phantoms/meas_MID36_steam_[10]_NOWSAT_FID53588.dat",
                        "/home/tangir/desktop/200228_phantoms/meas_MID55_steam_[25]_NOWSAT_FID53607.dat",
                        "/home/tangir/desktop/200228_phantoms/meas_MID60_slaser_[25]_NOWSAT_FID53612.dat"]

p.display_legends = ["a",
                     "b",
                     "c",
                     "d"]

p.run_pipeline_std()

# %% 09/03/2020 - Comparing sLASER w/o offset exc

get_ipython().magic('clear')
plt.close('all')

p = reco.pipeline()
p.data_filepaths = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID100_slaser_R_N=20+_1_longTE_SNR++++_FID54192.dat",
                    "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID54194.dat",
                    "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID54194.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID101_slaser_R_N=20+_1_longTE_SNR++++_FID54193.dat",
                        "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID103_slaser_R_N=20+_1_longTE_SNR++++_FID54195.dat",
                        "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID104_slaser_R_N=20+_1_longTE_SNR++++_FID54196.dat"]

p.display_legends = """
offset = 0
offset = -1.7ppm
offset = -1.7ppm for WS, 0 for REF
"""

p.phase_enable = True
p.phase_order = 1
p.zerofill_npts = 4096 * 6
p.realign_enable = False
p.analyse_snr_evol_enable = False
p.analyse_linewidth_evol_enable = False
p.apodize_enable = True
p.apodize_npts = 4096 * 4
p.apodize_damping_hz = 1
p.calibrate_enable = False
p.display_offset = 0.0

s, s_ref = p.run_pipeline_std()

# %% 11/03/2020 - Jojo's creatine phantom : Comparing different concentrations at 7T

get_ipython().magic('clear')
plt.close('all')

p = reco.pipeline()
p.data_filepaths = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID30_slaser_[10]_WSAT_FID53582.dat",
                    "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID35_steam_[10]_WSAT_FID53587.dat",
                    "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID54_steam_[25]_WSAT_FID53606.dat",
                    "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID59_slaser_[25]_WSAT_FID53611.dat",
                    "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID75_steam_[50]_WSAT_FID53627.dat",
                    "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID80_slaser_[50]_WSAT_FID53632.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID31_slaser_[10]_NOWSAT_FID53583.dat",
                        "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID36_steam_[10]_NOWSAT_FID53588.dat",
                        "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID55_steam_[25]_NOWSAT_FID53607.dat",
                        "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID60_slaser_[25]_NOWSAT_FID53612.dat",
                        "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID81_slaser_[50]_NOWSAT_FID53633.dat",
                        "/home/tangir/crmbm/acq_twix/test01_phantom_10/meas_MID82_steam_[50]_NOWSAT_FID53634.dat"]

p.display_legends = """
sLASER on 10mM Cre
STEAM on 10mM Cre
STEAM on 25mM Cre
sLASER on 25mM Cre
STEAM on 50mM Cre
sLASER on 50mM Cre
"""

p.phase_order = 1
p.phase_offset = 0.0  # np.pi
p.recombine_phasing = False

p.realign_enable = False

p.zerofill_enable = False

p.apodize_enable = False

p.calibrate_POI_range_ppm = [2.5, 3.5]
p.calibrate_POI_true_ppm = 3.03

p.display_offset = 0.0
p.analyse_linewidth_enable = False
p.analyse_snr_s_range_ppm = [2.5, 3.5]
p.analyse_linewidth_range_ppm = [4.6, 6]

p.remove_water_enable = True
p.remove_water_hsvd_range = [4.5, 6]
p.remove_water_hsvd_components = 20

p.data_process_only_this_data_index = [2]
s, s_ref = p.run_pipeline_std()

# turn off realistic simulations because of non adiabatic pulses here
s.sequence.pulse_rfc_real_shape_enable = False

fit_metabolites = [xxx.m_Cr_CH3, xxx.m_Cr_CH2]
