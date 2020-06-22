#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vitro data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.aliases as xxx
import mrs.reco as reco
import mrs.log as log
import numpy as np

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.DEBUG)

rdb = reco.data_db()

# %% 27/02/2019 - spinal cord phantom - nice WS STEAM, testing TWIX reco
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False

p.data_process_only_this_data_index = [0]
p.jobs["scaling"]["scaling_factor"] = 57412 / 0.000692
p.run()

p.data_process_only_this_data_index = [1]
p.jobs["scaling"]["scaling_factor"] = 1
p.run()

# %% 27/02/2019 - spinal cord phantom - OVS optim
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - spoiler optim
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - various optim
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - TM optim
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 27/02/2019 - spinal cord phantom - compare svs-st-vapor-643, eja-svs-steam and eja_svs_slaser
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 12/03/2019 - spinal cord phantom - TR optim
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - VOI traces
get_ipython().magic("clear")
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

p.run()

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - sLASER
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM versus sLASER
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 20/03/2019 - fatty_braino_002 - sLASER optim invitro in invivo B1/B0 conditions - VOI 1/2 profiles
get_ipython().magic("clear")
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

p.run()

# %% 20/03/2019 - fatty_braino_002 - sLASER optim invitro in invivo B1/B0 conditions, R factor - VOI 1/2 profiles
get_ipython().magic("clear")
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

p.run()

# %%20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - testing TWIX reco (1/2)
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["scaling"]["scaling_factor"] = 1 / 0.00005 * 1 / 1.2
p.jobs["phasing"]["order"] = 0

p.analyze_enable = False
p.run()

# %%20/03/2019 - fatty_braino_002, bad B1/B0 conditions - STEAM - compare previous TWIX with DICOM (2/2)
get_ipython().magic("clear")

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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["scaling"]["scaling_factor"] = 1 / 6650
p.jobs["phasing"]["order"] = 0

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - sym vs asym STEAM
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                # p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - sym vs asym STEAM - VOI profiles
get_ipython().magic("clear")
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
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - settling time in STEAM
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                # p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - settling time in STEAM - VOI profiles
get_ipython().magic("clear")
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

p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - OVS in STEAM
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                # p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - OVS in STEAM - VOI profiles
get_ipython().magic("clear")
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

p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - ramp time in STEAM
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                # p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - OVS in STEAM - VOI profiles
get_ipython().magic("clear")
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

p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - spoilers in STEAM
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 03/05/2019 - fatty_braino_003, bad B1/B0 conditions - spoilers in STEAM - VOI profiles
get_ipython().magic("clear")
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

p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - VAPOR BW in STEAM
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - testing 30/30/60 OVS in STEAM - VOI profiles
get_ipython().magic("clear")
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

p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - R/pulse length in sLASER - VOI profiles
get_ipython().magic("clear")
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

p.run()

# %% 07/05/2019 - fatty_braino_004, bad B1/B0 conditions - R/N pulses in sLASER - VOI profiles
get_ipython().magic("clear")
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

p.analyze_selectivity_range_list = [[800, 3550], [-10600, -7800], [-3650, 1450]]

p.run()

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
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        # p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        # p.jobs["calibrating"]
                        ]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.run()

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
snr_results = []
for d in list(p._analyze_results_dict.keys()):
    last_job_key = list(p._analyze_results_dict[d]["snr"].keys())[-1]
    snr_results.append(p._analyze_results_dict[d]["snr"][last_job_key])

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
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        # p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        # p.jobs["calibrating"]
                        ]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.run()

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
snr_results = []
for d in list(p._analyze_results_dict.keys()):
    last_job_key = list(p._analyze_results_dict[d]["snr"].keys())[-1]
    snr_results.append(p._analyze_results_dict[d]["snr"][last_job_key])

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
get_ipython().magic("clear")
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
p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        # p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        # p.jobs["calibrating"]
                        ]

p.jobs["phasing"]["order"] = 0
p.analyze_enable = True
p.jobs["scaling"]["scaling_factor"] = 80000000
p.data_process_only_this_data_index = [0]
p.run()

p.jobs["scaling"]["scaling_factor"] = 1
p.data_process_only_this_data_index = [1]
p.run()

# %% 21/05/2019 - fatty_braino_007, bad B1/B0 conditions - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.analyze_enable = False
p.run()

# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.analyze_enable = False
p.run()

# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions, stronger spoilers - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.analyze_enable = False
p.run()

# %% 21/05/2019 - fatty_braino_007, better B1/B0 conditions, even more stronger spoilers - R/N pulses in sLASER TE minimal
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.analyze_enable = False
p.run()

# %% 03/06/2019 - fatty_braino_head_coil_002, test non-WS steam
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["phasing"]["order"] = 0
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.analyze_enable = False
p.run()

# %% 05/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - VOI profiles
get_ipython().magic("clear")
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

p.analyze_selectivity_range_list = [[-6000, -2000], [-9000, -5000], [-5000, -2000]]
p.run()

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
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]
                        ]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.display_offset = 30000
p.run()

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
snr_results = []
for d in list(p._analyze_results_dict.keys()):
    last_job_key = list(p._analyze_results_dict[d]["snr"].keys())[-1]
    snr_results.append(p._analyze_results_dict[d]["snr"][last_job_key])

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
get_ipython().magic("clear")
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]


p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        # p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]
                        ]

p.jobs["phasing"]["order"] = 0
p.jobs["analyzing-snr"]["s_range_ppm"] = [1.9, 2.1]
p.jobs["displaying"]["range_ppm"] = [1, 6]

p.analyze_enable = True
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing non WS STEAM spectra
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0012_steam-shortte-snr
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr
"""

p.display_legends = """
steam ws 90
steam ws 30deg
steam ws 1deg
"""

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False

p.jobs["analyzing-snr"]["s_range_ppm"] = [1.5, 2]
p.jobs["displaying"]["range_ppm"] = [1, 6]
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing DCM reco with right ref scans...
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0011_steam-shortte-snr
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/braino_alcoholic_002/braino-with-alcohol-belt/20190612/01_0013_steam-shortte-snr
"""

p.display_legends = """
steam dcm not phased
"""

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False

p.jobs["phasing"]["order"] = 0
p.jobs["displaying"]["range_ppm"] = [1, 6]
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing TWIX reco with right ref scans...
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False

p.jobs["apodizing"]["damping_hz"] = 20
p.jobs["scaling"]["scaling_factor"] = 1
p.data_process_only_this_data_index = [0]
p.run()

p.jobs["scaling"]["scaling_factor"] = 1e8 * 50800 / 26300
p.data_process_only_this_data_index = [1]
p.run()

# %% 12/06/2019 - fatty_braino_head_coil_002, bad B1/B0 conditions with head coil - testing TWIX reco and non WS STEAM
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID91_steam_shortTE_SNR+_FID33707.dat
/home/tangir/crmbm/acq_twix/braino_alcoholic_002/meas_MID92_steam_shortTE_SNR+_FID33708.dat
"""

p.display_legends = """
VAPOR FA=30deg
VAPOR FA=1deg
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.jobs["phasing"]["POI_range_ppm"] = [4, 5]
p.jobs["phasing"]["order"] = 1
p.jobs["apodizing"]["damping_hz"] = 20

p.analyze_enable = False
p.run()

# %%16/05/2019 - fatty_braino_007, bad B1/B0 conditions - STEAM - testing TWIX reco
get_ipython().magic("clear")
plt.close("all")

for chan_to_turn_off in range(0, 9):

    print("chan_to_turn_off=" + str(chan_to_turn_off))
    p = reco.pipeline()
    p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID31_eja_svs_slaser_optim_FID32220.dat"]

    p.data_filepaths = ["/home/tangir/crmbm/acq_twix/fatty_braino_007/meas_MID30_eja_svs_slaser_optim_FID32219.dat"]

    p.display_legends = ["off channel = " + str(chan_to_turn_off)]

    p.job_list = [  p.jobs["phasing"],
                    p.jobs["scaling"],
                    # p.jobs["FID modulus"],
                    p.jobs["channel-combining"],
                    # p.jobs["concatenate"],
                    # p.jobs["zero-filling"],
                    # p.jobs["physio-analysis"],
                    # p.jobs["data-rejecting"],
                    p.jobs["realigning"],
                    p.jobs["averaging"],
                    p.jobs["noise-estimation"],
                    p.jobs["apodizing"],
                    # p.jobs["cropping"],
                    # p.jobs["water-removal"],
                    p.jobs["calibrating"],
                    p.jobs["displaying"]]

    p.jobs["phasing"]["order"] = 1
    p.jobs["channel-combining"]["phasing"] = False

    combine_weights = np.full([8, ], True)
    if(chan_to_turn_off < 8):
        combine_weights[chan_to_turn_off] = False
    p.jobs["channel-combining"]["weights"] = combine_weights

    p.run()

# %%27/11/2019 - phantom_test_inv, sLASER - inv test
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %%20/12/2019 - test_IR_sLASER - inv tests
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
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

p.job_list = [  # p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()

# %%20/12/2019 - test_IR_sLASER - head coil - WS/noWS/FID modulus tests
plt.close('all')

p = reco.pipeline()
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.data_process_only_this_data_index = [0]
p.run()

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.data_process_only_this_data_index = [1]
p.run()

# %% 22/01/2020 - test_IR_sLASER - head coil - WS tests
plt.close('all')

p = reco.pipeline()
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]
                        ]

p.jobs["analyzing-snr"]["s_range_ppm"] = [2.8, 3.2]
p.jobs["analyzing-lw"]["range_ppm"] = [2.8, 3.2]
p.jobs["calibrating"]["POI_range_ppm"] = [4, 5]
p.jobs["calibrating"]["POI_true_ppm"] = 4.7

p.data_process_only_this_data_index = []
p.run()

# %% 09/03/2020 - Comparing sLASER w/o offset exc
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                # p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                # p.jobs["cropping"],
                # p.jobs["water-removal"],
                # p.jobs["calibrating"],
                p.jobs["displaying"]]


p.jobs["phasing"]["order"] = 1
p.jobs["zero-filling"]["npts"] = 4096 * 6
p.jobs["apodizing"]["damping_hz"] = 1

p.display_offset = 0.0
p.analyze_enable = False
p.run()

# %% 11/03/2020 - Jojo's creatine phantom : Comparing different concentrations at 7T
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

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]
                        ]

p.jobs["phasing"]["order"] = 1
p.jobs["phasing"]["offset"] = 0.0  # np.pi
p.jobs["channel-combining"]["phasing"] = False

p.jobs["calibrating"]["POI_range_ppm"] = [2.5, 3.5]
p.jobs["calibrating"]["POI_true_ppm"] = 3.03

p.jobs["analyzing-snr"]["s_range_ppm"] = [2.5, 3.5]
p.jobs["analyzing-lw"]["range_ppm"] = [2.5, 3.5]

p.data_process_only_this_data_index = [2]
data_list = p.run()

# %% 26/05/2020 - Ethanol tests
plt.close('all')

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID95_steam_shortTE_SNR+_FID54187.dat",
                    "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID102_slaser_R_N=20+_1_longTE_SNR++++_FID54194.dat",
                    "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID105_slaser_R_N=20+_1_longTE_SNR++++_FID54197.dat",
                    "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID107_slaser_R_N=20+_1_longTE_SNR++++_FID54199.dat",
                    "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID109_slaser_R_N=20+_1_longTE_SNR++++_FID54201.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID96_steam_shortTE_SNR+_FID54188.dat",
                        "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID103_slaser_R_N=20+_1_longTE_SNR++++_FID54195.dat",
                        "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID106_slaser_R_N=20+_1_longTE_SNR++++_FID54198.dat",
                        "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID108_slaser_R_N=20+_1_longTE_SNR++++_FID54200.dat",
                        "/home/tangir/crmbm/acq_twix/test-slaser-ethanol/meas_MID110_slaser_R_N=20+_1_longTE_SNR++++_FID54202.dat"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["phasing (suspect)"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]
                        ]

p.display_legends = """
Ethanol - STEAM
Ethanol - sLASER 50/1 TE=50ms
Ethanol - sLASER 50/1 TE=90ms
Ethanol - sLASER 50/1 TE=130ms
Ethanol - sLASER 20/1 TE=130ms
"""

p.jobs["phasing"]["order"] = 1
p.jobs["zero-filling"]["npts"] = 16384 + 1024
p.jobs["cropping"]["final_npts"] = 16384 + 1024
p.jobs["apodizing"]["damping_hz"] = 15
p.jobs["calibrating"]["POI_range_ppm"] = [4.5, 5]
p.jobs["calibrating"]["POI_true_ppm"] = 4.7

p.analyze_enable = False
p.data_process_only_this_data_index = [3, 4]
p.run()
p.save(rdb)


# %% 03/06/2020 - Ethanol night tests
plt.close('all')

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID200_slaser_R_N=20+_1_longTE_SNR++++_FID56493.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID203_slaser_R_N=20+_1_longTE_SNR++++_FID56496.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID205_slaser_R_N=20+_1_longTE_SNR++++_FID56498.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID208_slaser_R_N=20+_1_longTE_SNR++++_FID56501.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID210_slaser_R_N=20+_1_longTE_SNR++++_FID56503.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID212_slaser_R_N=20+_1_longTE_SNR++++_FID56505.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID218_svs_se_30_FID56511.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID221_svs_se_30_FID56514.dat",
                    "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID228_eja_svs_press_FID56521.dat",
                    "/home/tangir/crmbm/acq_twix/glu_slaser_head/meas_MID263_slaser_R_N=20+_1_longTE_SNR++++_FID56556.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID204_slaser_R_N=20+_1_longTE_SNR++++_FID56867.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID206_slaser_R_N=20+_1_longTE_SNR++++_FID56869.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID208_slaser_R_N=20+_1_longTE_SNR++++_FID56871.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID226_slaser_R_N=20+_1_longTE_SNR++++_FID56889.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID228_slaser_R_N=20+_1_longTE_SNR++++_FID56891.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID230_slaser_R_N=20+_1_longTE_SNR++++_FID56893.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID252_slaser_R_N=10_2_longTE_SNR+++_FID56915.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID255_slaser_R_N=20+_1_longTE_SNR++++_FID56918.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID257_slaser_R_N=20+_1_longTE_SNR++++_FID56920.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID259_slaser_R_N=20+_1_longTE_SNR++++_FID56922.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID261_slaser_R_N=20+_1_longTE_SNR++++_FID56924.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID263_slaser_R_N=20+_1_longTE_SNR++++_FID56926.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID276_slaser_R_N=10_2_longTE_SNR+++_FID56939.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID280_slaser_R_N=20+_1_longTE_SNR++++_FID56943.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID283_slaser_R_N=20+_1_longTE_SNR++++_FID56946.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID285_slaser_R_N=20+_1_longTE_SNR++++_FID56948.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID287_slaser_R_N=20+_1_longTE_SNR++++_FID56950.dat",
                    "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID289_slaser_R_N=20+_1_longTE_SNR++++_FID56952.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID201_slaser_R_N=20+_1_longTE_SNR++++_FID56494.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID204_slaser_R_N=20+_1_longTE_SNR++++_FID56497.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID206_slaser_R_N=20+_1_longTE_SNR++++_FID56499.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID209_slaser_R_N=20+_1_longTE_SNR++++_FID56502.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID211_slaser_R_N=20+_1_longTE_SNR++++_FID56504.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID213_slaser_R_N=20+_1_longTE_SNR++++_FID56506.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID219_svs_se_30_FID56512.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID222_svs_se_30_FID56515.dat",
                        "/home/tangir/crmbm/acq_twix/ethanol_slaser_head/meas_MID230_eja_svs_press_FID56523.dat",
                        "/home/tangir/crmbm/acq_twix/glu_slaser_head/meas_MID264_slaser_R_N=20+_1_longTE_SNR++++_FID56557.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID205_slaser_R_N=20+_1_longTE_SNR++++_FID56868.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID207_slaser_R_N=20+_1_longTE_SNR++++_FID56870.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID209_slaser_R_N=20+_1_longTE_SNR++++_FID56872.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID227_slaser_R_N=20+_1_longTE_SNR++++_FID56890.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID229_slaser_R_N=20+_1_longTE_SNR++++_FID56892.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID231_slaser_R_N=20+_1_longTE_SNR++++_FID56894.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID253_slaser_R_N=10_2_longTE_SNR+++_FID56916.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID256_slaser_R_N=20+_1_longTE_SNR++++_FID56919.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID258_slaser_R_N=20+_1_longTE_SNR++++_FID56921.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID260_slaser_R_N=20+_1_longTE_SNR++++_FID56923.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID262_slaser_R_N=20+_1_longTE_SNR++++_FID56925.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID264_slaser_R_N=20+_1_longTE_SNR++++_FID56927.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID277_slaser_R_N=10_2_longTE_SNR+++_FID56940.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID282_slaser_R_N=20+_1_longTE_SNR++++_FID56945.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID284_slaser_R_N=20+_1_longTE_SNR++++_FID56947.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID286_slaser_R_N=20+_1_longTE_SNR++++_FID56949.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID288_slaser_R_N=20+_1_longTE_SNR++++_FID56951.dat",
                        "/home/tangir/crmbm/acq_twix/eth-gly-slaser/meas_MID290_slaser_R_N=20+_1_longTE_SNR++++_FID56953.dat"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                # p.jobs["phasing (suspect)"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]
                        ]

p.display_legends = """
Ethanol2 - 20/1 R=70
Ethanol2 - 20/1 R=20
Ethanol2 - 20/1 R=70 5ms
Ethanol2 - 20/1 TE=100 R=70
Ethanol2 - 20/1 TE=100 R=20
Ethanol2 - 20/1 TE=100 R=70 5ms
Ethanol2 - PRESS TE=50ms
Ethanol2 - PRESS TE=30ms
Ethanol2 - ejaPRESS TE=30ms
Glu - 20/1 R=70
Glycerol - 20/1 R=70
Glycerol - 20/1 R=20
Glycerol - 20/1 R=70 5ms
Ethanol - 20/1 R=70
Ethanol - 20/1 R=20
Ethanol - 20/1 R=70 5ms
Eth3 - 10/2
Eth3 - 20/1 6ms
Eth3 - 20/1 7ms
Eth3 - 20/1 8ms
Eth3 - 20/1 9ms
Eth3 - 20/1 10ms
Gly3 - 10/2
Gly3 - 20/1 6ms
Gly3 - 20/1 7ms
Gly3 - 20/1 8ms
Gly3 - 20/1 9ms
Gly3 - 20/1 10ms
"""

p.jobs["phasing"]["display"] = False
p.jobs["phasing"]["order"] = 0
p.jobs["channel-combining"]["phasing"] = False
p.jobs["zero-filling"]["npts"] = 16384 + 1024
p.jobs["cropping"]["final_npts"] = 16384 + 1024
p.jobs["apodizing"]["damping_hz"] = 5
p.jobs["calibrating"]["POI_range_ppm"] = [4.5, 5]
p.jobs["calibrating"]["POI_true_ppm"] = 4.7
p.jobs["displaying"]["range_ppm"] = [0, 5]

p.analyze_enable = False
p.data_process_only_this_data_index = np.arange(16,28)
p.run()
p.save(rdb)


# %% 03/06/2020 - Creatine tubes in water/agar
plt.close('all')

p = reco.pipeline()

p.data_filepaths = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID107_slaser_R_N=20+_1_longTE_SNR++++_FID58514.dat",
                    "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID109_slaser_R_N=20+_1_longTE_SNR++++_FID58516.dat",
                    "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID111_slaser_R_N=20+_1_longTE_SNR++++_FID58518.dat",
                    "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID248_slaser_R_N=20+_1_longTE_SNR++++_FID58655.dat",
                    "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID252_slaser_R_N=20+_1_longTE_SNR++++_FID58659.dat",
                    "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID254_slaser_R_N=20+_1_longTE_SNR++++_FID58661.dat"]

p.data_ref_filepaths = ["/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID108_slaser_R_N=20+_1_longTE_SNR++++_FID58515.dat",
                        "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID110_slaser_R_N=20+_1_longTE_SNR++++_FID58517.dat",
                        "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID112_slaser_R_N=20+_1_longTE_SNR++++_FID58519.dat",
                        "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID251_slaser_R_N=20+_1_longTE_SNR++++_FID58658.dat",
                        "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID253_slaser_R_N=20+_1_longTE_SNR++++_FID58660.dat",
                        "/home/tangir/crmbm/acq_twix/q-mrs-tests-cre/meas_MID255_slaser_R_N=20+_1_longTE_SNR++++_FID58662.dat"]

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                # p.jobs["data-rejecting"],
                # p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                # p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["phasing (suspect)"],
                p.jobs["displaying"]]

p.analyze_job_list = [  p.jobs["channel-combining"],
                        p.jobs["zero-filling"],
                        p.jobs["realigning"],
                        p.jobs["averaging"],
                        p.jobs["calibrating"]
                        ]

p.display_legends = """
Cre - 25M - TE=50ms
Cre - 25M - TE=90ms
Cre - 25M - TE=130ms
Cre - 10M - TE=50ms
Cre - 10M - TE=90ms
Cre - 10M - TE=130ms
"""

p.jobs["zero-filling"]["npts"] = 16384 + 1024
p.jobs["cropping"]["final_npts"] = 16384 + 1024
p.jobs["apodizing"]["damping_hz"] = 5
p.jobs["calibrating"]["POI_range_ppm"] = [4.5, 5]
p.jobs["calibrating"]["POI_true_ppm"] = 4.7
p.jobs["displaying"]["range_ppm"] = [0, 5]

p.analyze_enable = False
p.data_process_only_this_data_index = []
p.run()
p.save(rdb)
