#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of in vivo data.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import matplotlib.pylab as plt
import mrs.reco as reco
import mrs.log as log
import numpy as np
get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
log.setLevel(log.DEBUG)

# %% 20/06/2019 - 296_ym_p1_brainmoelle - Yasmin :)
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0014_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0012_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID79_slaser_R_N=20+_1_longTE_SNR++++_FID33878.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID81_slaser_R_N=5_5+_shortTE_SNR++_FID33880.dat
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0008_steam-shortte-snr/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID71_steam_shortTE_SNR+_FID33870.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0013_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0009_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID80_slaser_R_N=5_5+_shortTE_SNR++_FID33879.dat
/home/tangir/crmbm/acq/296_ym_p1_brainmoelle/296-ym-p1-brainmoelle/20190619/01_0006_steam-shortte-snr/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID69_steam_shortTE_SNR+_FID33868.dat
"""

p.display_legends = """
brain - sLASER R:N=5:5 (DICOM)
brain - sLASER R:N=25:1 (DICOM)
brain - sLASER R:N=25:1 (TWIX)
brain - sLASER R:N=5:5 (TWIX)
brain - STEAM (DICOM)
brain - STEAM (TWIX)
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

p.data_process_only_this_data_index = [2, 3]
p.run()
p.save()

# %% 25/06/2019 - 296_ym_p1_brainmoelle - FID modulus tests
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 32
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID73_slaser_R_N=20+_1_longTE_SNR++++_FID33872.dat
"""

p.display_legends = """
brain - sLASER no VAPOR + FID modulus
"""

p.job_list = [  # p.jobs["phasing"],
                # p.jobs["scaling"],
                p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = False
p.run()
p.save()

p.display_legends = """
sLASER no VAPOR + conventionnal
"""

p.job_list = [  p.jobs["phasing"],
                p.jobs["scaling"],
                # p.jobs["FID modulus"],
                p.jobs["channel-combining"],
                # p.jobs["concatenate"],
                p.jobs["zero-filling"],
                # p.jobs["physio-analysis"],
                p.jobs["data-rejecting"],
                p.jobs["realigning"],
                p.jobs["averaging"],
                p.jobs["noise-estimation"],
                p.jobs["apodizing"],
                p.jobs["cropping"],
                p.jobs["water-removal"],
                p.jobs["calibrating"],
                p.jobs["displaying"]]

p.analyze_enable = True
p.run()
p.save()

# %% 27/08/2019 - 308-rs-p1-moelle - Ocha
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline()

p.data_filepaths = """
/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0024_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID210_slaser_R_N=20+_1_longTE_SNR++++_FID38955.dat
"""

p.data_ref_filepaths = """
/home/tangir/crmbm/acq/308-rs-p1-moelle/20190827/01_0025_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID211_slaser_R_N=20+_1_longTE_SNR++++_FID38956.dat
"""

p.display_legends = """
brain - sLASER 20:1 resp trig
brain - sLASER 20:1 resp trig (TWIX)
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

p.jobs["data-rejecting"]["auto"] = True
p.data_process_only_this_data_index = [1]
test=p.run()
p.save()
