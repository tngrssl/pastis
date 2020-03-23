#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script used to store calls for the reconstruction of 31P in vitro data.

@author: Tangi Roussel
"""

from __future__ import division
import matplotlib.pylab as plt
import mrs.reco as reco

from IPython import get_ipython
import warnings
warnings.filterwarnings("ignore", ".*GUI is implemented*")
plt.close('all')
get_ipython().magic('clear')

# %% 08/08/2019 - phosphates_mix_001 - levure/lessive/javel/serum

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq/phosphates_mix_001/phosphates-mix/20190808/01_0006_fid/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_001/phosphates-mix/20190808/01_0007_fid/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_001/phosphates-mix/20190808/01_0008_fid/original-primary_e09_0001.dcm
"""

p.display_legends = """
un peu de lessive
beaucoup
beaucoup beaucoup
"""

p.phase_enable = True
p.phase_POI_range_ppm = [-5, 5]

p.zerofill_enable = False
p.calibrate_enable = False
p.analyse_snr_enable = False
p.analyse_linewidth_enable = False
p.remove_water_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 60

p.display_magnitude_mode = True
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()


# %% 08/08/2019 - phosphates_mix_002 - levure/lessive/javel/serum

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0004_fid/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0002.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0003.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0004.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0005.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0006.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0007.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0008.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0009.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0006_svs-st-vapor-643-optim/original-primary_e09_0010.dcm
"""

p.display_legends = """
non-localized
FA test
FA test
FA test
FA test
FA test
FA test
FA test
FA test
FA test
FA test
"""

p.phase_enable = True
p.phase_POI_range_ppm = [-5, 5]

p.zerofill_enable = False
p.calibrate_enable = False
p.analyse_snr_enable = False
p.analyse_linewidth_enable = False
p.remove_water_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 60

p.display_magnitude_mode = True
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()


# %% 08/08/2019 - phosphates_mix_002 - levure/lessive/javel/serum, tests STEAM, sLASER, ISIS

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0010_slaser-r-n/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0012_bow-isis-15/original-primary_e09_0001.dcm
"""

p.display_legends = """
sLASER
ISIS
"""

p.phase_enable = True
p.phase_POI_range_ppm = [-5, 5]

p.zerofill_enable = False
p.calibrate_enable = False
p.analyse_snr_enable = False
p.analyse_linewidth_enable = False
p.remove_water_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 110

p.display_magnitude_mode = True
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()


# %% 08/08/2019 - phosphates_mix_002 - test raw data ISIS

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq/phosphates_mix_002/phosphates-mix/20190809/01_0012_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq_twix/phosphates_mix_002/meas_MID110_bow_isis_15_FID37704.dat
"""

p.data_ref_filepaths = """
"""

p.display_legends = """
dcm ISIS
raw ISIS
"""

p.display_amp_factor_list = [3, 1e8]

p.phase_enable = False
p.phase_display = False
p.phase_POI_range_ppm = [-5, 5]
p.phase_POI_range_ppm = [-2, +2]

# remove 1h channel from data
p.recombine_phasing = False
p.recombine_weights = [False, True]

p.realign_enable = False
p.zerofill_enable = False
p.calibrate_enable = False
p.remove_water_enable = False

p.analyse_snr_enable = True
p.analyse_snr_s_range_ppm = [-1, +1]
p.analyse_snr_n_range_ppm = [-20, -10]
p.analyse_snr_area_integrate = False
p.analyse_snr_magnitude_mode = True

p.analyse_linewidth_enable = True
p.analyse_linewidth_range_ppm = [-1, +1]
p.analyse_linewidth_magnitude_mode = True
p.analyse_linewidth_display = True

p.apodize_enable = True
p.apodize_damping_hz = 110

p.display_magnitude_mode = False
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()


# %% 13/08/2019 - phosphates_mix_003, optim exc pulse BW

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0016_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0017_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0018_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0019_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0020_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0021_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0022_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0023_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0024_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0025_bow-isis-15/original-primary_e09_0001.dcm
"""

p.data_ref_filepaths = """
"""

p.display_legends = """
1 kHz
2 kHz
3 kHz
4 kHz
5 kHz
6 kHz
7 kHz
8 kHz
9 kHz
10 kHz
"""

p.display_amp_factor_list = []

p.phase_enable = False
p.phase_display = False
p.phase_POI_range_ppm = [-5, 5]
p.phase_POI_range_ppm = [-2, +2]

# remove 1h channel from data
p.recombine_enable = True
p.recombine_phasing = False
p.recombine_weights = [False, True]

p.realign_enable = False
p.zerofill_enable = False
p.calibrate_enable = False
p.remove_water_enable = False

p.analyse_snr_enable = True
p.analyse_snr_s_range_ppm = [13, 15]
p.analyse_snr_n_range_ppm = [-20, -10]
p.analyse_snr_area_integrate = False
p.analyse_snr_magnitude_mode = True

p.analyse_linewidth_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 110

p.display_magnitude_mode = True
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()


# %% 13/08/2019 - phosphates_mix_003, optim exc pulse B1

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0026_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0027_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0028_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0029_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0030_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0031_bow-isis-15/original-primary_e09_0001.dcm
/home/tangir/crmbm/acq/phosphates_mix_003/phosphates-mix/20190813/01_0032_bow-isis-15/original-primary_e09_0001.dcm
"""

p.display_legends = """
50V
100V
150V
200V
250V
300V
320V
"""

p.display_amp_factor_list = []

p.phase_enable = False
p.phase_display = False
p.phase_POI_range_ppm = [-5, 5]
p.phase_POI_range_ppm = [-2, +2]

# remove 1h channel from data
p.recombine_enable = True
p.recombine_phasing = False
p.recombine_weights = [False, True]

p.realign_enable = False
p.zerofill_enable = False
p.calibrate_enable = False
p.remove_water_enable = False

p.analyse_snr_enable = True
p.analyse_snr_s_range_ppm = [13, 15]
p.analyse_snr_n_range_ppm = [-20, -10]
p.analyse_snr_area_integrate = False
p.analyse_snr_magnitude_mode = True

p.analyse_linewidth_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 110

p.display_magnitude_mode = True
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()


# %% 21/01/2020 - phantom_31p, optim B1 on levure tube

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID179_svs_st_vapor_643_optim_FID50085.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID180_svs_st_vapor_643_optim_FID50086.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID181_svs_st_vapor_643_optim_FID50087.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID182_svs_st_vapor_643_optim_FID50088.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID183_svs_st_vapor_643_optim_FID50089.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID184_svs_st_vapor_643_optim_FID50090.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID185_svs_st_vapor_643_optim_FID50091.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID186_svs_st_vapor_643_optim_FID50092.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID187_svs_st_vapor_643_optim_FID50093.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID188_svs_st_vapor_643_optim_FID50094.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID189_svs_st_vapor_643_optim_FID50095.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID191_svs_st_vapor_643_optim_FID50097.dat
"""

p.display_legends = """
300
275
250
225
200
175
150
125
100
75
50
25
"""

p.data_process_only_this_data_index = [2, 12]

p.display_amp_factor_list = []

p.phase_enable = True
p.phase_display = False
p.phase_POI_range_ppm = [-1, 1]

# remove 1h channel from data
p.recombine_enable = True
p.recombine_phasing = False
p.recombine_weights = [False, True]

p.realign_enable = False
p.zerofill_enable = False
p.calibrate_enable = False
p.remove_water_enable = False

p.analyse_snr_enable = True
p.analyse_snr_s_range_ppm = [-1, 1]
p.analyse_snr_n_range_ppm = [-20, -10]
p.analyse_snr_area_integrate = True
p.analyse_snr_magnitude_mode = False

p.analyse_linewidth_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 500

p.display_magnitude_mode = True
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()


# %% 21/01/2020 - phantom_31p, optim B1 on lessive

get_ipython().magic('clear')


plt.close("all")

p = reco.pipeline()
p.data_coil_nChannels = 1
p.ppm0 = 0.0
p.data_filepaths = """
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID209_svs_st_vapor_643_optim_FID50115.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID210_svs_st_vapor_643_optim_FID50116.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID211_svs_st_vapor_643_optim_FID50117.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID212_svs_st_vapor_643_optim_FID50118.dat
/home/tangir/crmbm/acq_twix/phantom_31p/meas_MID213_svs_st_vapor_643_optim_FID50119.dat
"""

p.display_legends = """
250
200
150
100
50
"""

p.data_process_only_this_data_index = []

p.display_amp_factor_list = []

p.phase_enable = True
p.phase_display = False
p.phase_POI_range_ppm = [-1, 1]

# remove 1h channel from data
p.recombine_enable = True
p.recombine_phasing = False
p.recombine_weights = [False, True]

p.realign_enable = False
p.zerofill_enable = False
p.calibrate_enable = False
p.remove_water_enable = False

p.analyse_snr_enable = True
p.analyse_snr_s_range_ppm = [-1, 1]
p.analyse_snr_n_range_ppm = [-20, -10]
p.analyse_snr_area_integrate = True
p.analyse_snr_magnitude_mode = False

p.analyse_linewidth_enable = False

p.apodize_enable = True
p.apodize_damping_hz = 100

p.display_magnitude_mode = True
p.display_range_ppm = [-25, 25]
s = p.run_pipeline_std()
