#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A user script to run shim and voltage prediction.

@author: Tangi Roussel
"""
# %% init
from IPython import get_ipython
import pickle
import mrs.predict as predict

p = predict.predict_shim()

p.data_filepaths = """
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID109_svs_st_vapor_643_optim_trig_FID29472.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID112_eja_svs_slaser_optim_trig_FID29475.dat
/home/tangir/crmbm/acq_twix/291-vs-moelle-spectro-p1/meas_MID145_eja_svs_slaser_optim_trig_FID29508.dat
/home/tangir/crmbm/acq_twix/296_ym_p1_brainmoelle/meas_MID155_slaser_R_N=20+_1_longTE_SNR++++_FID34189.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID105_slaser_R_N=20+_1_longTE_SNR++++_FID38626.dat
/home/tangir/crmbm/acq_twix/307-AP-P1-MOELLE/meas_MID124_slaser_R_N=10_2_longTE_SNR+++_FID38645.dat
/home/tangir/crmbm/acq_twix/308-rs-p1-moelle/meas_MID166_slaser_R_N=20+_1_longTE_SNR++++_FID38911.dat
/home/tangir/crmbm/acq_twix/310-mg-p1-moelle/meas_MID142_slaser_R_N=20+_1_longTE_SNR++++_FID39214.dat
/home/tangir/crmbm/acq_twix/311-sl-p1-moelle/meas_MID124_slaser_R_N=20+_1_longTE_SNR++++_FID39736.dat
/home/tangir/crmbm/acq_twix/313-ft-p1-moelle/meas_MID65_slaser_R_N=20+_1_longTE_SNR++++_FID41497.dat
/home/tangir/crmbm/acq_twix/314-yt-p1-moelle/meas_MID81_slaser_R_N=20+_1_longTE_SNR++++_FID41679.dat
/home/tangir/crmbm/acq_twix/316-ap-p1-moelle/meas_MID40_slaser_R_N=20+_1_longTE_SNR++++_FID42199.dat
/home/tangir/crmbm/acq_twix/319-fc-p1-moelle/meas_MID138_slaser_R_N=20+_1_longTE_SNR++++_FID43716.dat
/home/tangir/crmbm/acq_twix/328-af-p1-moelle/meas_MID66_slaser_R_N=20+_1_longTE_SNR++++_FID45773.dat
/home/tangir/crmbm/acq_twix/329-pi-p1-moelle/meas_MID170_slaser_R_N=20+_1_longTE_SNR++++_FID46234.dat
/home/tangir/crmbm/acq_twix/336-nb-p1-moelle/meas_MID75_slaser_R_N=20+_1_longTE_SNR++++_FID48206.dat
/home/tangir/crmbm/acq_twix/338-ro-p1-moelle/meas_MID114_slaser_R_N=20+_1_longTE_SNR++++_FID48496.dat
/home/tangir/crmbm/acq_twix/347-re-p1-moelle/meas_MID229_slaser_R_N=10_2_longTE_SNR+++_FID50588.dat
/home/tangir/crmbm/acq_twix/300-pm-p2-moelle/meas_MID69_slaser_R_N=20+_1_longTE_SNR++++_FID50927.dat
"""

p.build_library()
p.print_library_stats()
p.learn_from_library()
p.display()
p.save_model_to_disk("shim_vref_predidction_model.pkl")

# %% testing

get_ipython().magic('clear')

with open("shim_vref_predidction_model.pkl", 'rb') as f:
    [p] = pickle.load(f)

p.print_library_stats()
p.data_included_for_ref_voltage_prediction["Patient height"] = True
p.data_included_for_ref_voltage_prediction["Patient weight"] = True
p.data_included_for_ref_voltage_prediction["Patient age"] = True
p.data_included_for_ref_voltage_prediction["Patient sex"] = True
p.data_included_for_ref_voltage_prediction["VOI position"] = True

p.data_included_for_shim_prediction["Patient height"] = True
p.data_included_for_shim_prediction["Patient weight"] = True
p.data_included_for_shim_prediction["Patient age"] = False
p.data_included_for_shim_prediction["Patient sex"] = True
p.data_included_for_shim_prediction["Reference voltage"] = False
p.data_included_for_shim_prediction["VOI position"] = True
p._model_poly_order = 1

p.learn_from_library()
p.predict()
p.save_shims_to_disk("predicted_shims.txt")
