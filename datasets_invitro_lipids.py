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

get_ipython().magic("clear")
plt.close("all")

get_ipython().magic("matplotlib auto")
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.max_open_warning'] = 1000
plt.rcParams['font.size'] = 9
log.setLevel(log.INFO)

# template to use here (see definitions below)
reco_template = "lipids"

# %% "brain" reconstruction template
template_name = "lipids"

p = reco.pipeline()
p.settings["storage_file"] = "/home/tangir/crmbm/acq_db/%s.pkl" % template_name
p.settings["POI_range_ppm"] = [4.5, 5.2]
p.settings["POI_shift_range_ppm"] = [1, 2]
p.settings["POI_shift_true_ppm"] = 1.3
p.settings["POI_SNR_range_ppm"] = [1, 2]
p.settings["POI_LW_range_ppm"] = [4.5, 5.2]
p.settings["allowed_apodization"] = 1.0
p.settings["display_range_ppm"] = [0, 6]
p.settings["display"] = True

p.job_list = [  p.job["phasing"],
                p.job["scaling"],
                p.job["channel_combining"],
                p.job["noise_estimation"],
                p.job["averaging"],
                p.job["calibrating"],
                # p.job["water_removal"],
                # p.job["cropping"],
                p.job["apodizing"],
                p.job["displaying"]
                ]

p.analyze_enable = False
p.job["scaling"]["scaling_factor_dcm"] = 1 / 100
p.job["cropping"]["final_npts"] = 1024

p.save_template(template_name)

# %% 23/04/2021 steam / histo tests on foie gras sample
get_ipython().magic("clear")
plt.close("all")

p = reco.pipeline(reco_template)

# replaced
# /crmbm/data_cemerem/data
# by
# /run/user/12000/gvfs/smb-share:server=139.124.150.244,share=data_cemerem/data

p.dataset[0]["legend"] = "Histo T2 estimation 23/04/2021"
p.dataset[0]["dcm"]["files"] = ["/run/user/12000/gvfs/smb-share:server=139.124.150.244,share=data_cemerem/data/users/pdaude/vida/vat-test-1-mrs/20210423/9286_0027_histo-bh/original-primary-spectroscopy-none_e09_0001.dcm"]

p.dataset[1]["legend"] = "Histo T2 estimation 23/04/2021"
p.dataset[1]["dcm"]["files"] = ["/run/user/12000/gvfs/smb-share:server=139.124.150.244,share=data_cemerem/data/users/pdaude/vida/vat-test-1-mrs/20210423/9286_0027_histo-bh/original-primary-spectroscopy-none_e09_0002.dcm"]

p.dataset[2]["legend"] = "Histo T2 estimation 23/04/2021"
p.dataset[2]["dcm"]["files"] = ["/run/user/12000/gvfs/smb-share:server=139.124.150.244,share=data_cemerem/data/users/pdaude/vida/vat-test-1-mrs/20210423/9286_0027_histo-bh/original-primary-spectroscopy-none_e09_0003.dcm"]

p.dataset[3]["legend"] = "Histo T2 estimation 23/04/2021"
p.dataset[3]["dcm"]["files"] = ["/run/user/12000/gvfs/smb-share:server=139.124.150.244,share=data_cemerem/data/users/pdaude/vida/vat-test-1-mrs/20210423/9286_0027_histo-bh/original-primary-spectroscopy-none_e09_0004.dcm"]

p.dataset[4]["legend"] = "Histo T2 estimation 23/04/2021"
p.dataset[4]["dcm"]["files"] = ["/run/user/12000/gvfs/smb-share:server=139.124.150.244,share=data_cemerem/data/users/pdaude/vida/vat-test-1-mrs/20210423/9286_0027_histo-bh/original-primary-spectroscopy-none_e09_0005.dcm"]

p.run()
p.save_datasets()
