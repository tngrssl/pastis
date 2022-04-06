# PASTIS

## PASTIS

PASTIS stands for (P)rocessing (AS)sessment (T)echnique for (I)mproved (S)pectroscopy ;)

PASTIS can be used to process and quantify single-voxel Magnetic Resonance SPectroscopy (MRS) data. It can also simulate MRS data using various MRS sequences for different B0 fields. It was originally developed to reconstruct, process and quantify spinal cord MRS data at 7 T and has therefore special features related to motion detection and compensation. For these reasons, it was later on used for cardiac MRS but could work for any organ really :)

PASTIS relies a lot on the suspect package. More info here: https://github.com/openmrslab/suspect

## Credits

PASTIS was originally written by Tangi Roussel. If you are using PASTIS or part of it in your work, please cite the original paper:

Respiratory-triggered quantitative MR spectroscopy of the human cervical spinal cord at 7 T.
Roussel T, Le Fur Y, Guye M, Viout P, Ranjeva JP and Callot V. Magn Reson Med. 2022 https://doi.org/10.1002/mrm.29182

## Features

* Read raw data from Siemens (TWIX MR Syngo VB17 & VE11)
* Read raw data from Bruker (fid & rawdata files)
* Read dicom data (standard dicom, Siemens MR Syngo VB17, Siemens MR Syngo XA20)
* Read and write NIFTI MRS files
* VOI overlay on anatomical image
* Data processing
	* Automatic phasing
	* Automatic channel combination
	* Zero-filling
	* Fully automatic data discard for SNR and FWHM enhancement
	* Automatic frequency realignment
	* Apodization
	* Peak HSVD removal
	* Signal and noise estimation
	* Linewidth estimation
	* Spectral display
* Data simulation
	* Based on GAMMA library
	* PRESS, STEAM, sLASER sequences
	* Including real RF pulse shapes for sLASER
	* Fully editable metabolite basis set
	* Macromolecular baseline modelization
	* Linear combination time-domain model
	* Possible to save/load metabolite simulation signals
* Quantification
	* Based on the previous model
	* Dynamic model options
		* Parameters can be freezed
		* Parameters can be linked to each other
	* Jacobian matrix information
	* Cramér-Rao Bounds estimations
* Quantification using suspect's LCModel wrapper
* Dataframe-based storage in pkl files
* Reconstruction and quantification of diffusion-weighted MRS data (experimental)

## Installation

### Python package

PASTIS was developped on Python 3.8 and can be installed with *pip* using this command line:

```console
$ pip install -e path/to/the/folder/pastis
```

## Usage

In practice, PASTIS is accessible via the Python package called **mrs**.

### Data reconstruction

Let's say you acquired some MRS data on a MRI scanner and the raw data is stored in the following file:
```
./demo/data/spinal_cord_sLASER_shortTE_WS.dat
```

Let's start with package imports. The **reco** submodule is required for anything that has to do with data reconstruction.
```python
import pastis.reco as reco
```

PASTIS also uses a logger to output some information in the console (see the **log** submodule). You can choose the amount of debug messages to output with the setLevel (more info about Python loggers here: https://docs.python.org/3/howto/logging.html)
```python
import pastis.log as log
log.setLevel(log.INFO)
```

To read the raw data file, you write:
```python
my_mrs_signal = reco.MRSData2("./demo/data/spinal_cord_sLASER_longTE_WS.dat")
```

The object **my_mrs_signal** now contains the acquired data signal but also a lot of information regarding the acquisition sequence and the patient information. You can for example print the number of averages used during the acquisition:
```python
>>> print("NA = %d" % my_mrs_signal.sequence.na)
NA = 64
```

You can check the data signal dimensions:
```python
>>> print(my_mrs_signal.shape)
(64, 8, 8192)
```

And see that it is actually raw data: 64 averages acquired using a 8-channel Rx coil. Let's start by **phasing**:
```python
my_mrs_signal = my_mrs_signal.correct_phase_3d()
```

Now let's intensity **scale**, **combine the channels**, **average** all individual spectra and **crop** to 4096 time points. See below if you want more information on these methods.
```python
my_mrs_signal = my_mrs_signal.correct_intensity_scaling_nd()
my_mrs_signal = my_mrs_signal.correct_combine_channels_3d()
my_mrs_signal = my_mrs_signal.correct_average_2d()
my_mrs_signal = my_mrs_signal.correct_crop_1d(4096)
```

You can **display** the spectrum using:
```python
my_mrs_signal.display_spectrum_1d()
```

### Data quantification

The **fit** submodule is required for anything that has to do with data quantification. Some useful aliases are declared in the submodule **xxx**.
```python
import pastis.fit as fit
import pastis.aliases as xxx
```

Let's create a fitting tool for our dataset:
```python
my_fit_tool = fit.fit_pastis(my_mrs_signal)
```

Before running the fit, we need to configure the fitting tool. First, let's define which metabolite to fit.
```python
# before the fit, some peak area integration is performed: choose what peak to integrate
my_fit_tool.metabolites_area_integration = [xxx.m_NAA_CH3, xxx.m_Cr_CH3, xxx.m_Cho_CH3]
# and their respective integration width (ppm)
my_fit_tool.area_integration_peak_ranges = [0.1, 0.1, 0.1]

# choose which metabolites to include in the fit basis set
my_fit_tool.metabolites = [xxx.m_Water,
                        xxx.m_LipA,
                        xxx.m_LipB,
                        xxx.m_LipC,
                        xxx.m_NAA,
                        xxx.m_NAAG,
                        xxx.m_Cr_CH3,
                        xxx.m_Cr_CH2,
                        xxx.m_PCr,
                        xxx.m_GPC,
                        xxx.m_PC,
                        xxx.m_mI,
                        xxx.m_Glu,
                        xxx.m_Gln,
                        xxx.m_Tau]
```

Now, let's set the initial paramter values, the parametes ranges for the fit.
```python
# --- preparing minimum, maximum and initial fitting parameter sets ---

# create default minimum and maximum parameter sets
my_fit_tool.params_min = my_fit_tool.params_min.set_default_min().add_macromolecules_min()
my_fit_tool.params_max = my_fit_tool.params_max.set_default_max().add_macromolecules_max()

# initial parameter values for fit
my_fit_tool.params_init = (my_fit_tool.params_min + my_fit_tool.params_max) / 2.0

# --- concentrations

# fit ranges for concentration for all metabolites and macromolecules
my_fit_tool.params_min[:, xxx.p_cm] = 0.0
my_fit_tool.params_max[:, xxx.p_cm] = 200.0
# increase maximum concentration for all macromolecules
my_fit_tool.params_max[xxx.m_All_MMs, xxx.p_cm] = 1000.0

# set initial concentrations to zero
my_fit_tool.params_init[:, xxx.p_cm] = 0.0
# start fit with 0.1 concentrations for metabolites in basis set
my_fit_tool.params_init[my_fit_tool.metabolites, xxx.p_cm] = 0.1

# --- linewidths

# minimal damping for all metabolites and macromolecules
my_fit_tool.params_min[:, xxx.p_dd] = 5

# increase damping ranges for all macromolecules (these are usually broad peaks)
my_fit_tool.params_min[xxx.m_All_MMs, xxx.p_dd] = 150
my_fit_tool.params_max[xxx.m_All_MMs, xxx.p_dd] = 300

# start fit with minimal linewidth (= narrow peaks)
my_fit_tool.params_init[:, xxx.p_dd] = my_fit_tool.params_min[:, xxx.p_dd] * 1.1

# --- frequency shifts

# frequency shift ranges for all metabolites and MMs
my_fit_tool.params_min[:, xxx.p_df] = -10.0
my_fit_tool.params_max[:, xxx.p_df] = 10.0
# start fit with no frequency shifts
my_fit_tool.params_init[:, xxx.p_df] = 0.0

# -- phase shifts

# phase shift ranges for all metabolites and MMs
my_fit_tool.params_min[:, xxx.p_dp] = -0.1
my_fit_tool.params_max[:, xxx.p_dp] = +0.1
# start fit with no phase shifts
my_fit_tool.params_init[:, xxx.p_dp] = 0.0
```

And finally, let's set the constraints:
```python
# --- set relations between parameters ---

# first, let's lock all the metabolites
my_fit_tool.params_linklock[:] = xxx.ll_FIXED

# except the metabolites included in the fit basis set
my_fit_tool.params_linklock[my_fit_tool.metabolites, :] = xxx.ll_FREE

# we want a global phase for all metabolites and macromolecules
# in practice we will specify that all phase shifts are obeying the phase of Creatine
my_fit_tool.params_linklock[my_fit_tool.metabolites, xxx.p_dp] = xxx.ll_SLAVE1
my_fit_tool.params_linklock[xxx.m_Cr_CH3, xxx.p_dp] = xxx.ll_MASTER1

# --- an exception for water ---

# all water fitting parameters are free
my_fit_tool.params_linklock[xxx.m_Water, :] = xxx.ll_FREE

# water max concentration is increased
my_fit_tool.params_max[xxx.m_Water, xxx.p_cm] = 1000
```

We then run the fit and wait a few seconds. Some nice interactive plot will show up :)
```python
# run the fit
my_fit_tool.run()

```

You can print the final estimated parameters as following:
```python
>>> my_fit_tool.params_fit.print()
(INFO)  params.print: displaying parameters...
(INFO) ----------------------------------------
[#] Metabolite [cm]       [dd]       [df]       [dp]
# 0 Ala        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
# 1 Asc        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
# 2 Asp        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
# 3 Cho        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
# 4 Cho_CH3    ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
# 5 Cr_CH3     (39.1)⁰    (51.1)⁰    (-0.5)⁰    (-0.0)⁻²
# 6 Cr_CH2     (67.0)⁰    (100.0)⁰   (-6.9)⁰    (-0.0)⁺²
# 7 EA         ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
# 8 Eth        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
# 9 GABA       ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#10 Glc        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#11 Gln        (32.2)⁰    (33.2)⁰    (-1.7)⁰    (-0.0)⁺²
#12 Glu        ( 4.1)⁰    ( 5.3)⁰    (-0.1)⁰    (-0.0)⁺²
#13 Gly        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#14 GPC        ( 3.1)⁰    (14.5)⁰    (-0.8)⁰    (-0.0)⁺²
#15 Gsh        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#16 Lac        ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#17 mI         (77.2)⁰    (58.2)⁰    (-9.9)⁰    (-0.0)⁺²
#18 NAA        (28.5)⁰    (44.5)⁰    (-5.2)⁰    (-0.0)⁺²
#19 NAA_CH3    ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#20 NAAG       (28.3)⁰    (59.0)⁰    (-4.9)⁰    (-0.0)⁺²
#21 PCr        (18.1)⁰    (100.0)⁰   (-9.9)⁰    (-0.0)⁺²
#22 PC         (14.0)⁰    (72.4)⁰    (-7.6)⁰    (-0.0)⁺²
#23 PE         ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#24 sI         ( 0.0)⁺¹   ( 5.5)⁺¹   ( 0.0)⁺¹   ( 0.0)⁺¹
#25 Tau        ( 3.6)⁰    ( 6.0)⁰    (-1.4)⁰    (-0.0)⁺²
#26 Water      (709.5)⁰   (34.5)⁰    ( 3.9)⁰    (-0.1)⁰
(INFO) ----------------------------------------
```

### Pipeline

If you want to perform such tasks on a large amount of datasets, I suggest you to use the **reco.pipeline** class for the data reconstruction. An jupyter-notebook example is available in the **demo** folder.

## Data processing jobs
WIP!
