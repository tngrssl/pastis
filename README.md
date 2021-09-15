# PASTIS

## PASTIS

PASTIS stands for (P)rocessing (AS)sessment (T)echnique for (I)mproved (S)pectroscopy ;)

PASTIS can be used to process and quantify single-voxel Magnetic Resonance SPectroscopy (MRS) data. It can also simulate MRS data using various MRS sequences for different B0 fields.

PASTIS relies a lot on the suspect package. More info here: https://github.com/openmrslab/suspect

## Credits

PASTIS was originally written by Tangi Roussel. If you are using PASTIS or part of it in your work, please cite the original paper:

Respiratory-triggered quantitative MR spectroscopy of the human cervical spinal cord at 7 T.
Roussel T, Le Fur Y, Guye M, Viout P, Ranjeva JP and Callot V
(submitted)

## Features

* Read raw data from Siemens (TWIX MR Syngo VB17 & VE11)
* Read raw data from Bruker (fid files)
* Read dicom data (standard dicom, Siemens MR Syngo VB17, Siemens MR Syngo XA20)
* VOI overlay on anatomical image
* Data processing
	* Automatic phasing
	* Automatic channel combination
	* Zero-filling
	* Fully automatic data rejection for SNR and FWHM enhancement
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

PASTIS was developped on Python 3.7 and can be installed with *pip* using this command line:

```
pip install -e path/to/the/folder/pastis
```

### GAMMA

GAMMA is a system package used for NMR simulations. It can be installed on Debian-based system using this command-line:

```
apt-get install gamma
```

More info about GAMMA on https://scion.duhs.duke.edu/vespa/gamma/

### GAMMA alternative



## Usage

Coming soon...

### Data processing jobs
