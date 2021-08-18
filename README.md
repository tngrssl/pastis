# PASTIS

## PASTIS

PASTIS stands for (P)rocessing (AS)sessment (T)echnique for (I)mproved (S)pectroscopy ;)

PASTIS can be used to process and quantify single-voxel Magnetic Resonance SPectroscopy (MRS) data. It can also simulate MRS data using various MRS sequences for different B0 fields.

PASTIS relies a lot on the suspect package. More info here: https://github.com/openmrslab/suspect

## Requirements (see requirements.txt)

* Python 3.7
* Python packages:
    * pygamma numpy suspect pandas pymapvbvd xlrd termcolor

* GAMMA library installed on system (https://scion.duhs.duke.edu/vespa/gamma/)
    * apt-get install gamma (Debian/Ubuntu)

## Features

* Read raw data from Siemens (TWIX MR Syngo VB17)
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
