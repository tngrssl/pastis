# SHRIMP

## SHRIMP

SHRIMP stands for Spectroscopy Half-baked Resonance using Insane Magnetic Processing. That's a temporary acronym ;)

SHRIMP can be used to process and quantify single-voxel Magnetic Resonance SPectroscopy (MRS) data. It can also simulate MRS data using various MRS sequences for different B0 fields.

## Dependencies

- Python >= 3
- suspect
- pyGAMMA
- GAMMA library

## Features

- Read raw data from Siemens (TWIX format)
- Read dicom data
- Data processing
	- Automatic phasing
	- Automatic channel combination
	- Zero-filling
	- Fully automatic data rejection
	- Automatic frequency realignment
	- Apodization
	- Peak HSVD removal
	- Signal and noise estimation
	- Linewidth estimation
	- Spectral display
- Data simulation
	- Based on GAMMA library
	- PRESS, STEAM, sLASER sequences
	- Including real RF pulse shapes for sLASER
	- Fully editable metabolite basis set
	- Macromolecular baseline modelization
	- Linear combination time-domain model
- Quantification
	- Based on the previous model
	- Dynamic model options
		- Parameters can be freezed
		- Parameters can be linked to each other
	- Jacobian matrix information
	- Cramér-Rao Bounds estimations

## TODO

- Recode the way we input filepaths, legends etc. in dict way:
p.dataset_list["legend"]["raw-data"]["water-suppressed"] = "/..."
p.dataset_list["legend"]["raw-data"]["non-water-suppressed"] = "/..."
p.dataset_list["legend"]["dicoms"]["water-suppressed"] = "/..."
p.dataset_list["legend"]["dicoms"]["non-water-suppressed"] = "/..."
p.dataset_list["legend"]["physio"] = "/..."
p.dataset_list["legend"]["physio"] = "/..."

 + write code to convert all existing entries !

- better handling of POI ppm for processing: automatic and single
- data rejection fails on SC/P2, why ?
- find a way to save templates for reco pipeline (pkl?)
- Use a gitable non-binary format for the metabolite db and basis sets (thinking of FODS)
- HDF5/ISMRMD format save/load: for now, only save the MRS signal with a few attributes. If needed, will complete implementation
- Implement partial volume stuff for absolute quantification
- Merge with Yasmin's code for DW MRS
+ check # TODO tags in the code :)
