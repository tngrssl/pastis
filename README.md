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
	- Semi/fully automatic data rejection
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

- Use a gitable non-binary format for the metabolite db and basis sets (thinking of FODS)
- HDF5 format save/load: for now, only save the MRS signal with a few attributes. If needed, will complete implementation
- Implement partial volume stuff for absolute quantification
- Merge with Yasmin's code for DW MRS
- implement FQN fit criteria + real-time indicator
- peak-linewidth for prefit
- adapt reco/fit pipelines for heart/muscle
- extract header from dicoms: patient name, etc.
- keep track of data rejection stats in pipeline
+ check # TODO tags in the code :)