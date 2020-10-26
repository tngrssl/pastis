# SHRIMP

## SHRIMP

SHRIMP stands for Spectroscopy Half-baked Resonance using Insane Magnetic Processing. That's a temporary acronym ;)

SHRIMP can be used to process and quantify single-voxel Magnetic Resonance SPectroscopy (MRS) data. It can also simulate MRS data using various MRS sequences for different B0 fields.

## Dependencies

- Python 3.7
- Python packages:
-- pyGAMMA (https://scion.duhs.duke.edu/vespa/gamma/wiki/PyGamma)
-- suspect
-- numpy
-- pandas
-- pymapvbvd

- GAMMA library installed on system (https://scion.duhs.duke.edu/vespa/gamma/)

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
- Dataframe-based final stats

## TODO

- slASER GAMMA implementation: impact for fit on not? check in vitro, etc.

- Parallelize stuff: data rejecting approaches

- Find a gitable non-binary format to store the metabolite db and basis sets (I was thinking of FODS but it is not really implemented in conda...)

- HDF5/ISMRMD format save/load: for now, only save the MRS signal with a few attributes. If needed, will complete implementation

- Implement partial volume stuff for absolute quantification

- Merge with Yasmin's code for DW MRS
