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

- Parallelize stuff: data rejecting approaches

- Find a gitable non-binary format to store the metabolite db and basis sets (I was thinking of FODS but it is not really implemented in conda...)

- HDF5/ISMRMD format save/load: for now, only save the MRS signal with a few attributes. If needed, will complete implementation

- Implement partial volume stuff for absolute quantification

- Merge with Yasmin's code for DW MRS

## Questions and things to try for the paper

- easy: check SNR threshold for which we need to but NA=2 or 4 moving average during processing

- slASER GAMMA implementation: impact for fit on not? check in vitro, etc.

- Process new data 349_AP_P2_MOELLE

- Show group-averaged rejected spectrum

- Evaluate data rejection approach with peak properties estimated before and after

- Compare theory and effective acquisition durations

- Check WS and noWS water FWHM

- Did 2nd order shim help? Check notebooks

- Water LW looks visually good in SC compared to brain right? So why is the data so different?

- Data rejection fails on SC/P2, why ?

- Find why MRS quality data varies from one subject to another (THIS IS CRITICAL)
-- Things to check:
    -- B0 (LW)
    -- B1 (Vref)
    -- BMI stuff (height/weight)
    -- unperfect breathing gating (check respiratory stability?)
    -- susceptiblity effects (check LW variation STD)
    -- lungs are closer to SC for small subjects (check voxel position, subject height)
- When doing this, normalize SNR per volume, NA, include NA rejected ? already done? voxel size, position ?

- SNR according to NA plot, normalized or not

- check for intra-subject repro if any data

- What does CRB errors depend on? scatter plots CRBS(met) = f(SNR, FWHM)
-- Use this to choose which metabolites to include, how many NA are really needed

- Ideally, implement 2 fit strategies:
1. NAA/Cho/Cre/mI like at 3T
2. add some other fun metabolites... GABA?
Depending on above

- CRBs: relative / absolute? work according to the last quantification review!

- Percentage of acceptable or quantifiable MR spectra?
