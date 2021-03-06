## TODO

* ISMRMD format read/write: need to include some meta information related to MRS (dwell time, ppm0, sequence, etc.) but cannot figure it out...

* Make sure that the gamma library is not needed for pygamma to run. If confirmed, removed all code related to storing/loading of metabolite basis sets in PKL file?

* Improve data file handling by reading folders of dicom and twix files and putting everything in order. Problem is to find some kind of UID to link raw data (TWIX) to reconstructed data (DICOM). Some ReferenceImage# field used to do the job in VB17 version, not anymore >=VE11...

* Extract VOI orientation from dcm/twix files and improve voxel_size, CSDE and VOI anatomy methods

* Parallelize stuff for data rejecting approaches

* Find a gitable non-binary format to store the metabolite db and basis sets (I was thinking of FODS but it is not really implemented in conda...)

* Implement partial volume stuff for absolute quantification

* Check and update metabolite chemical shifts and J-coupling values with recent literature (important if going for 1H UHF brain MRS for example):

Article (Govindaraju2000)
Govindaraju, V.; Young, K. & Maudsley, A. A.
Proton NMR chemical shifts and coupling constants for brain metabolites
NMR Biomed, 2000, 13, 129-153

Article (Kreis2012)
Kreis, R. & Bolliger, C. S.
The need for updates of spin system parameters, illustrated for the case of γ-aminobutyric acid
NMR Biomed, 2012, 25, 1401-1403

Article (Govind2015)
Govind, V.; Young, K. & Maudsley, A. A.
Corrigendum: proton NMR chemical shifts and coupling constants for brain metabolites. Govindaraju V, Young K, Maudsley AA, NMR Biomed. 2000; 13: 129-153.
NMR in biomedicine, 2015, 28, 923-924

InBook (Govind2016)
Govind, V.
1H-NMR Chemical Shifts and Coupling Constants for Brain Metabolites
eMagRes, American Cancer Society, 2016, 1347-1362

Book (Graaf2019)
de Graaf, R. A.
In Vivo NMR Spectroscopy: Principles and Techniques
Wiley-Interscience, 2019
