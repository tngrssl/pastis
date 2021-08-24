## TODO

* enable log file, do not use root logger but give it a name "pastis"

* remove useless log messages using pause

* fix reco_templates issue

* clean-up git, remove datasets* files, notebooks and all unnecessary files from git, including their history (https://stackoverflow.com/questions/35115585/remove-files-completely-from-git-repository-along-with-its-history)

* leave 2 demo files: reco and fit examples

* complete readme with more information, code snippets, examples

* improve scaling job to make it automatic (normalization to a CONSTANT)

* check and update metabolite chemical shifts and J-coupling values with recent literature (important if going for 1H UHF brain MRS for example):

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

* extract VOI orientation from dcm/twix files and improve voxel_size, CSDE and VOI anatomy methods

* Parallelize stuff: data rejecting approaches

* Find a gitable non-binary format to store the metabolite db and basis sets (I was thinking of FODS but it is not really implemented in conda...)

* HDF5/ISMRMD format save/load: for now, only save the MRS signal with a few attributes. If needed, will complete implementation

* Implement partial volume stuff for absolute quantification
