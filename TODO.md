## TODO

- extract VOI orientation from dcm/twix files and improve voxel_size, CSDE and VOI anatomy methods

- weird slASER GAMMA implementation: impact for fit on not? check in vitro, brain, etc.

- Parallelize stuff: data rejecting approaches

- Find a gitable non-binary format to store the metabolite db and basis sets (I was thinking of FODS but it is not really implemented in conda...)

- HDF5/ISMRMD format save/load: for now, only save the MRS signal with a few attributes. If needed, will complete implementation

- Implement partial volume stuff for absolute quantification
