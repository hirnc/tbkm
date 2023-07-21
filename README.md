# tbkm
Programs used to analyze <i>Tricholoma bakamatsutake</b>, <i>T. matsutake</b>, and the related species.

# Prerequisite
- The HDF5 library
  Available from https://www.hdfgroup.org/solutions/hdf5/
  Tested with HDF5-1.14.1
- BamTools
  Available from https://github.com/pezmaster31/bamtools
  Tested with version 2.5.2

# Compiling
- Install all prerequisites. Modify Makefile if needed.
- `make' will produce two executables, gff_coverage and create_read_count_matrix, in the current directory
- Refer to the on-screen help (with -h option) for the details
