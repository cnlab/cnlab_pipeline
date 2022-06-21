This directory contains scripts to convert the DICOMS to BIDS-formatted niftis in bids_data.

`dicom2bids_heudiconv_bbprime.ipynb` converts the DICOMS to BIDS-formatted niftis and includes options to loop through participants or create job files that can be used on the slurm cluster (still in development).

`fieldmap_inclusion.ipynb` is used to append the `IntendedFor` parameter to the BIDS fmap JSON files to link fmaps and functional sequences.
