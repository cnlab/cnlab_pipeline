# cnlab_pipeline

### Python code for nipype pipelines used for standard lab analysis protocols

Each directory consists of a README.md, core scripts and parameter descriptions, and templates model specification files.

--------

* `BIDS` - Scripts to convert raw DICOMS to BIDS-compatible niftis
* `FC` - Functional connectivity code
* `FMRIPREP` - Scripts to run fMRIPrep
* `L1` - Python modules containing nipype workflow functions for first-level analysis
* `L2` - Python modules containing nipype workflow functions for second-level analysis 
* `MVPA` - Scripts to run MVPA analyses
* `PE` - Scripts to compute pattern expression values
* `PPI` - Python wrapper code to run PPI using SPM
* `ROI` - Scripts to extract parameter estimates from ROIs
* `archive` - archived scripts

Below is the file structure for the repository:
```
├── BIDS
│   ├── README.md
│   └── templates
├── FC
│   ├── README.md
│   ├── extract_signal.ipynb
│   ├── make_adjacency_matrix-FOR_PROJECT1.ipynb
│   ├── make_adjacency_matrix.ipynb
│   ├── make_adjacency_matrix_WITH_UPDATED_SIMILARITY_METRICS.ipynb
│   ├── make_adjacency_matrix_for_synchrony.ipynb
│   └── templates
├── FMRIPREP
│   ├── README.md
│   └── templates
├── L1
│   ├── README.md
│   ├── first_level.pyc
│   ├── l1analysis_SPM.py
│   ├── l1analysis_description.json
│   ├── l1analysis_notebook.ipynb
│   └── templates
├── L2
│   ├── README.md
│   ├── l2analysis_SPM.py
│   ├── l2analysis_description.json
│   ├── l2analysis_notebook.ipynb
│   └── templates
├── MVPA
│   ├── README.md
│   └── templates
├── PE
│   ├── README.md
│   └── templates
├── PPI
│   ├── BIDStest_ONR_gPPIv13_wrapper.ipynb
│   ├── BIDStest_ONR_gPPIv13_wrapper.m
│   ├── README.md
│   ├── batch8example_gPPIv13_wrapper.ipynb
│   └── templates
├── ROI
│   ├── README.md
│   └── templates
└── archive
```