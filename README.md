[![DOI](https://zenodo.org/badge/885963051.svg)](https://doi.org/10.5281/zenodo.14161759)

# Differentiable Information Imbalance in Molecular Systems
Source code of data pipelines "[Automatic feature selection and weighting in molecular systems using Differentiable Information Imbalance](https://doi.org/10.48550/arXiv.2411.00851)".

## Installation

1) Ideally, create this `conda` environment to run these notebooks in. To install the environment 'dii_ms', which includes the version of DADApy of the publication:
```bash
conda env create -f environment.yml
```
If you are a `pip` user feel free to install the dependencies in `requirements.txt` into a VENV or anywhere else.
Alternatively, install the dependencies listed in the `environment.yml` using the environment manager of your choice.

2) Make a jupyter kernel from this environment for the subsequent analysis in jupyter notebooks:
```bash
python -m ipykernel install --user --name dii_ms
```

3) Analyze the jupyter notebooks. Data files necessary for calculations, and of precalculated results, were too large to store here and can be downloaded from OSF, fron the corresponding folders (see below).


## Source Data
All source data is located in an OSF Repository: https://osf.io/swtg5/.
Download it and add the contained data to the individual directories.

### CLN025 and Gaussian_random_variables_and_monomials
- Download all data files and folders from OSF and copy the missing large data files into the corresponding folder of this project.
- You can choose to use pre-calculated results for the plots or do the calculations. They may take several hours.
- Ignore warning messages about maxk. The DADApy project is under active development and you may chose to use newer versions instead of the here specified release 0.3.2.

### liquid_water_mlp_fitpredict
- Place the contents of the OSF repository in the directory `liquid_water_mlp_fitpredict/data`.
- There is already generated ACSF and SOAPs in `water_phase_store`.
    To generate new descriptors run `python make_descriptors.py`.
- Run `water_descriptors.ipynb` to remake the plot by loading the pregenerated results from `water_descriptors.ipynb`.
    Alternatively, you can perform feature selection on your own by executing the middle cells.
- `build_bp_inputfile.py` allows for automatic creation of input files from calculated weights.
    It is not cleanly written so use at your own risk.



