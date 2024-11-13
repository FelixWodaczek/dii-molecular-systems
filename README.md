# Differentiable Information Imbalance in Molecular Systems
Source code of data pipelines "(Automatic feature selection and weighting in molecular systems using Differentiable Information Imbalance)[https://doi.org/10.48550/arXiv.2411.00851]".

## Installation

1) Ideally, create this `conda` environment to run these notebooks in. To install the environment 'dii_ms', which includes the version of DADApy of the publication:
```bash
conda env create -f environment.yml
```
Alternatively, install the dependencies listed in the `environment.yml` using the environment manager of your choice (eg. `pip`).

2) Make a jupyter kernel from this environment for the subsequent analysis in jupyter notebooks:
```bash
python -m ipykernel install --user --name dii_ms
```

3) Analyze the jupyter notebooks. Data files necessary for calculations, and of precalculated results, were too large to store here and can be downloaded from OSF, fron the corresponding folders (see below).


## Source Data
All source data is located in an (OSF Repository)[https://osf.io/swtg5/].
Download it and add the contained data to the individual directories.

## TODOs
- [ ] Remove old files from water_phases and check if imports work correctly
- [ ] in final iteration add compiled versions of all jupyter notebooks
