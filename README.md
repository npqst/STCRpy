# STCRpy

Author: Nele Quast

## Pre-release alpha version (v0.1.24)
We anticipate changes to this repository. Please feel free to reach out with any comments or feedback.
Manuscript in preparation, please cite: 

**Quast, N. , Deane, C., & Raybould, M. (2025). STCRpy: a software suite for TCR:pMHC structure parsing, interaction profiling, and machine learning dataset preparation (0.1.24). Zenodo. https://doi.org/10.5281/zenodo.14962269**


## Installation

The core stcrpy package can be installed from PyPi. 
After installing stcrpy, the anarci HMM models must be built to enable annotation.
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple stcrpy
ANARCI --build_models
```

Certain extensions and functionality have additional dependencies: 

### PyTorch dataset generation
```
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install torch_geometric
```

### PyMOL visualisations 
```
conda install -c conda-forge -c schrodinger pymol-bundle numpy -y
```


## Building and deploying
In case build and twine arent installed: 
```
pip install --upgrade build
pip install --upgrade twine
```

To build the package run: (this will build the wheel and compress the source code to `stcrpy/dist`)
```
python -m build
```

To deploy the package to the test-pypi server: 
```
python3 -m twine upload --repository testpypi dist/*
```
!! Note that to re-deploy the code the version number in `setup.py` must be updated so that pypi can stage the package. 


## Debugging
To install a development version of stcrpy and enable local debugging: 
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple anarci-mhc
ANARCI --build_models
pip install -e .
```



