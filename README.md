# stcrpy

Author: Nele Quast


## Installation

The core stcrpy package can be installed from PyPi: 
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple stcrpy
```

Certain extensions and functionality have additional dependencies: 

### Interaction Profiling
```
pip install --no-deps plip
```

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
To install a development version of stvrpy and enable local debugging: 
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple anarci-mhc
ANARCI --build_models
pip install -e .
```



