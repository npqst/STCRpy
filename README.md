# STCRpy

Author: Nele Quast


## Installation

The core STCRpy package can be installed from PyPi: 
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple STCRpy
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
conda install -c conda-forge -c schrodinger pymol-bundle numpy=1.26.0 -y
```


## Building and deploying
In case build and twine arent installed: 
```
pip install --upgrade build
pip install --upgrade twine
```

To build the package run: (this will build the wheel and compress the source code to `STCRpy/dist`)
```
python -m build
```

To deploy the package to the test-pypi server: 
```
python3 -m twine upload --repository testpypi dist/*
```
!! Note that to re-deploy the code the version number in `setup.py` must be updated so that pypi can stage the package. 


