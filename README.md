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