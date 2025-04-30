

<img src="./stcrpy_logo.png" alt="drawing" width="300"/>


# STCRpy 
[![stcrpy installation](https://github.com/npqst/STCRpy/actions/workflows/conda-workflow.yml/badge.svg)](https://github.com/npqst/STCRpy/actions/workflows/conda-workflow.yml)
[![stcrpy unittests](https://github.com/npqst/STCRpy/actions/workflows/unittest-workflow.yml/badge.svg)](https://github.com/npqst/STCRpy/actions/workflows/unittest-workflow.yml)
[![stcrpy_docs](https://readthedocs.org/projects/stcrpy/badge/?version=latest)](https://stcrpy.readthedocs.io/en/latest/)


Structural TCR python (STCRpy) is a software suite for analysing and processing T-cell receptor structures. 

Please feel free to reach out with any comments or feedback.

Under review, please cite: 

**Quast, N. , Deane, C., & Raybould, M. (2025). STCRpy: a software suite for TCR:pMHC structure parsing, interaction profiling, and machine learning dataset preparation. BioRxiv. https://doi.org/10.1101/2025.04.25.650667**

<img src="./stcrpy_main_fig.png" alt="drawing" width="1500"/>



# Installation

## TL;DR installation
```
pip install stcrpy
pip install plip
conda install -c conda-forge pymol-open-source  numpy -y
ANARCI --build_models           # this step will take a few minutes
```

## Step by step installation
We recommend installing STCRpy in a [conda](https://www.anaconda.com/docs/getting-started/miniconda/install#macos-linux-installation) (or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)) environment using python 3.9 to 3.12: 
```
conda create -n stcrpy_env python==3.12 -y
conda activate stcrpy_env
```

The core functionality of STCRpy can be installed as follows:
```
pip install stcrpy
```

After installing stcrpy, the anarci HMM models must be built to enable annotation.
```
ANARCI --build_models           # this step will take a few minutes
```

To enable interaction profiling, install PLIP (Adasme et. al., 2021):
```
pip install plip
```

To enable pymol visualisations, install pymol open source locally within the environment. Unfortunately, pymol currently needs to be installed even if you already have a pymol version. Be sure to install pymol within a managed conda (or mamba) environment to prevent interference with any existing versions. 
```
conda install -c conda-forge pymol-open-source  numpy -y
```

To generate pytorch and pytorch-geometric compatible datasets: 
```
pip install pytorch --index-url https://download.pytorch.org/whl/cpu
pip install torch_geometric
```
Note that this installs the CPU version of pytorch, for GPU / CUDA versions install according to the [pytorch installation docs](https://pytorch.org/get-started/locally/).

The EGNN example also uses `einops`. To install: 
```
pip install einops
```

# Documentation
STCRpy [documentation](https://stcrpy.readthedocs.io/en/latest/) is hosted on ReadtheDocs.

# Examples
STCRpy generates and operates on TCR structure objects. The majority of the API can be accessed through functions of the format: [`tcr.some_stcrpy_function()`](https://stcrpy.readthedocs.io/en/latest/stcrpy.tcr_processing.html#stcrpy.tcr_processing.TCR.TCR). TCR objects are associated with their MHC and antigen if these are presented in the structure. 

A notebook with examples can be found under [examples/STCRpy_examples.ipynb](./examples/STCRpy_examples.ipynb)

First import STCRpy:
```
import stcrpy
```

### To fetch a TCR structure from STCRDab or the PDB: 
```
tcr = stcrpy.fetch_TCR("8gvb")
```
This will return a TCR strcuture or object, or, if there are multiple copies of TCR crystal structures in the PDB file, will return a list containing TCR structure objects. It may be useful to unpack the list into distinct objects, or use python generators to operate on the lists. 

### To load a TCR structure from a PDB or MMCIF file:
```
tcr = stcrpy.load_TCR("filename.{pdb, cif}")
```

### To load multiple TCR structures from a list of files at once:
```
multiple_tcrs = stcrpy.load_TCRs([file_1, file_2, file_3])
```

### To save a TCR object to PDB or MMCIF files: 
```
tcr.save(filename.{pdb, cif})           # save the TCR and it's associated MHC and antigen
tcr.save(filename.{pdb, cif}, TCR_only=True)           # save the TCR only
```

### To calculate the TCR to pMHC geometry:
```
tcr.calculate_geometry()            # change the 'mode' keyword argument to change the geometry calculation method. See paper / documentation for details.
```

### To score the TCR to pMHC geometry:
```
tcr.score_docking_geometry()
```

### To profile interactions: 
```
tcr.profile_peptide_interactions()          # interaction profiling parameters can be adjusted, see documentation for details
```

### To visualise interactions:
```
tcr.visualise_interactions()
```

### To run full analysis on a set of TCR structures: 
```
from stcrpy.tcr_methods.tcr_batch_operations import analyse_tcrs
germlines_and_alleles_df, geometries_df, interactions_df = analyse_tcrs(list_or_dict_of_files)
```

### To generate graph datasets:
```
dataset = TCRGraphDataset(
            root=PATH_TO_DATASET,
            data_paths=PATH_TO_TCR_FILES
        )
```





