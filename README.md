

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
conda install pymol-open-source -y
ANARCI --build_models           # this step will take a few minutes
```

## Step by step installation
We recommend installing STCRpy in a [conda](https://www.anaconda.com/docs/getting-started/miniconda/install#macos-linux-installation) (or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)) environment using python 3.9 to 3.12. You can also use a python virtual environment if you do not need pymol visualisations. 

<details> <summary>conda</summary>

```
conda create -n stcrpy_env python==3.12 -y
conda activate stcrpy_env
```

</details>
<details> <summary>mamba</summary>

```
mamba create -n stcrpy_env python==3.12 -y
mamba activate stcrpy_env
```

</details>
<details> <summary>venv</summary>

```
python -m venv stcrpy_env
source stcrpy_env/bin/activate
```

</details>


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
conda install pymol-open-source -y
```

To generate pytorch and pytorch-geometric compatible datasets (see the [pytorch docs](https://pytorch.org/get-started/locally/) for hardware specific instructions): 
```
pip install stcrpy[ml_datasets]
```

> Note that the installs for pytorch can be platform specific.
> If errors are ecountered here it is best to manually install the depedencies following the [pytorch installation docs](https://pytorch.org/get-started/locally/).
> For example:
> ```
> pip install torch --index-url https://download.pytorch.org/whl/cpu
> pip install torch_geometric
> ```
> This installs the CPU version of pytorch (for GPU / CUDA versions follow the install [pytorch installation docs](https://pytorch.org/get-started/locally/)).
>
> The EGNN example also uses `einops`. Which can be manually installed as follows:
> ```
> pip install einops
> ```

### Docker container installation

This directory contains [Docker](https://www.docker.com/) configuration for running STCRpy in a containerized environment with uv dependency management. To build the docker image run:

```
docker build -t stcrpy .
```

And then run and example with:

```
mkdir output
docker run --rm \
  -v $(pwd)/output:/app/output \
  stcrpy python -c "
import stcrpy
tcr = stcrpy.fetch_TCRs('8gvb')[0]
tcr.profile_peptide_interactions()
tcr.get_interaction_heatmap(plotting_kwargs={'save_as': '/app/output/heatmap.png'})
"
```

If all goes well, there should be a heatmap.png file in the "./output" folder.

# Documentation
STCRpy [documentation](https://stcrpy.readthedocs.io/en/latest/) is hosted on ReadtheDocs.

# Examples
STCRpy generates and operates on TCR structure objects. The majority of the API can be accessed through functions of the format: `tcr.some_stcrpy_function()`. ([See TCR object docs here](https://stcrpy.readthedocs.io/en/latest/stcrpy.tcr_processing.html#stcrpy.tcr_processing.TCR.TCR)). TCR objects are associated with their MHC and antigen if these are presented in the structure. 

A notebook with examples can be found under [examples/STCRpy_examples.ipynb](./examples/STCRpy_examples.ipynb)

First import STCRpy:
```
import stcrpy
```

### To fetch a TCR structure from STCRDab or the PDB: 
```
multiple_tcrs = stcrpy.fetch_TCRs("8gvb")
```
This will return a list of all of the TCR structures found in the PDB file, represented as TCR structure objects.

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

### To calculate TCR prediction metrics such as RMSD, interface RMSD (of the TCR:pMHC interface) or DockQ scores:

```
# RMSD
from stcrpy.tcr_metrics import RMSD

rmsd_calculator = RMSD()
rmsd = rmsd_calculator.calculate_rmsd(pred_tcr, reference_tcr, save_alignment=False)      # Calculates the RMSD of each region of the TCR. To check the alignment set save_alignment to True.

# To calculate RMSD for a set of predictions against a set of reference structures from files: 
files = list(zip(prediction_files, reference_files))
rmsd_df = rmsd_calculator.rmsd_from_files(files)



# Interface RMSD of TCR:pMHC interface
from stcrpy.tcr_metrics import InterfaceRMSD

interface_rmsd_calculator = InterfaceRMSD()
irmsds = interface_rmsd_calculator.get_interface_rmsd(tcr, reference_tcr)

# DockQ
from stcrpy.tcr_metrics.tcr_dockq import TCRDockQ

dockq_calculator = TCRDockQ()               # by default this will merge the TCR and pMHC chains and calculate DockQ of the complete TCR:pMHC interface. To calculate DockQ scores per chain, use TCR_pMHC_interface=False
dockq_results = dockq_calculator.tcr_dockq(tcr, reference_tcr, save_merged_complex=False)           # to investigate the merged TCR:pMHC structure set save_merged_complex=True 

```

### Torsion angles and internal coordinates
STCRpy builds upon the Biopython PDB module, and you can calculate the internal coordinates, such as backbone torsion angles, using the [`internal_coordinates` function](https://biopython.org/docs/dev/api/Bio.PDB.internal_coords.html).

```
# internal coordinate calculations should be made per chain
for c in tcr.get_chains():
    c.atom_to_internal_coordinates()            # calculate the internal coordinates 

# internal coordinates can be accessed per residue: 
res = next(tcr.get_residues())
res.internal_coord.get_angle("psi")         # retrieve angles via angle keys
```

### Domain angles between TCR chains
STCRpy can be used to calculate the geometry and angles between the TCR variable domains of abTCRs and gdTCRs. This follows the ABangle implementation [(Dunbar et al. 2013)](https://academic.oup.com/peds/article/26/10/611/1509255).
```
tcr.get_TCR_angles()

# returns dictionary of TCR domain angles and measurements.
# For example: 
# {
#    'BA': -56.72234454750631,
#    'BC1': 122.55277240895967,
#    'AC1': 73.96532018128327,
#    'BC2': 82.63524566165464,
#    'AC2': 99.60327202896609,
#    'dc': np.float64(15.606353954437227)
# }

```


# Symmetry mate handling
Some TCR:pMHC crystals are formed of repeating cell units in which the TCR and the antigen do not directly contact. 
STCRpy generates symmetry mates in these cases to pair pMHC with TCRs in the structure. 
Note that symmetry mate generation requires pymol to be installed. By default, symmetry mate generation is enabled, however, it can be toggled by setting: 
`include_symmetry_mates=False` in `get_tcr_structure`. 

## Example: 
```
tcr_6ulr_paired_antigen = stcrpy.fetch_TCRs("6ulr")
tcr_6ulr_no_antigen = stcrpy.fetch_TCRs("6ulr", include_symmetry_mates=False)       # does not generate symmetry mates

```



# Limitations

## Connected peptide chains
STCRpy is currently not configured to handle cases where the antigen peptide is connected to the TCR or MHC chain - this is primarily because the parsing pipeline operates on chain objects and it can be tricky to consistently separate the peptide segment from the remainder of the TCR chain. A known case is PDB code 6MNO. 

## Gamma-Delta TCR geometry
STCRpy supports gamma-delta TCR parsing, interaction profiling and visusalisation, but is not currently configured to calculate gd-TCR geometry.

## MHC Class II geometry scoring
STCRpy can be used to calculate and characterise the geometries of TCRs to MHC class II antigen, however, due to the smaller number of complexes we have not fit parametric distributions to the geometry features, which means it is not possible to calculate a geometry score.


