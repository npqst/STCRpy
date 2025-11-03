# Docker Setup for STCRpy

This directory contains Docker configuration for running STCRpy in a containerized environment with uv dependency management.

## Quick Start

### Build and run the basic STCRpy container:
```bash
docker-compose up -d stcrpy
docker-compose exec stcrpy bash
```

### Build and run with ML dependencies (PyTorch, PyTorch Geometric):
```bash
docker-compose --profile ml up -d stcrpy-ml
docker-compose --profile ml exec stcrpy-ml bash
```

## Services

### `stcrpy` (default)
Basic STCRpy installation with core dependencies:
- BioPython
- PLIP for interaction profiling
- ANARCI for sequence annotation
- All core STCRpy functionality

### `stcrpy-ml` (ML profile)
STCRpy with machine learning dependencies:
- PyTorch
- PyTorch Geometric
- All features from base image
- Use for graph neural network applications

### `stcrpy-batch` (batch profile)
For batch processing tasks:
```bash
docker-compose --profile batch run stcrpy-batch python your_script.py
```

## Directory Structure

The following directories are mounted as volumes:
- `./data` → `/app/data` - Input data files
- `./examples` → `/app/examples` - Example scripts and notebooks
- `./output` → `/app/output` - Analysis outputs

Create these directories before running:
```bash
mkdir -p data output
```

## Usage Examples

### Interactive Python session:
```bash
docker-compose exec stcrpy python
```

```python
import stcrpy
tcr = stcrpy.fetch_TCRs("8gvb")
```

### Run a script:
```bash
docker-compose exec stcrpy python examples/STCRpy_examples.py
```

### Process local PDB files:
```bash
# Place your PDB files in ./data/
docker-compose exec stcrpy python -c "
import stcrpy
tcr = stcrpy.load_TCR('/app/data/your_file.pdb')
tcr.calculate_geometry()
tcr.save('/app/output/processed.pdb')
"
```

## Building Images

### Build all images:
```bash
docker-compose build
```

### Build specific image:
```bash
docker-compose build stcrpy
docker-compose build stcrpy-ml
```

## Managing Containers

### Stop containers:
```bash
docker-compose down
```

### View logs:
```bash
docker-compose logs -f stcrpy
```

### Remove all containers and volumes:
```bash
docker-compose down -v
```

## Notes

- **uv**: This setup uses [uv](https://github.com/astral-sh/uv) for fast Python package management
- **ANARCI models**: Built automatically during image creation (takes a few minutes)
- **PyMOL**: Not included in Docker images (requires GUI support). For visualization, export files and use local PyMOL installation
- **PLIP**: ✅ **Fully functional!** Interaction profiling is enabled using system OpenBabel packages (python3-openbabel) and PLIP source code
- **Performance**: First build may take 10-15 minutes due to ANARCI model building

## Troubleshooting

### ANARCI build fails:
The ANARCI model building is automatic but may fail on some systems. If needed, rebuild manually:
```bash
docker-compose exec stcrpy ANARCI --build_models
```

### Permission issues with volumes:
Ensure your user has write permissions to `./data` and `./output` directories.

### ML dependencies fail:
If PyTorch installation fails, you may need to specify the CPU version:
```bash
docker-compose exec stcrpy-ml uv pip install torch --index-url https://download.pytorch.org/whl/cpu
```
