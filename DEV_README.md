# STCRpy Developer's README

These instructions are for anyone seeking to build and deploy, or extend and debug STCRpy. 

## Building and deploying
In case build and twine aren't installed: 
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

### Building the documentation
```
./docs/build_docs.sh
```
The bash script includes a push to github which will trigger the documentation build.

## Debugging
To install a development version of stcrpy and enable local debugging: 
```
git clone git@github.com:npqst/STCRpy.git
pip install -e .
ANARCI --build_models
```
To debug and develop features related to interaction profiling, visualisations in pymol or dataset generation, you will also need to install the relevant packages as described above. 


## Some reported gotchas
### ipykernel not working after installing pymol
* Symptoms: Jupyter notebook kernel is crashing or not starting, pymol is not importing. If the error refers to **sqlite** and looks something like the below, it is likely a dependency problem between ipykernel and pymol. 
```
ImportError: $PATH_TO_CONDA/envs/$ENV_NAME/lib/python3.12/lib-dynload/_sqlite3.cpython-312-x86_64-linux-gnu.so: undefined symbol: sqlite3_deserialize.
```

* Try the following: 
```
conda install "sqlite>=3.38"
```
You may need to reinstall ipykernel: 
```
conda uninstall ipykernel
conda install ipykernel
```


### VS code debugger not working after pymol installation
Sometimes installing pymol seems to cause the VS code debugger with break points to crash. To resolve this uninstall pymol. 
```
conda uninstall pymol
```