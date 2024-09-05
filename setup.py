from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="stcrpy",
    version="0.1.18",
    description="Set of methods to parse, annotate, and calculate features of TCR structures",
    license="BSD 3-clause license",
    maintainer="Nele Quast",
    long_description=long_description,
    long_description_content_type="text/markdown",
    maintainer_email="quast@stats.ox.ac.uk",
    include_package_data=True,
    packages=find_packages(".", exclude=("test", "test.*")),
    package_data={"stcrpy": ["tcr_geometry/reference_data/*.pdb"]},
    install_requires=[
        "biopython",
        "numpy==1.26.0",  # required for PLIP utils, which aren't compatible beyond this numpy version
        "lxml",
        "openbabel-wheel",
        "rdkit",
        "anarci-mhc==0.0.6",
        "pandas",
        "matplotlib",
        "scipy",
        "scikit-learn",
    ],
)
