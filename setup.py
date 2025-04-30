from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="stcrpy",
    version="1.0.3",
    description="Set of methods to parse, annotate, and calculate features of TCR structures",
    license_files=["./LICENCE", "./stcrpy/tcr_geometry/TCRCoM_LICENCE"],
    maintainer="Nele Quast",
    long_description=long_description,
    long_description_content_type="text/markdown",
    maintainer_email="quast@stats.ox.ac.uk",
    include_package_data=True,
    packages=find_packages(".", exclude=("test", "test.*")),
    package_data={"stcrpy": ["tcr_geometry/reference_data/*.pdb"]},
    install_requires=[
        "biopython",
        "numpy==1.26.4",  # required for pymol and scipy mutual compatability?
        "lxml",
        "openbabel-wheel==3.1.1.21",
        "rdkit",
        "anarci-mhc",
        "pandas",
        "matplotlib",
        "scipy",
        "requests",
        "scikit-learn",
    ],
)
