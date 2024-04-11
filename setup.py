from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='TCRpy',
    version='0.1',
    description='Set of methods to parse, annotate, and calculate features of TCR structures',
    license='BSD 3-clause license',
    maintainer='Nele Quast',
    long_description=long_description,
    long_description_content_type='text/markdown',
    maintainer_email='quast@stats.ox.ac.uk',
    include_package_data=True,
    packages=find_packages(include=('TCRpy'), exclude=('test', 'test.*')),
    install_requires=[
        'numpy',
    ],
)