#!/bin/bash

## This bash script should be used to update the documentation

# cd to docs if not already there
if [ "$(basename "$PWD")" != "docs" ]; then
    cd docs || { echo "Failed to change directory to docs"; exit 1; }
fi

sphinx-apidoc -o ./source ../stcrpy       # generate sphinx source from documentation

make clean      # clean up existing builds
make html       # build html files

cp -r build/html/* .build_html/         # copy static html to prebuilt repository

git add ./.build_html                
git commit -m 'updated docs'
git push                            # git push triggers new deploy by readthedocs

