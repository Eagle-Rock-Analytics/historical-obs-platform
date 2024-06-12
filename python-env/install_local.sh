#!/bin/bash
set -e
set -x

#curl -L -O "https://repo.anaconda.com/archive/Anaconda3-2024.02-1-$(uname)-$(uname -m).sh"
#bash Anaconda3-2024.02-1-$(uname)-$(uname -m).sh
source ~/.bashrc

PYTHON_VERSION=3.11
ENVIRONMENT_NAME=hist-obs

conda create -y --quiet -c conda-forge -n ${ENVIRONMENT_NAME} python=${PYTHON_VERSION}
conda activate $ENVIRONMENT_NAME
conda install -c conda-forge mamba --quiet --yes
mamba env update -f requirements.yml
ipython kernel install --name "${ENVIRONMENT_NAME}" --user
