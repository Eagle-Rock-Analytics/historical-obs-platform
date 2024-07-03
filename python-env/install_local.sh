#!/bin/bash
set -e
set -x

# In case there is not a base anaconda distribution, download and install it. 
#curl -L -O "https://repo.anaconda.com/archive/Anaconda3-2024.02-1-$(uname)-$(uname -m).sh"
#bash Anaconda3-2024.02-1-$(uname)-$(uname -m).sh

# Activate base environment:
# This step depends on how you activate your base conda environment (either by conda init, .bashrc, .zshrc, .tcshrc, etc.
source ~/.bashrc
# source ~/.zshrc 


# Define python version and environment name
PYTHON_VERSION=3.9
ENVIRONMENT_NAME=hist-obs

# Create conda environment
conda create -y --quiet -c conda-forge -n ${ENVIRONMENT_NAME} python=${PYTHON_VERSION}
conda activate $ENVIRONMENT_NAME

# install mamba (mamba is a reimplementation of the conda package manager in C++
# https://github.com/mamba-org/mamba
conda install -c conda-forge mamba --quiet --yes

# install conda packages (channels managed in YAML file)
# pip dependencies also managed in YAML file
mamba env update -f requirements.yml

# Install kernel for jupyter and ipython
ipython kernel install --name "${ENVIRONMENT_NAME}" --user

# Clean conda and mamba (Remove unused packages and caches)
conda clean --all -y
mamba clean --all -y

