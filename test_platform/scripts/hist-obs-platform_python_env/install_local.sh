#!/bin/bash
set -e
set -x

PYTHON_VERSION=3.9.18
ENVIRONMENT_NAME=era_py${PYTHON_VERSION/./}

#save list of conda envs as temp file
conda env list > $HOME/.my_conda_envs

#if the current environment name is not in your environments, install it
if [[ "$(cat $HOME/.my_conda_envs)" != *"${ENVIRONMENT_NAME} "* ]]; then

  conda create -y --quiet -c conda-forge -n ${ENVIRONMENT_NAME} python=${PYTHON_VERSION}

  source activate $ENVIRONMENT_NAME

  # install mamba
  conda install -c conda-forge mamba --quiet --yes

  # install conda packages (channels managed in YAML file)
  # pip dependencies also managed in YAML file
  mamba env update -f requirements.yml

#if it is, just update the environment
else

  source activate $ENVIRONMENT_NAME
  
  mamba env update -f requirements.yml

  # activate the jupyter kernel
  ipython kernel install --name "${ENVIRONMENT_NAME}" --user

fi

#remove temp file
rm $HOME/.my_conda_envs
