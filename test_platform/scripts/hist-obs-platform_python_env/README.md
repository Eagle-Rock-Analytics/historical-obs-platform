# hist-obs-platform_python_env
Eagle Rock Analytics Anaconda Python Environment for hist-obs-platform

## Conda environment for climate science in Eagle Rock Analytics hist-obs-platform

This repository contains resources to facilitate installing a climate-science-ready conda environment from scratch on several machine environments.  It uses mamba to deal with the complex dependencies required by various climate-related packages (e.g., cartopy and geopandas).

* `requirements.yml` - contains a list of packages useful for climate science, intended to be installed via `mamba env update -f requirements.yml`
* `install_local.sh` - a bash script that automatically creates an environment named `era_py311` on the local machine. The anaconda base environment needs to be loaded for it to work.
