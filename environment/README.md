# 📦 Setting Up the Conda Environment

This project uses a Conda environment defined in [`environment.yml`](./environment.yml). To set it up, follow these steps:

## 🛠️ 1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/)

You likely already have Conda installed if you have worked with conda environments before. To check that you have Conda installed and available in your terminal:

```bash
conda --version
```

## 📥 2. Create the environment

Run the following command from the `environment/` directory of the repository:

```bash
conda create -n hist-obs -y
```

This will create a new (but empty) environment named `hist-obs` 

## 🧪 3. Activate the environment

```bash
conda activate hist-obs
```

## 🛠️ 4. Install mamba in the environment 
Mamba is a fast replacement for conda that uses a C++ core for much quicker environment solving and package installation. We recommend using mamba instead of conda to speed up the setup process, especially due to the complexity of the dependency stack for geospatial packages. 

```bash
conda install mamba -c conda-forge -y
```

## 🛠️ 5. Install the required packages in the environment
Mamba will install all the packages in the `environment.yml` file (which includes python) into the `hist-obs` environment. 

```bash
mamba env update --file environment.yml --prune -y
```

## 🧹 6. Clean up unnecessary packages and caches
```bash
conda clean --all -y
```
You're now ready to run code using the hist-obs environment!
