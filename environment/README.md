## 📦 Setting Up the Conda Environment

This project uses a Conda environment defined in [`environment.yml`](./environment.yml). To set it up, follow these steps:

### 🛠️ 1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/)

Make sure you have Conda installed and available in your terminal:

```bash
conda --version  

### 📥 2. Create the environment

Run the following command from the root of the repository:

```bash
conda env create -f environment.yml   

This will create a new environment named hist-obs with all necessary dependencies.