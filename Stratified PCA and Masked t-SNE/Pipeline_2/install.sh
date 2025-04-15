#!/bin/bash

# Export environment variable
export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True

# Install pip packages
pip install matplotlib==3.3.4
pip install numpy==1.20.3
pip install pandas==1.3.5
pip install PyYAML==5.1.2
pip install scipy==1.5.3
pip install seaborn==0.11.2
pip install xgboost==1.1.1
pip install uncertainty-calibration==0.0.7
pip install scikit-allel==1.3.1
pip install scikit-learn==1.0.1
pip install sklearn-crfsuite==0.3.6
pip install ray

# Install mamba/conda packages
mamba install -y tqdm==4.62.3
mamba install -y -c anaconda openblas
mamba install -y -c bioconda bcftools
mamba install -y bioconda::gcta
mamba install -y conda-forge::r-base
mamba install -y conda-forge::r-tidyverse
mamba install -y conda-forge::r-rtsne
