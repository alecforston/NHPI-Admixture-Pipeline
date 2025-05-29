# How to Use Pipeline 2
We created an environment called Capstone_Final. This environment is very finicky because you have to install an
older version of python to get the program Gnomix to work. 

## Running the Script
After installing dependencies, running the *African_PCA_Before_Masking.sh* file will generate all of the files and figures in pipeline 2.

**Run the command:**

```bash
sbatch African_PCA_Before_Masking.sh
```

Make sure to run this command in the directory:
~/path_to_your_files

That command being said, we have done two things to make it easier to validate our data reproducibility
including commenting out the line 

```python
python
~/groups/grp_Capstone_Pipeline/nobackup/archive/gnomix/gnomix.py
African_Subset_ld_filtered_data_LD_Update_VCF.vcf
~/groups/grp_Capstone_Pipeline/nobackup/archive/Matt_Code_Reproducibility/Pipeline_2/data_2
2 True
~/groups/grp_Capstone_Pipeline/nobackup/archive/gnomix/pretrained_gnomix_models/chr2/model_chm_2.pkl”
```

because it takes ~30 hours to run. This script takes in our phased and linkage disequilibrium filtered vcf file called
African_Subset_ld_filtered_data_LD_Update_VCF.vcf and produces three files in the
“~/groups/grp_Capstone_Pipeline/nobackup/archive/Matt_Code_Reproducibility/Pipeline_2/data_2”
directory called “query_file_phased.vcf”, “query_results.fb”, and
“query_results.msp” which are the inputs for the imputation step.

Imputing “query_file_phased.vcf” and  “query_results.msp” will result in the output.vcf
file in the folder: ~/groups/grp_Capstone_Pipeline/nobackup/archive/Matt_Code_Reproducibility/Pipeline_2/

Additionally, we ran our R script that builds our figures in
order to make it easier for you to visually see what our pipeline produces.

Rscript Plot.R

This script generates called "PCA_Not_Masked.png", "PCA_Masked.png", "T_SNE_Masked.png"
that can be downloaded and viewed on your own personal computer.

## Dependencies

```bash
conda create -n Capstone_Final python=3.7
```

After doing this, we then installed all of the needed dependencies using the following commands:

```bash
pip install matplotlib==3.3.4

pip install numpy==1.20.3

pip install pandas==1.3.5

pip install PyYAML==5.1.2

pip install scipy==1.5.3

pip install seaborn==0.11.2

pip install xgboost==1.1.1

mamba install tqdm==4.62.3

export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True

pip install uncertainty-calibration==0.0.7

pip install scikit-allel==1.3.1

pip install scikit-learn==1.0.1

pip install sklearn-crfsuite==0.3.6

mamba install -c anaconda openblas 

mamba install -c bioconda bcftools  

pip install ray

mamba install bioconda::gcta 

mamba install conda-forge::r-base

mamba install conda-forge::r-tidyverse

mamba install conda-forge::r-rtsne
```

You also need to download the program Impute5 using the command
in this directory:

```bash
 ~/path_to_your_files wget"https://www.dropbox.com/scl/fo/ukwimchnvp3utikrc3hdo/ALnlW6hpad9EjQ5Z3r6ffMw/impute5_v1.2.0.zip?rlkey=n2zty39bdst5j5tycd0sf89ee&e=1&dl=1"
-O impute5.zip

unzip impute5.zip
```
