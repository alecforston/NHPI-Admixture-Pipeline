#!/bin/bash --login

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=250G   # memory per CPU core
#SBATCH -J "Gnomix"   # job name


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
conda activate Environment_name

python ~/path_to_gnomix/gnomix.py African_Subset_ld_filtered_data_LD_Update_VCF.vcf ~/path_to_file/data_2 2 True ~/path_to_model/model_chm_2.pkl
