#!/bin/bash --login

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10240M   # memory per CPU core

#SBATCH -J "Chunking And Imputing Final"   # job name
#SBATCH --mail-user=alecsf@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
conda activate Capstone_Final

##Converts VCF File to Plink Format 
plink --vcf ../african_subset_chr2.vcf.gz --make-bed --out Plink_African_Subset_Check_LD_Update

#Generates a list of SNPs that should be pruned called pruned_sites_African_Subset_LD_Update.prine.in
plink --bfile Plink_African_Subset_Check_LD_Update --indep-pairwise 50 10 .2 --out pruned_sites_African_Subset_LD_Update


##This takes the Plink file Plink_African_Subset_Check_LD_Update and the SNPs to be pruned in the file pruned_sites_African_Subset_LD_Update.prune.in
##and produces a Plink file called African_Subset_ld_filtered_data_LD_Update
plink --bfile Plink_African_Subset_Check_LD_Update --extract pruned_sites_African_Subset_LD_Update.prune.in --make-bed --out African_Subset_ld_filtered_data_LD_Update

##This command takes the Plink file input African_Subset_ld_filtered_data_LD_Update and recodes it to a VCF called African_Subset_ld_filtered_data_LD_Update_VCF
plink --bfile African_Subset_ld_filtered_data_LD_Update --recode vcf --out African_Subset_ld_filtered_data_LD_Update_VCF

#These commands take the Plink file African_Subset_ld_filtered_data_LD_Update and then produces pca eiganvalues and eiganvectors called GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated.eigenval and GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated.eigenvec. This is for our data before we have performed masking. 
gcta --bfile African_Subset_ld_filtered_data_LD_Update --make-grm --out GCTA_African_Subset_ld_filtered_data_LD_Update
gcta --grm GCTA_African_Subset_ld_filtered_data_LD_Update --pca 20 --out GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated

##This runs Gnomix on the vcf file called African_Subset_ld_filtered_data_LD_Update_VCF.vcf and prodces .msp map file and phased vcf file in data_2 directory. This command figures out which SNPs need to be masked.
#python ~/path_to_gnomix/gnomix.py African_Subset_ld_filtered_data_LD_Update_VCF.vcf /data_2 2 True ~/path_to_model/model_chm_2.pkl

# impute data, run T-sne, and produce figures
sbatch Masking_and_Figures.sh
