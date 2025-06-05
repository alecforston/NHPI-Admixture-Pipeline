#!/bin/bash

conda activate Capstone_Final

# Convert VCF File to Plink Format 
plink --vcf ../african_subset_chr2.vcf.gz --make-bed --out Plink_African_Subset_Check_LD_Update

# Generate a list of SNPs that should be pruned called pruned_sites_African_Subset_LD_Update.prine.in
plink --bfile Plink_African_Subset_Check_LD_Update --indep-pairwise 50 10 .2 --out pruned_sites_African_Subset_LD_Update

##This takes the Plink file Plink_African_Subset_Check_LD_Update and the SNPs to be pruned in the file pruned_sites_African_Subset_LD_Update.prune.in
##and produces a Plink file called African_Subset_ld_filtered_data_LD_Update
plink --bfile Plink_African_Subset_Check_LD_Update --extract pruned_sites_African_Subset_LD_Update.prune.in --make-bed --out African_Subset_ld_filtered_data_LD_Update

# Takes the Plink file input African_Subset_ld_filtered_data_LD_Update and recode it to a VCF called African_Subset_ld_filtered_data_LD_Update_VCF
plink --bfile African_Subset_ld_filtered_data_LD_Update --recode vcf --out African_Subset_ld_filtered_data_LD_Update_VCF

# Take the Plink file African_Subset_ld_filtered_data_LD_Update and produce pca eiganvalues and eiganvectors called GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated.eigenval and GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated.eigenvec.
# This is for our data before we have performed masking. 
gcta --bfile African_Subset_ld_filtered_data_LD_Update --make-grm --out GCTA_African_Subset_ld_filtered_data_LD_Update
gcta --grm GCTA_African_Subset_ld_filtered_data_LD_Update --pca 20 --out GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated

##This runs Gnomix on the vcf file called African_Subset_ld_filtered_data_LD_Update_VCF.vcf and prodces .msp map file and phased vcf file in data_2 directory. This command figures out which SNPs need to be masked.
python ~/path_to_gnomix/gnomix.py African_Subset_ld_filtered_data_LD_Update_VCF.vcf /data_2 2 True ~/path_to_model/model_chm_2.pkl

# Impute data, run T-sne, and produce figures
# Clean the .msp headers:
python ~/groups/grp_Capstone_Pipeline/nobackup/archive/Matt_Code_Reproducibility/Pipeline_2/msp_cleaner1.py

# Clean the .msp header again and keep the relevant data
python ~/groups/grp_Capstone_Pipeline/nobackup/archive/Matt_Code_Reproducibility/Pipeline_2/msp_cleaner2.py

# Clean the vcf file py1
python ~/groups/grp_Capstone_Pipeline/nobackup/archive/Matt_Code_Reproducibility/Pipeline_2/vcf_cleaner1.py

# Clean the vcf file py2
python ~/groups/grp_Capstone_Pipeline/nobackup/archive/Matt_Code_Reproducibility/Pipeline_2/vcf_cleaner2.py

# Download Impute5
wget "https://www.dropbox.com/scl/fo/ukwimchnvp3utikrc3hdo/ALnlW6hpad9EjQ5Z3r6ffMw/impute5_v1.2.0.zip?rlkey=n2zty39bdst5j5tycd0sf89ee&e=1&dl=1" -O impute5.zip
unzip impute5.zip

# Get References
awk 'NR==1 {for (i=7; i<=NF; i++) header[i]=$i; next} {for (i=7; i<=NF; i++) count[i]+=$i==3} END {for (i in count) print header[i], count[i]}' data_2/query_results_cleaned_final.msp  | sort -k2 -nr | head -n 35 | awk '{print $1}' | uniq > data_2/linkage_reference_individuals.txt

# Load bcftools
module load bcftools/1.19-7gzjaqo

# Make Reference VCF File
bcftools view -S data_2/linkage_reference_individuals.txt -o data_2/phased_reference.vcf data_2/cleaned_phased.vcf

# Remove References From Masked VCF File
bcftools view -S data_2/linkage_reference_individuals.txt -o data_2/main.vcf -O v data_2/cleaned_phased.vcf

# Convert References to Binary Files
bcftools view -Ob -o data_2/phased_reference.bcf data_2/phased_reference.vcf

# Index Reference Binary File
bcftools index data_2/phased_reference.bcf

# Convert Main Data to Binary File
bcftools view -Ob -o data_2/main.bcf data_2/main.vcf

# Index Main Data file
bcftools index data_2/main.bcf

# Chunker Step: To impute, the chromosome must be cut into chunks
impute5_v1.2.0/imp5Chunker_v1.2.0_static --h data_2/phased_reference.bcf --g data_2/main.bcf --o data_2/phased_CHUNKER_ALL.txt --r 2

# Xcftools step
impute5_v1.2.0/xcftools_static view --i data_2/phased_reference.bcf -o data_2/phased_reference_xcf.bcf -O sh -r 2 -T 8 -m 0.03125

#Imputation Based on the chunks
CHUNK_FILE="data_2/phased_CHUNKER_ALL.txt"

# Path to the Impute5 binary
IMPUTE5_BIN="impute5_v1.2.0/impute5_v1.2.0_static"

REFERENCE_FILE="data_2/phased_reference_xcf.bcf"
TARGET_FILE="data_2/main.bcf"

# Output directory
OUTPUT_DIR="data_2/imputed_phased_chunks"
mkdir -p "$OUTPUT_DIR"

# Read the chunk file line by line
while IFS=$'\t' read -r chunk_id chr_id buffered_region imp_region length num_target_markers num_ref_markers; do
    # Extract start/end positions from regions
    imp_start=$(echo "$imp_region" | cut -d':' -f2 | cut -d'-' -f1)
    imp_end=$(echo "$imp_region" | cut -d':' -f2 | cut -d'-' -f2)
    buf_start=$(echo "$buffered_region" | cut -d':' -f2 | cut -d'-' -f1)
    buf_end=$(echo "$buffered_region" | cut -d':' -f2 | cut -d'-' -f2)

    output_file="${OUTPUT_DIR}/imputed_chunk_${chunk_id}.bcf"

    echo "Running imputation for chunk ${chunk_id}..."
    echo "  Imputation region: ${chr_id}:${imp_start}-${imp_end}"
    echo "  Buffer region:     ${chr_id}:${buf_start}-${buf_end}"
    echo "  Output file:       ${output_file}"

    # Impute each chunk
    "$IMPUTE5_BIN" \
        --h "$REFERENCE_FILE" \
        --g "$TARGET_FILE" \
        --o "$output_file" \
        --r "${chr_id}:${imp_start}-${imp_end}" \
        --buffer-region "${chr_id}:${buf_start}-${buf_end}"

done < "$CHUNK_FILE"

# Make a list of the chunked files
ls data_2/imputed_phased_chunks/imputed_chunk_*.bcf > data_2/bcf_list.txt

# Combine chunked files into one .bcf file
bcftools concat -Ob -o data_2/combined.bcf -f data_2/bcf_list.txt

# Convert bcf to vcf
bcftools view -O v -o output.vcf data_2/combined.bcf

##Convert vcf to Plink
plink --vcf output.vcf --make-bed --out Masked_Plink

##perform PCA on masked Plink file
gcta --bfile Masked_Plink --make-grm --out GCTA_Masked_Plink
gcta --grm GCTA_Masked_Plink --pca 20 --out GCTA_Masked_Plink_Final


##Make Figures and perform T-SNE
Rscript Plot.R
