## Masking
.msp cleaner:
query_file_phased.vcf  query_results.fb  query_results.msp
#### Clean the .MSP File (msp_clean_job.sh / mspCleaner.py)
python mspCleaner
##### Header: (msp_cleaner_header.py)

```python
file_path = "../data_2/query_results.msp"
output_path = "query_clean_header.msp"

# Read file into a list of lines
with open(file_path, "r") as file:
    lines = file.readlines()

# Skip the original first line
lines = lines[1:]

# Clean the new first line
first_line_tokens = lines[0].strip().split()
cleaned_tokens = [token.split('_')[-1] if '_' in token else token for token in first_line_tokens]
cleaned_first_line = ' '.join(cleaned_tokens)

# Replace the first line
lines[0] = cleaned_first_line + '\n'

# Write to output file (no pandas, no quotes)
with open(output_path, "w") as out_file:
    out_file.writelines(lines)
```

##### New
```python
import pandas as pd
from io import StringIO
import re

file_path = "query_clean_header.msp"

with open(file_path, 'r') as f:
    lines = f.readlines()

# Clean and split the header
header_line = lines[0].strip().replace('n snps', 'n_snps').replace('.0', '')
columns = header_line.split()

# Load the rest of the file with regex-based whitespace splitting
data_str = ''.join(lines[1:])
df = pd.read_csv(
    StringIO(data_str),
    sep=r'\s+',
    engine='python',
    names=columns,
)

# Drop columns ending in ".1"
columns_to_drop = [col for col in df.columns if re.search(r'\.1', col)]
df.drop(columns=columns_to_drop, inplace=True)

df.to_csv("query_results_cleaned_final.msp", sep='\t', index=False)
```
##### Old
```python
import pandas as pd
from io import StringIO

file_path = "query_results.msp"

with open(file_path, 'r') as f:
    lines = f.readlines()

# Clean and split the header

header_line = lines[0].strip().replace('n snps', 'n_snps')
columns = header_line.split()

# Load the rest of the file with regex-based whitespace splitting
data_str = ''.join(lines[1:])
df = pd.read_csv(
    StringIO(data_str),
    sep=r'\s+',
    engine='python',
    names=columns,
)

# Drop columns ending in ".1"
columns_to_drop = [col for col in df.columns if col.endswith('.1')]
df.drop(columns=columns_to_drop, inplace=True)

df.to_csv("query_results_cleaned.msp", sep='\t', index=False)
```
#### Clean the VCF File
```python
file_path = "../data_2/query_file_phased.vcf"
output_path = "chr2_phased.vcf"

# Read file into a list of lines
with open(file_path, "r") as file:
    lines = file.readlines()

# Skip the original first line
#lines = lines[1:]

# Clean the new first line
first_line_tokens = lines[9].strip().split()
cleaned_tokens = [token.split('_')[-1] if '_' in token else token for token in first_line_tokens]
cleaned_first_line = ' '.join(cleaned_tokens)

# Replace the first line
lines[9] = cleaned_first_line + '\n'

# Write to output file
with open(output_path, "w") as out_file:
    out_file.writelines(lines)
```

```python
file_path = "../data_2/query_file_phased.vcf"
output_path = "chr2_phased.vcf"

with open(file_path, "r") as file:
    lines = file.readlines()

# Clean the header
first_line_tokens = lines[9].strip().split()
cleaned_tokens = [token.split('_')[-1] if '_' in token else token for token in first_line_tokens]
# Drop the last column
cleaned_tokens = cleaned_tokens[:-1]
# Create the cleaned line
cleaned_first_line = ' '.join(cleaned_tokens)
# Replace the line in the list
lines[9] = cleaned_first_line + '\n'

# Write to output file
with open(output_path, "w") as out_file:
    out_file.writelines(lines)
```


#### Mask the VCF (masking_job.sh)
##### New Command
python masking.py ../data_2/query_file_phased.vcf masked_ALL.chr2.vcf query_results_cleaned_final.msp 3
##### Bash job file:
```bash
#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=200G   # memory per CPU core
#SBATCH -J "masking"   # job name
#SBATCH --mail-user=alecsf@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
python masking.py ../data_2/query_file_phased.vcf masked_ALL.chr2.vcf query_results_cleaned_final.msp 3
```
##### Python File:
```python
import pandas as pd
import sys
import os

def print_usage():
    print("Usage: python masking2.py <vcf_file> <output_file> <msp_file> <masking_number>")
    print("\nArguments:")
    print("  <vcf_file>       : The input VCF file to be processed.")
    print("  <output_file>    : The output VCF file where results will be saved. This file will be overwritten.")
    print("  <msp_file>       : The masking position file containing SNP masking information.")
    print("  <masking_number> : The integer corresponding to the group value in the .msp file. All other numbers will be masked.\n")

# Check if the correct number of arguments are passed
if len(sys.argv) != 5:
    print_usage()
    sys.exit(1)

vcf_file = sys.argv[1]
output_file = sys.argv[2]
msp_file = sys.argv[3]

try:
    masking_number = int(sys.argv[4])
except ValueError:
    print("Error: masking_number must be an integer.")
    sys.exit(1)

# Check if input files exist
if not os.path.exists(vcf_file):
    print(f"Error: The file '{vcf_file}' does not exist.")
    sys.exit(1)

if not os.path.exists(msp_file):
    print(f"Error: The file '{msp_file}' does not exist.")
    sys.exit(1)

with open(vcf_file, 'r') as vcf_in, open(output_file, 'w') as vcf_out:
    masking_df = pd.read_csv(msp_file, sep="\s+")
    columns = masking_df.columns.tolist()
    # Increment index each time the POS is greater than the epos
    msp_index = 0

    for line in vcf_in:
        if not line.strip():  # Skip empty lines
            continue
        if line.startswith("##"):  # Keep VCF metadata
            vcf_out.write(line)
            continue
        elif line.startswith("#"):  # Process header line
            vcf_headers = line.strip().split()
            sample_index = {sample: idx for idx, sample in enumerate(vcf_headers[9:])}
            vcf_out.write(line)
            continue

        fields = line.strip().split()

		# Ensure the line has enough fields (should be at least 9 fields)
        if len(fields) < 9:
            continue

		pos = int(fields[1])

        if pos > masking_df['epos'].iloc[msp_index] and msp_index < len(masking_df):
            msp_index += 1

        if masking_df['spos'].iloc[msp_index] <= pos <= masking_df['epos'].iloc[msp_index]:
            for sample in vcf_headers[9:]:
                if sample in masking_df.columns and masking_df[sample].iloc[msp_index] != masking_number:
                    idx = sample_index[sample] + 9  # Adjust for VCF column indexing
                    fields[idx] = ".|."

        if msp_index == len(masking_df):
            break                        

		# Write the modified or original fields to the output file
        vcf_out.write("\t".join(fields) + "\n")

with open(output_file, 'r+') as file:
    lines = file.readlines()
    if len(lines) > 0:
        file.seek(0)
        file.writelines(lines[:-1])  # Remove the last line
        file.truncate()
```

Modern Pacific Islanders are often admixed, possessing European and occasionally Native American and African ancestries. Because ancestries introduced via colonial settlement did not necessarily follow the same island settlement process (or founder sizes and dates) as the original Polynesian settlement, such ancestries need to be distinguished, necessitating an ancestry-specific approach. 

For this reason we removed European chromosomal segments, as well as African and Native American, from the Pacific island samples. 

We refer to the remaining (unmasked) chromosomal segments as Polynesian ancestry chromosomal segments, and we refer to analyses that use only these segments as Polynesian ancestry-specific analyses. Non-Polynesian populations are included as references and are not masked. The main two tables use masking. 
## Impute 5
### Process:

#### Download Impute5
```bash
wget "https://www.dropbox.com/scl/fo/ukwimchnvp3utikrc3hdo/ALnlW6hpad9EjQ5Z3r6ffMw/impute5_v1.2.0.zip?rlkey=n2zty39bdst5j5tycd0sf89ee&e=1&dl=1" -O impute5.zip
unzip impute5.zip
```
#### Copy over files
```bash
cd ../../impute5_v1.2.0/data
cp ../../gnomix/data/chr2_phased.vcf chr2_phased.vcf
cp ../../gnomix/data/query_results_cleaned_final.msp query_results_cleaned_final.msp
```
#### Clean vcf file:
```python
input_file = "chr2_phased.vcf"

with open(input_file, 'r') as f:
        content = f.read()
# Replace spaces with tabs
content_with_tabs = content.replace(' ', '\t')
with open(input_file, 'w') as f:
        f.write(content_with_tabs)
```
#### Get References:
```bash
awk 'NR==1 {for (i=7; i<=NF; i++) header[i]=$i; next} {for (i=7; i<=NF; i++) count[i]+=$i==3} END {for (i in count) print header[i], count[i]}' query_results.msp | sort -k2 -nr | head -n 803 | awk '{print $1}' | sed 's/\.[0-9]*$//' | uniq > reference_individuals.txt
```
New:
```bash
awk 'NR==1 {for (i=7; i<=NF; i++) header[i]=$i; next} {for (i=7; i<=NF; i++) count[i]+=$i==3} END {for (i in count) print header[i], count[i]}' query_results_cleaned_final.msp  | sort -k2 -nr | head -n 35 | awk '{print $1}' | uniq > linkage_reference_individuals.txt
```
#### Make Reference VCF File (30 minutes)
```bash
bcftools view -S reference_individuals.txt -o impute_reference.vcf ALL.chr2.vcf
```
New:
```bash
bcftools view -S linkage_reference_individuals.txt -o phased_reference.vcf chr2_phased.vcf
```
#### Remove References From Masked VCF File (30 minutes)

```bash
bcftools view -S ^reference_individuals.txt -o filtered_masked_ALL.chr2.vcf -O v masked_ALL.chr2.vcf
```
New:
```bash
bcftools view -S ^linkage_reference_individuals.txt -o filtered_chr2_phased.vcf -O v chr2_phased.vcf
```
#### Convert to Binary Files (30 minutes)

```bash
bcftools view -Ob -o impute_reference.bcf impute_reference.vcf
bcftools index impute_reference.bcf
bcftools view -Ob -o filtered_masked_ALL.chr2.bcf filtered_masked_ALL.chr2.vcf
bcftools index filtered_masked_ALL.chr2.bcf
```

New:
```bash
bcftools view -Ob -o phased_reference.bcf phased_reference.vcf
bcftools index phased_reference.bcf
bcftools view -Ob -o filtered_chr2_phased.bcf filtered_chr2_phased.vcf
bcftools index filtered_chr2_phased.bcf
```
#### Move back to Impute
```bash
cd ../../impute5_v1.2.0/data/
```
#### Copy the Files Over
```bash
cp ../../gnomix/data/filtered_masked_ALL.chr2.bcf .
cp ../../gnomix/data/filtered_masked_ALL.chr2.bcf.csi .
cp ../../gnomix/data/impute_reference.bcf .
cp ../../gnomix/data/impute_reference.bcf.csi .
```
Other files to potentially copy:
```bash
cp ../../gnomix/data/query_results.msp .
cp ../../gnomix/data/ALL.chr2.vcf .
```
#### Chunk & Impute
##### Old Chunker
```bash
../imp5Chunker_v1.2.0_static --h ./impute_reference.bcf --g ./filtered_masked_ALL.chr2.bcf --o ./IFM_CHUNKER_ALL.txt --r 2

```
##### New Chunker
```bash
../imp5Chunker_v1.2.0_static --h ./phased_reference.bcf --g ./filtered_chr2_phased.bcf --o ./phased_CHUNKER_ALL.txt --r 2

```
##### Old
```bash
#!/bin/bash

# Extract dynamic start position from query_results.msp
START=$(head -n 2 ../../gnomix/data/query_results.msp | tail -n 1 | awk -F'\t' '{print $2}')

# Extract dynamic end position from VCF
END=$(tail -n 1 ../../gnomix/data/ALL.chr2.vcf | awk -F'\t' '{print $2}')

# Add a buffer of 100000 around the region
BUFFER_START=$((START - 1000))
BUFFER_END=$((END + 1000))

# Define the chromosome
CHR=2

../xcftools_static view --i impute_reference.bcf -o reference_xcf.bcf -O sh -r 20 -T 8 -m 0.03125

../impute5_v1.2.0_static --h impute_reference_xcf.bcf --g filtered_masked_ALL.chr2.bcf --o imputed_filtered_masked_ALL.chr2.bcf  --r ${CHR}:${START}-${END} --buffer-region ${CHR}:${BUFFER_START}-${BUFFER_END}

```
##### Unused
```bash
#!/bin/bash

# Extract dynamic start position from query_results.msp
START=$(head -n 2 query_results.msp | tail -n 1 | awk -F'\t' '{print $2}')

# Extract dynamic end position from VCF
END=$(tail -n 1 ALL.chr2.vcf | awk -F'\t' '{print $2}')

# Add a buffer of 100000 around the region
BUFFER_START=$((START - 1000))
BUFFER_END=$((END + 1000))

# Define the chromosome
CHR=2

# Run impute5
../impute5_v1.2.0_static \
--h impute_reference_xcf.bcf 
--g filtered_masked_ALL.chr2.bcf
--o imputed_filtered_masked_ALL.chr2.bcf 
--r ${CHR}:${START}-${END} \
--buffer ${CHR}:${BUFFER_START}-${BUFFER_END}
```
##### Chunker & Impute Old
```bash
#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=200G   # memory per CPU core
#SBATCH -J "Chunking And Imputing Final"   # job name
#SBATCH --mail-user=alecsf@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#Chunker Step
../imp5Chunker_v1.2.0_static --h ./impute_reference.bcf --g ./filtered_masked_ALL.chr2.bcf --o ./IFM_CHUNKER_ALL.txt --r 2

#xcftools step
../xcftools_static view --i impute_reference.bcf -o impute_reference_xcf.bcf -O sh -r 2 -T 8 -m 0.03125

#Imputation Based on the chunks
# Path to your chunk file
CHUNK_FILE="IFM_CHUNKER_ALL.txt"

# Path to the impute5 binary
IMPUTE5_BIN="../impute5_v1.2.0_static"

REFERENCE_FILE="impute_reference_xcf.bcf"
TARGET_FILE="filtered_masked_ALL.chr2.bcf"

# Output directory
OUTPUT_DIR="imputed_chunks"
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

    "$IMPUTE5_BIN" \
        --h "$REFERENCE_FILE" \
        --g "$TARGET_FILE" \
        --o "$output_file" \
        --r "${chr_id}:${imp_start}-${imp_end}" \
        --buffer-region "${chr_id}:${buf_start}-${buf_end}"

done < "$CHUNK_FILE"
```
../impute5_v1.2.0_static \
--h impute_reference_xcf.bcf 
--g filtered_masked_ALL.chr2.bcf
--o imputed_filtered_masked_ALL.chr2.bcf 
--r 20:10554-243188367 
--buffer 20:10054-243188867
##### Chunker & Impute New
```bash
#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=200G   # memory per CPU core
#SBATCH -J "Phased Chunking And Imputing Final"   # job name
#SBATCH --mail-user=alecsf@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#Chunker Step
../imp5Chunker_v1.2.0_static --h ./phased_reference.bcf --g ./filtered_chr2_phased.bcf --o ./phased_CHUNKER_ALL.txt --r 2

#xcftools step
../xcftools_static view --i phased_reference.bcf -o phased_reference_xcf.bcf -O sh -r 2 -T 8 -m 0.03125

#Imputation Based on the chunks
CHUNK_FILE="phased_CHUNKER_ALL.txt"

# Path to the impute5 binary
IMPUTE5_BIN="../impute5_v1.2.0_static"

REFERENCE_FILE="phased_reference_xcf.bcf"
TARGET_FILE="filtered_chr2_phased.bcf"

# Output directory
OUTPUT_DIR="imputed_phased_chunks"
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

    "$IMPUTE5_BIN" \
        --h "$REFERENCE_FILE" \
        --g "$TARGET_FILE" \
        --o "$output_file" \
        --r "${chr_id}:${imp_start}-${imp_end}" \
        --buffer-region "${chr_id}:${buf_start}-${buf_end}"

done < "$CHUNK_FILE"
```
#### Combine Chunked bcf files
```bash
ls imputed_phased_chunks/imputed_chunk_*.bcf > bcf_list.txt
bcftools concat -Ob -o combined.bcf -f bcf_list.txt
```
#### Recombine VCF with Reference VCF
```bash
bcftools merge -o admixture_accounted_ALL.Chr2.bcf -O v imputed_filtered_masked_ALL.chr2.bcf impute_reference.vcf
```


### Other
https://jmarchini.org/software/#impute-5
Sample Run:
impute5 --h reference_xcf.bcf --m chr20.b37.gmap.gz --g target.bcf --r
20:2000000-7000000 --o imputed.bcf



impute5_v1.2.0 --h /test/reference.bcf --
g /test/target.bcf --m
/test/chr20.b37.gmap.gz --o home/imputed.bcf
--r 20:1000000-4000000


./impute5_v1.2.0_static --h ./test/reference.bcf --g ./test/target.bcf --m ./test/chr20.b37.gmap.gz --o ./imputed.bcf --r 20:1100000-4000000 --buffer 20:1000000-4100000


./impute5_v1.2.0_static 
--h ./test/reference.bcf 
--g ./test/target.bcf 
--m ./test/chr20.b37.gmap.gz 
--o ./imputed.bcf
--r 20:1100000-4000000 
--buffer 20:1000000-4100000

### Command:
./impute5_v1.2.0_static 
--h ./test/reference.bcf 
--g ./test/target.bcf 
--o ./imputed.bcf
--r 20:1100000-4000000 
--buffer 20:1000000-4100000


remove references from main file:
bcftools view -S ^reference_individuals.txt -o filtered_masked_ALL.chr2.vcf -O v masked_ALL.chr2.vcf

### --h
Reference

how to get the top  awk 

```bash
awk 'NR==1 {for (i=7; i<=NF; i++) header[i]=$i; next} {for (i=7; i<=NF; i++) count[i]+=$i==3} END {for (i in count) print header[i], count[i]}' query_results.msp | sort -k2 -nr | head -n 803 | awk '{print $1}' | sed 's/\.[0-9]*$//' | uniq > reference_individuals.txt
```

The indexes of the individuals with the most African DNA will be used as our sample



bcftools view -S reference_individuals.txt -o impute_reference.vcf ALL.chr2.vcf

^ Takes about 20 minutes __make_reference_vcf.sh__

c
^ Converting all to bcf and indexing takes 30 minutes __convert_vcf_to_bcf.sh__

cp ../../gnomix/data/ALL.chr2.bcf



../impute5_v1.2.0_static --h reference_xcf.bcf --g ALL.chr2.bcf --o imputed.bcf --r 20:10554-243188367 --buffer-region 20:10054-243188867

### --o
imputed{chromosome}.bcf
### --r 
chr: input

Get the Lower Int:

	head -n 2 query_results.msp | tail -n 1 |awk -F'\t' '{print $2}' 
__10554__

Get the Upper Int

	tail -n 1 ALL.chr2.vcf | awk -F'\t' '{print $2}'
__243188367__
{chr}:{lower}-{upper}
### --buffer
{chr}:{lower- 500}-{upper - 500}
### Help
Input:
  --h arg                             Haplotype reference panel in VCF/BCF/XCF
                                      format (must have .vcf[.gz]/.bcf
                                      extension). The file must be indexed
                                      (tabix/imp5 index).
  --g arg                             File containing target haplotypes for a
                                      study cohort that you want to impute in
                                      VCF/BCF/XCF format (must have
                                      .vcf[.gz]/.bcf extension). The file must
                                      be indexed (csi/tabix index).
  --m arg                             Fine-scale recombination map for the
                                      region to be analyzed. If not specified,
                                      a constant recombination rate of 1cM per
                                      Mb is used.
  --r arg                             Region to be imputed. Example --r
                                      chr20:1000000-5000000 (region within
                                      chromosome 20). Buffer-region expand this
                                      region.
  --buffer-region arg                 Region to be imputed with a buffer at
                                      both the borders. Buffer-region must
                                      include the imputed region (--r).
  --sparse-maf arg (=0.03125)         (Expert setting) Rare variant threshold
  --haploid                           Specifies the target samples are haploid
                                      in the region (e.g. for males in non-PAR
                                      region on Chromosome X).

Model parameters:
  --ne arg (=1000000)                 Effective population size.
  --err-imp arg (=9.99999975e-05)     (Expert setting) Imputation HMM error
                                      rate

State selection:
  --max-pbwt-depth arg (=10)          Max depth of PBWT indexes to condition on
  --min-pbwt-depth arg (=2)           Min depth of PBWT indexes to condition on
  --pbwt-cm arg (=0.0199999996)       Frequency of the selection algorithm in
                                      cM
  --Kpbwt arg (=1500)                 Max number of states to condition on
  --neigh-select                      Use only the positional prefix array to
                                      select states.

Test statistics parameters:
  --surfbat [=arg(=paired)]           Outputs FORMAT/SAP field (standard
                                      haploid dosages and surrogate family
                                      haploid dosages for surrogate family
                                      based association testing) and paired
                                      p-value derived from chisquare
                                      distribution with one degree of freedom.
                                      . Use only if if the SNP array data has
                                      been phased on the whole chromosome as
                                      not supported during ligation.
  --surfbat-maf arg (=0.00999999978)  Surfbat p-value restricted to sites with
                                      MAF above specified parameter
  --surfbat-info arg (=0.300000012)   Surfbat p-value restricted to sites with
                                      INFO score above specified parameter

Output:
  --o arg (=impute5.out.bcf)          Specifies output file name. Accepted
                                      extensions: [.vcf[.gz],.bcf,.bgen]. If
                                      the format is in BGEN format, by default
                                      the file in unphased. To output a phased
                                      BGEN file, specify the --out-ap-field
                                      option.
  --l arg                             Location of the log file to be written.
                                      If not specified, only console output is
                                      generated.
  --no-out-index                      Skip the calculation of the csi index
  --no-out-gp-field                   Do not output FORMAT/GP field (genotype
                                      posterior probabilities) if output is in
                                      VCF/BCF format.
  --no-out-ds-field                   Do not output FORMAT/DS field (genotype
                                      dosage) if output is in VCF/BCF format.
  --out-buffer                        Output SNP array variants in the buffer
                                      (necessary if the SNP array data has been
                                      phased in chunks rather than on the whole
                                      chromosome)
  --out-ap-field                      Output FORMAT/AP field (ALT haplotype
                                      probabilities) if output is in VCF/BCF
                                      format or outputs phased BGEN file. Use
                                      only if if the SNP array data has been
                                      phased on the whole chromosome as not
                                      supported during ligation.
  --bgen-bits arg (=8)                (Expert setting) Only used toghether when
                                      the output is in BGEN file format.
                                      Specifies the number of bits to be used
                                      for the encoding probabilities of the
                                      output BGEN file. If the output is in the
                                      .vcf[.gz]/.bcf format, this value is
                                      ignored. Accepted values: 1-32
  --bgen-compr arg (=zstd)            (Expert setting) Only used toghether when
                                      the output is in BGEN file format.
                                      Specifies the compression of the output
                                      BGEN file. If the output is in the
                                      .vcf[.gz]/.bcf format, this value is
                                      ignored. Accepted values: [no,zlib,zstd]

Other parameters:
  --help                              Produce help message, listing all the
                                      accepted arguments
  --contigs-fai arg                   If specified, header contig names and
                                      their lengths are copied from the
                                      provided fasta index file (.fai) instead
                                      of being taken from the target panel
                                      (default behavior). This allows to create
                                      files with all the contigs in the header
                                      (in the case the contigs in the reference
                                      panel are limited to a single chromosome)
                                      and therefore quickly merge
                                      chromosome-level files with bcftools
                                      merge --naive
  --threads arg (=1)                  Number of threads.
  --seed arg (=42)                    Seed for RNG.
  --estimate-mem-usage                Experimental. Estimates the dynamic
                                      memory usage.

All five options are mandatory and their descriptions are:
### --h (required)
#### VCF Reference File
--h specifies the haplotype reference panel in VCF/BCF/XCF format (must
have .vcf[.gz]/.bcf extension). The file must be indexed. The reference panel should
be phased and non-missing at every position.
### gMAP file (not required)
 --m specifies the fine-scale recombination map for the region to be analyzed. Maps
for humans can be found HERE. In the case this parameter is not defined, a constant
recombination rate is assumed.

pos	chr	cM
12994	2	0.000000
15491	2	0.000848
15672	2	0.000908

### -- g (required)
--g specifies the file containing target haplotypes for a study cohort that you want to
impute in VCF/BCF format (must have .vcf[.gz]/.bcf extension). The file must be
indexed (tabix index). The target dataset should be phased and non-missing in the
set of markers specified. Markers that are only present in the reference panel and not
in the target set, are imputed.
### --r (required)
specifies the target region or chromosome to be imputed . Buffer parameters will
expand this region, if specified.

20:2000000-7000000
### --o (required)
specifies the output filename. A proper extension is mandatory.
IMPUTE5 considers as genotype markers the markers in the intersection bet