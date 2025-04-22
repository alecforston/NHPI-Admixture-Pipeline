# How to use Pipeline 1
We created a conda environment with the latest versions of both PLINK and pandas, and activated the environment. 

## Installing Dependencies
conda create –n environment_name
conda activate environment_name
conda install bioconda::plink
conda install anaconda::pandas

In our working directory, we made two subdirectories, named “inputs” and “intermediates.” The names of these subdirectories
are included in the script, so these directories must be in the directory where the script is run. These subdirectories 
should already be in the Pipeline_1 directory.

For the sake of clarity, the intermediates subdirectory should be empty before running pipeline.sh. The inputs subdirectory 
should have 4 files: hapmap subset .bed, .bim and .fam files, which are PLINK output files, and a hapmap subset “3 columns” file,
which includes the FID, IID and self-identified ancestry of every individual in the cohort. 

## Running the Script

We ran the pipeline.sh script, which does all the work. The arguments of the script are as follows, excluding the name of the
script: “p” specifies that the input files are PLINK output files. Next is the prefix for the plink input files. The final argument
is the path to the “3 columns” file described above. Thus, to run it with the provided input files, the line is as following:

sbatch pipeline.sh p hapmap_subset
inputs/hapmap_subset_3cols.txt

If the script is able to run properly, the intermediates directory should have several files for each of the starting clusters.
To be specific, each cluster should have a .bed, .bim, .fam, .grm.bin, .grm.id, .grm.N.bin, .hh, .log, .txt, .eigenval and .eigenvec file.
Furthermore, in the working directory there should be three new files (in addition to the slurm file): all_individuals.txt,
reassigned_individuals.txt, and individuals_count_by_cluster.txt. 

all_individuals.txt should have 3-4 columns for each row. If the individual was reassigned to a new cluster, their new cluster assignment is
in the fourth column. The other columns are the FID, IID and their starting cluster.

reassigned_individuals.txt should have exactly 4 columns for each row, those being each individual’s FID, IID, starting cluster, and ending
cluster assignment. When we ran the script, we ended up with 43 reassigned individuals, only 6 of which were assigned to a cluster that is 
not the new “Other” cluster.

individuals_count_by_cluster.txt should have one row for each of the starting clusters, with one extra row for the new “Other” cluster. After
running the script, this file should look exactly as follows:

CEU     150

Other   37

MEX     89

CHB     129

YRI     186
