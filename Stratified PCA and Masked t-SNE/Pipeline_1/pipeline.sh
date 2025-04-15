#!/usr/local/bin/bash

#SBATCH --time=12:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=262144M   # memory per CPU core
#SBATCH --mail-user=nates99@byu.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#create a 1-dimensional associative array that maps each individual to their cluster assignment

declare -A individual_to_cluster

option=$1

if [ "$option" == "v" ]; then
	plink --bfile $2 --make-bed --out $2

elif [ "$option" == "p" ]; then
	prefix=$2
	bed="inputs/$prefix.bed"
	bim="inputs/$prefix.bim"
	fam="inputs/$prefix.fam"
fi

mapfile -t clusters < <(awk '{print $3}' "$3" | sort -u)

for cluster in "${clusters[@]}"; do
	echo "" > intermediates/$cluster.txt
done

while read -r FID IID cluster; do
       echo "$FID $IID" >> intermediates/$cluster.txt
done < $3

#stratified PLINK on each cluster, followed by stratified PCA

for cluster in "${clusters[@]}"; do
	plink --bfile inputs/$prefix --keep intermediates/$cluster.txt --make-bed --out intermediates/$cluster
	gcta --bfile intermediates/$cluster --make-grm --out intermediates/$cluster
	gcta --grm intermediates/$cluster --pca 10 --out intermediates/$cluster
done

#build a string with which to run regroup, then run regroup

regroup_args=()

for cluster in "${clusters[@]}"; do
	regroup_args+=("$cluster" "intermediates/${cluster}.eigenvec")
done

python3 regroup.py "${regroup_args[@]}"



