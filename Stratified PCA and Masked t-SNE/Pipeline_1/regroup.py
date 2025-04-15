import sys
import pandas as pd
import numpy as np
from collections import Counter

eigenvec_names = sys.argv[1::2]
eigenvecs = sys.argv[2::2]

cluster_to_means = {}
cluster_to_stds = {}
cluster_to_outlier_indices = {}
cluster_to_outlier_ids = {}
all_individuals = []

for index, eigenvec in enumerate(eigenvecs):

    current_name = eigenvec_names[index]

    eigenvec_df = pd.read_csv(eigenvec, sep='\s+', header=None)

    pcs = eigenvec_df.iloc[:, 2:].to_numpy()

    pc_means = np.mean(pcs, axis=0)
    pc_stds = np.std(pcs, axis=0)

    cluster_to_means[current_name] = pc_means
    cluster_to_stds[current_name] = pc_stds

    pc1_pass = np.abs(pcs[:, 0] - pc_means[0]) < 2 * pc_stds[0]
    pc2_pass = np.abs(pcs[:, 1] - pc_means[1]) < 3 * pc_stds[1]
    pc3_pass = np.abs(pcs[:, 2] - pc_means[2]) < 3 * pc_stds[2]

    outliers = ~(pc1_pass & pc2_pass & pc3_pass)

    outlier_indices = np.where(outliers)[0]

    outlier_ids = eigenvec_df.iloc[outlier_indices, :2]

    cluster_to_outlier_indices[current_name] = outlier_indices
    cluster_to_outlier_ids[current_name] = outlier_ids
    
    for i, row in eigenvec_df.iterrows():
        fid, iid = row[0], row[1]
        all_individuals.append((fid, iid, current_name))


#    outlier_ids.to_csv(eigenvec_names[index] + "_outliers1.txt", sep="\t", index=False, header=False) #this line isn't entirely necessary, just to make sure what I was doing was working to this point

updated_means = {}
updated_stds = {}

for index, eigenvec in enumerate(eigenvecs):
    current_name = eigenvec_names[index]
    eigenvec_df = pd.read_csv(eigenvec, sep='\s+', header=None)
    
    non_outlier_df = eigenvec_df.drop(cluster_to_outlier_indices[current_name])
        
    
    non_outliers_pcs = non_outlier_df.iloc[:, 2:].to_numpy()
    non_outliers_means = np.mean(non_outliers_pcs, axis=0)
    non_outliers_stds = np.std(non_outliers_pcs, axis=0)
    
    updated_means[current_name] = non_outliers_means
    updated_stds[current_name] = non_outliers_stds

#Now for each outlier, try each of the other clusters and see if it fits into that cluster based on the criteria
reassigned = {}

for index, name in enumerate(eigenvec_names):
    current_outlier_ids = cluster_to_outlier_ids[name]
    current_outlier_indices = cluster_to_outlier_indices[name]
    
    if len(current_outlier_ids) == 0:
        continue 

    current_eigenvec_df = pd.read_csv(eigenvecs[index], sep='\s+', header=None)
    current_outlier_pcs = current_eigenvec_df.iloc[current_outlier_indices, 2:].to_numpy()

    for i, outlier_pc in enumerate(current_outlier_pcs):
        individual_id = current_outlier_ids.iloc[i, 1]
        family_id = current_outlier_ids.iloc[i, 0]
        found_match = False
        
        for new_cluster in eigenvec_names:
            if name != new_cluster:
                current_means = updated_means[new_cluster]
                current_stds = updated_stds[new_cluster]

                pc1_pass = np.abs(outlier_pc[0] - current_means[0]) < 2 * current_stds[0]
                pc2_pass = np.abs(outlier_pc[1] - current_means[1]) < 3 * current_stds[1]
                pc3_pass = np.abs(outlier_pc[2] - current_means[2]) < 3 * current_stds[2]

                if pc1_pass and pc2_pass and pc3_pass:
                    print(f"Outlier {current_outlier_ids.iloc[i, 1]} fits in cluster {new_cluster}")
                    reassigned[(family_id, individual_id)] = new_cluster
                    found_match = True
                    break
        if not found_match:
            reassigned[(family_id, individual_id)] = "Other"
            print(f"Outlier {current_outlier_ids.iloc[i, 1]} does not fit in any cluster")

all_individuals_output = "all_individuals.txt"
reassigned_individuals = "reassigned_individuals.txt"
individuals_count_by_clusters = "individuals_count_by_cluster.txt"
original_clusters = {}

for i in range(len(all_individuals)):
    fid, iid, original_cluster = all_individuals[i]
    new_cluster = reassigned.get((fid, iid), original_cluster)
    #if new_cluster != original_cluster:
        #print(f"Individual {iid} was moved from {original_cluster} to {new_cluster}")
    original_clusters[(fid, iid)] = original_cluster
    all_individuals[i] = (fid, iid, new_cluster)

with open(all_individuals_output, "w") as Output:
    for fid, iid, new_cluster in all_individuals:
        original_cluster = original_clusters[(fid, iid)]
        if original_cluster != new_cluster:
            Output.write(f"{fid}\t{iid}\t{original_cluster}\t{new_cluster}\n")
        else:
            Output.write(f"{fid}\t{iid}\t{original_cluster}\n")
#print(f"Moved individuals are in {reassigned_individuals}")

with open(reassigned_individuals, "w") as Output:
    for fid, iid, new_cluster in all_individuals:
        original_cluster = original_clusters[(fid, iid)]
        if new_cluster != original_cluster:
            Output.write(f"{fid}\t{iid}\t{original_cluster}\t{new_cluster}\n")
#print(f"Reassigned individuals saved to {all_individuals_output}")

with open(individuals_count_by_clusters, "w") as Output:
    cluster_count = {}
    for _, _, cluster in all_individuals:
        if cluster in cluster_count:
            cluster_count[cluster] += 1
        else:
            cluster_count[cluster] = 1
    for cluster, count in cluster_count.items():
        Output.write(f"{cluster}\t{count}\n")
#print(f"Cluster count is save")
