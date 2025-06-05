library(tidyverse)
library(Rtsne)

set.seed(100)
# Read data from a TSV file
#your_data <- as.data.frame(read.table("/Users/tannervanorden/Desktop/Bioinformatics Capstone/GCTA_Plink_AfricanSubset_Results_LD_Pruned.eigenvec", header = FALSE))
#EigenValues <- read.table("/Users/tannervanorden/Desktop/Bioinformatics Capstone/GCTA_Plink_AfricanSubset_Results_LD_Pruned.eigenval", header = FALSE, sep = '\t')
your_data <- as.data.frame(read.table("GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated.eigenvec", header = FALSE))
EigenValues <- read.table("GCTA_Plink_AfricanSubset_Results_LD_Pruned_Updated.eigenval", header = FALSE, sep = '\t')

eigenvalues = EigenValues$V1

#your_data <- your_data %>%
#  mutate(Location = c(rep("Utah", 165), rep("Han Chinese", 137), rep("Mexican", 86), rep("Yoruba", 203)))

#Location = c(rep("ACB", 96), rep("ASW", 61), rep("BEB", 86), rep("CDX", 93), rep("CEU", 99), rep("CHB", 103), rep("CHS", 105), rep("CLM", 94), rep("ESN", 99), rep("FIN", 99), rep("GBR", 91), rep("GIH", 103), rep("GWD", 113), rep("IBS", 107), rep("ITU", 102), rep("JPT", 104), rep("KHV", 99), rep("LWK", 99), rep("MSL", 85), rep("MXL", 64), rep("PEL", 85), rep("PJL", 96), rep("PUR", 104), rep("STU", 102), rep("TSI", 107), rep("YRI", 108))
Location = c(rep("African Caribbean in Barbados", 96), rep("African Ancestry in Southwest US", 61), rep("Esan in Nigeria", 99), rep("Western Gambian", 113), rep("Luhya in Webuye, Kenya", 99), rep("Mende in Sierra Leone", 85), rep("Yoruba in Ibadan, Nigeria", 89), "Remove", rep("Yoruba in Ibadan, Nigeria", 18)) #Second most outlier is point 78 biggest outlier is 90

# Calculate the sum of all eigenvalues
your_data$Location <- factor(Location)
your_data = your_data %>% slice(-643) %>% slice(-631)
Location = Location[c(-643, -631)]
SUM_TOTAL_EIGENVECS <- sum(eigenvalues)

# Calculate the percentage of variance explained by PC1 and PC2
variance_explained_pc1 <- (eigenvalues[1] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc2 <- (eigenvalues[2] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc3 <- (eigenvalues[3] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc4 <- (eigenvalues[4] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc5 <- (eigenvalues[5] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc6 <- (eigenvalues[6] / SUM_TOTAL_EIGENVECS) * 100
# Plot the PCA using ggplot2
pca <- ggplot(data = your_data, aes(x = V3, y = V4, color = Location)) + geom_point(alpha = 0.7, size = 3) +
  scale_color_brewer(palette = "Set1") + 
  labs(x = paste("PC1 (", round(variance_explained_pc1, 2), "%)"), y = paste("PC2 (", round(variance_explained_pc2, 2), "%)")) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.3), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = c(.95, 0.9),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))
ggsave("PCA_Not_Masked.png", plot = pca)

set.seed(100)
# Read data from a TSV file
#your_data <- as.data.frame(read.table("/Users/tannervanorden/Desktop/Bioinformatics Capstone/GCTA_Plink_AfricanSubset_Results_LD_Pruned.eigenvec", header = FALSE))
#EigenValues <- read.table("/Users/tannervanorden/Desktop/Bioinformatics Capstone/GCTA_Plink_AfricanSubset_Results_LD_Pruned.eigenval", header = FALSE, sep = '\t')
your_data <- as.data.frame(read.table("GCTA_Masked_Plink_Final.eigenvec", header = FALSE))
EigenValues <- read.table("GCTA_Masked_Plink_Final.eigenval", header = FALSE, sep = '\t')

eigenvalues = EigenValues$V1

#your_data <- your_data %>%
#  mutate(Location = c(rep("Utah", 165), rep("Han Chinese", 137), rep("Mexican", 86), rep("Yoruba", 203)))

#Location = c(rep("ACB", 96), rep("ASW", 61), rep("BEB", 86), rep("CDX", 93), rep("CEU", 99), rep("CHB", 103), rep("CHS", 105), rep("CLM", 94), rep("ESN", 99), rep("FIN", 99), rep("GBR", 91), rep("GIH", 103), rep("GWD", 113), rep("IBS", 107), rep("ITU", 102), rep("JPT", 104), rep("KHV", 99), rep("LWK", 99), rep("MSL", 85), rep("MXL", 64), rep("PEL", 85), rep("PJL", 96), rep("PUR", 104), rep("STU", 102), rep("TSI", 107), rep("YRI", 108))
Location = c(rep("African Caribbean in Barbados", 94), rep("African Ancestry in Southwest US", 61), rep("Esan in Nigeria", 88), rep("Western Gambian", 105), rep("Luhya in Webuye, Kenya", 99), rep("Mende in Sierra Leone", 79), rep("YRI", 100)) 


# Calculate the sum of all eigenvalues
your_data$Location <- factor(Location)
SUM_TOTAL_EIGENVECS <- sum(eigenvalues)

# Calculate the percentage of variance explained by PC1 and PC2
variance_explained_pc1 <- (eigenvalues[1] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc2 <- (eigenvalues[2] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc3 <- (eigenvalues[3] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc4 <- (eigenvalues[4] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc5 <- (eigenvalues[5] / SUM_TOTAL_EIGENVECS) * 100
variance_explained_pc6 <- (eigenvalues[6] / SUM_TOTAL_EIGENVECS) * 100
# Plot the PCA using ggplot2
pca <- ggplot(data = your_data, aes(x = V3, y = V4, color = Location)) + geom_point(alpha = 0.7, size = 3) +
  scale_color_brewer(palette = "Set1") + 
  labs(x = paste("PC1 (", round(variance_explained_pc1, 2), "%)"), y = paste("PC2 (", round(variance_explained_pc2, 2), "%)")) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.3), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = c(.95, 0.9),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))
ggsave("PCA_Masked.png", plot = pca)




eigenvec <- your_data[, -c(1,2, ncol(your_data))]

tsne_results <- Rtsne(eigenvec, theta = .5, perplexity = 50, exaggeration_factor = 10, 
                      max_iter = 15000, verbose = TRUE)

data = as.data.frame(tsne_results$Y)


T_SNE <- ggplot(data = data, aes(x = V1, y = V2, color = Location)) + geom_point(alpha = 0.7, size = 4) +
  scale_color_brewer(palette = "Set1") + scale_x_continuous(limits = c(-10, 20)) +
  labs(x = "t-SNE 1", y = "t-SNE 2 ") +
  theme_minimal() + 
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.3), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 17),
        plot.title = element_text(hjust = 0.5, size = 19),
        legend.position = c(.4, .3),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 17))

ggsave("T_SNE_Masked.png", plot = T_SNE)
T_SNE
