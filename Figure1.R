library(adegenet)
library(tibble)
library(plotly)
library(ggplot2)
library(Rtsne)
library(viridis)


set.seed(100)
# Read data from a TSV file
your_data <- as.data.frame(read.table("/Users/benecharry/Desktop/BioCapstone/GCTA_Plink_AfricanSubset_Results.eigenvec", header = FALSE))
EigenValues <- read.table("/Users/benecharry/Desktop/BioCapstone/GCTA_Plink_AfricanSubset_Results.eigenval", header = FALSE, sep = '\t')

eigenvalues = EigenValues$V1
kmeans_result <- kmeans(your_data[, c("V3")], centers = 3)

#your_data <- your_data %>%
#  mutate(Location = c(rep("Utah", 165), rep("Han Chinese", 137), rep("Mexican", 86), rep("Yoruba", 203)))

#Location = c(rep("ACB", 96), rep("ASW", 61), rep("BEB", 86), rep("CDX", 93), rep("CEU", 99), rep("CHB", 103), rep("CHS", 105), rep("CLM", 94), rep("ESN", 99), rep("FIN", 99), rep("GBR", 91), rep("GIH", 103), rep("GWD", 113), rep("IBS", 107), rep("ITU", 102), rep("JPT", 104), rep("KHV", 99), rep("LWK", 99), rep("MSL", 85), rep("MXL", 64), rep("PEL", 85), rep("PJL", 96), rep("PUR", 104), rep("STU", 102), rep("TSI", 107), rep("YRI", 108))
Location = c(rep("African Caribbean in Barbados", 96), rep("African Ancestry in Southwest US", 61), rep("Esan in Nigeria", 99), rep("Western Gambian", 113), rep("Luhya in Webuye, Kenya", 99), rep("Mende in Sierra Leone", 85), rep("Yoruba in Ibadan, Nigeria", 108))
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
pca <- ggplot(data = your_data, aes(x = V3, y = V4, color = Location)) + 
  geom_jitter(size = 3, alpha = 0.7, width = 0.02, height = 0.02) +
  scale_color_viridis_d(option = "turbo") +
  labs(x = paste("PC1 (", round(variance_explained_pc1, 2), "%)"), y = paste("PC2 (", round(variance_explained_pc2, 2), "%)"), title = "PCA Plot of African and African Diaspora Populations from the 1000 Genomes") +
  theme_minimal() + 
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.3), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14))
pca