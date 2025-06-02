# -------------------------------------------------------------------------------------------------
# Remove all stored objects & unload non-base packages
rm(list = ls())

detachAllPackages <- function() {
  base_pkgs <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
  loaded_pkgs <- search()[grep("package:", search())]
  pkgs_to_detach <- setdiff(loaded_pkgs, paste("package:", base_pkgs, sep = ""))
  
  if (length(pkgs_to_detach) > 0) {
    lapply(pkgs_to_detach, function(pkg) detach(pkg, character.only = TRUE))
  }
}

detachAllPackages()

# -------------------------------------------------------------------------------------------------
# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(readxl)
library(dplyr)
library(vegan)
library(ape)
library(microbiome)
library(microbiomeutilities)
library(ggpubr)
library(data.table)
library(tidyverse)
library(BiocParallel)

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

# -------------------------------------------------------------------------------------------------
# Import data from Excel
otu_mat <- read_excel("colombia_metagenome.xlsx", sheet = "OTU_matrix")
tax_mat <- read_excel("colombia_metagenome.xlsx", sheet = "Taxonomy_table")
samples_df <- read_excel("colombia_metagenome.xlsx", sheet = "Sample_name")

# Define row names
otu_mat <- otu_mat %>% remove_rownames() %>% column_to_rownames(var="ASV")
tax_mat <- tax_mat %>% remove_rownames() %>% column_to_rownames(var="ASV")
row.names(samples_df) <- samples_df$Sample_name

# Convert to matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

# Create Phyloseq objects
ASV <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
SAM <- sample_data(samples_df)
rownames(SAM) <- SAM$Sample_name

# Merge Phyloseq objects
carbom <- phyloseq(ASV, TAX, SAM)

# -------------------------------------------------------------------------------------------------
# Subset Santander fermentation samples (removing unwanted groups)
unwanted_samples <- c("soil", "leaf", "pod", "hand", "tool", "fly", "fermentation_box")
carbom.sub <- subset_samples(carbom, !Sample_description3 %in% unwanted_samples)

# Filter to only include "Santander" location
carbom.sub <- subset_samples(carbom.sub, Location == "Santander")

# -------------------------------------------------------------------------------------------------
# Data filtering
carbom.sub.prun <- prune_samples(sample_sums(carbom.sub) >= 5000, carbom.sub)
carbom.sub.prun <- prune_taxa(taxa_sums(carbom.sub.prun) > 2, carbom.sub.prun)

# -------------------------------------------------------------------------------------------------
# Rarefaction (normalize sequencing depth)
carbom.rarefied <- rarefy_even_depth(
  carbom.sub.prun, 
  rngseed = 1, 
  sample.size = min(sample_sums(carbom.sub.prun)), 
  replace = FALSE
)

# -------------------------------------------- Select Bacteria Taxa -----------------------------------------------------
carbom.rarefied.sub = carbom.rarefied

# Select specific taxa (Bacteria)
carbom.rarefied.sub <- subset_taxa(carbom.rarefied.sub, Domain =="Bacteria")

## Get Bray-Curtis distances for Bacteria
ps1 = carbom.rarefied.sub
library(tidyr)

# Calculate Bray-Curtis distance matrix
set.seed(1)
bray_dist_matrix_bacteria <- phyloseq::distance(ps1, method = "bray") %>% 
  as.matrix() %>% 
  as_tibble(rownames = "A") %>%
  pivot_longer(-A, names_to="B", values_to="bacteria")


# -------------------------------------------- Select Fungi Taxa -----------------------------------------------------
carbom.rarefied.sub = carbom.rarefied

# Select specific taxa (Fungi)
carbom.rarefied.sub <- subset_taxa(carbom.rarefied.sub, Domain =="Fungi")

# Get Bray-Curtis distances for Fungi
ps1 = carbom.rarefied.sub

# Calculate Bray-Curtis distance matrix
set.seed(1)
bray_dist_matrix_fungi <- phyloseq::distance(ps1, method = "bray") %>% 
  as.matrix() %>% 
  as_tibble(rownames = "A") %>%
  pivot_longer(-A, names_to="B", values_to="fungi")


# ----------------------------------------------------------------------------------------------------------------------
# Merge distance data matrix
combined = inner_join(bray_dist_matrix_bacteria, bray_dist_matrix_fungi, by=c("A","B"))

# Select only values below the diagonal from the data matrix
combined_filter = combined %>% filter(A < B)

# ----------------------------------------------------------------------------------------------------------------------
## Mantel Test
library("tibble")

# Create lower triangle distance matrix
bacteria_mantel = combined %>% 
  select(A, B, bacteria) %>%
  pivot_wider(names_from = B, values_from = bacteria) %>%
  column_to_rownames("A") %>% 
  as.dist()

fungi_mantel = combined %>% 
  select(A, B, fungi) %>%
  pivot_wider(names_from = B, values_from = fungi) %>%
  column_to_rownames("A") %>% 
  as.dist()

# Run Mantel test
mantel(bacteria_mantel, fungi_mantel, method="pearson", permutations = 10000)

# ----------------------------------------------------------------------------------------------------------------------
# Plot
p = ggplot(combined_filter, aes(x=bacteria, y=fungi)) + 
  geom_point(shape = 21, fill = '#fe499a', size = 3, alpha = 0.8, stroke = 0.3, color = "black") +
  geom_smooth(method=lm, size=0.4, colour="#a51672", fill = "#dedde6")
A1 = p + theme_classic() + 
  ylab("Fungi Bray-Curtis Dissimilarity") + 
  xlab("Bacteria Bray-Curtis Dissimilarity") + 
  ggtitle("Santander") + 
  theme(legend.position="bottom") +
  theme(axis.text.y = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=10, color = "black"),
        axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 1, size=10, color = "black"),
        axis.title.y = element_text(size=10, face="bold"),
        axis.title.x = element_text(size=10, face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=10, face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  annotate("text", size=4, label = ~italic(p)~{"="}~"9.999e-05", x = 0.75, y = 0.15) +
  annotate("text", size=4, label = ~italic(r)~{"="}~0.5618, x = 0.75, y = 0.08) +
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A1
# Save the plot
ggsave("Fig1g.pdf", plot = A1, height=2.7, width=2.7, units='in')

