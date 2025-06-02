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
carbom.sub.prun <- prune_samples(sample_sums(carbom.sub) >= 1000, carbom.sub)
carbom.sub.prun <- prune_taxa(taxa_sums(carbom.sub.prun) > 2, carbom.sub.prun)

# -------------------------------------------------------------------------------------------------
# Rarefaction (normalize sequencing depth)
carbom.rarefied <- rarefy_even_depth(
  carbom.sub.prun, 
  rngseed = 1, 
  sample.size = min(sample_sums(carbom.sub.prun)), 
  replace = FALSE
)

carbom.rarefied.sub <- subset_taxa(carbom.rarefied, Domain == "Bacteria")

# -------------------------------------------------------------------------------------------------
# PERMANOVA Analysis
ps1 <- carbom.rarefied.sub

# Calculate Beta Diversity
set.seed(1)
bray.pca.dist <- phyloseq::distance(ps1, method = "bray")  # Bray-Curtis distance matrix
metadf <- data.frame(sample_data(ps1))  # Convert sample data to dataframe

# PERMANOVA analysis: test main effects with interactions
adonis2(
  bray.pca.dist ~ Time2 * Harvest + Replicate, 
  data = metadf, 
  permutations = 9999
)  # Adonis test

# -------------------------------------------------------------------------------------------------
# Bray-Curtis Distances with CAP Ordination
# Effect of Time while constraining Harvest
ps1@sam_data$Time2 <- factor(ps1@sam_data$Time2, 
                             levels = c("0hr", "24hr", "48hr", "72hr", "96hr"))

# Define Colors and Labels
B <- c('0hr'='0h',
       '24hr'='24h',
       '48hr'='48h',
       '72hr'='72h',
       '96hr'='96h')

# Define Color Palette
mypal <- c('#2d4800', '#528501', '#a9d11c', '#eaabcd', '#c0047b')

ordcap <- ordinate(ps1, "CAP", "bray", ~ Time2 + Condition(Harvest + Replicate))
distBC_CAP.plot <- plot_ordination(ps1, ordcap, "Sample_name", color = "Time2")
distBC_CAP.plot$layers <- distBC_CAP.plot$layers[-1]  # Remove default layer

# Customize Plot
final_plot <- distBC_CAP.plot + 
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_point(aes(fill = Time2, color = Time2), shape = 21, size = 7, 
             color = "black", stroke = 0.6) +
  scale_fill_manual(values = mypal, label = B) +
  theme_classic() +
  guides(color = "none", shape = "none", fill = guide_legend(title = "Time")) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.07, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.05),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.35),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    legend.position = "bottom"
  ) +
  annotate("text", size = 4.5, label = ~italic(p)~"<"~"0.0001", x = -0.5, y = 1.4) +
  annotate("text", size = 4.5, label = ~italic(R)^"2"~"="~0.76315, x = -0.5, y = 1.7) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 6)))
final_plot

# Save Plot
ggsave("Fig1c_ii.pdf", plot = final_plot, height = 3.8, width = 3.8, units = "in")

