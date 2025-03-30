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

carbom.rarefied.sub <- subset_taxa(carbom.rarefied, Domain == "Fungi")

# -------------------------------------------------------------------------------------------------
# Subset Family Data
ps1.com <- carbom.rarefied.sub
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "NA"  # Edit unclassified taxa
ps1.com@phy_tree <- NULL  # Remove phy_tree

# Aggregate by Family level with specific detection and prevalence thresholds
ps1.com.fam <- microbiome::aggregate_rare(ps1.com, "Family", detection = 10/100, prevalence = 25/100)

Santander = ps1.com.fam

#remove a specific taxa
Santander <- subset_taxa(Santander, Family !="Unknown")

#Reorder samples
sample.order=c('fermentation_0hrs','fermentation_24hrs','fermentation_48hrs','fermentation_72hrs','fermentation_96hrs')

#Rename samples
B <- c('fermentation_0hrs'='0h',
       'fermentation_24hrs'='24h',
       'fermentation_48hrs'='48h',
       'fermentation_72hrs'='72h',
       'fermentation_96hrs'='96h')

####santander colouring
mypal=c("Aspergillaceae"="#587c7c",
        "Debaryomycetaceae"="#003f5e",
        "Dipodascaceae"="#007c84",
        "Glomerellaceae"="#e8d666",
        "Nectriaceae"="#9ea615",
        "Other"="#f4f4f4",
        "Phaffomycetaceae"="#bbe0ce",
        "Pichiaceae"="#fedcc1",
        "Pyriculariaceae"="#f7a08c",
        "Saccharomycetaceae"="#f15640",
        "Trichomonascaceae"="black",
        "Ustilaginaceae"="#e5e1e0",
        "Chaetomiaceae"="#7c5ea8",
        "Mycosphaerellaceae"="#078eb6",
        "Sordariaceae"="#83cdbe",
        "Clavicipitaceae"="#934171",
        "Malasseziaceae"="#592852",
        "Schizosaccharomycetaceae"="#927555",
        "Cryptococcaceae"="#ee217c",
        "Sclerotiniaceae"="#36804c",
        "Unikaryonidae"="#aeff00"
)

# Create the bar plot
F1 <- ggplot(data = psmelt(Santander), aes(x = Sample_description3, y = Abundance)) +
  geom_bar(aes(fill = Family, color = Family), stat = "identity", position = "fill", width = 0.9, size = 0.01) +
  theme_minimal() +
  scale_fill_manual(values = mypal) +
  scale_color_manual(values = mypal) +
  
  # Customize the plot theme
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.25, size = 11, color = "black", margin = margin(t = -2)),
    axis.title.y = element_text(size = 11, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 7, face = c("plain", "italic")),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.line = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.1, "cm"),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  
  # Customize axis labels and legends
  labs(x = "Sample", y = "Relative abundance") +
  scale_x_discrete(limits = sample.order, labels = B) +
  guides(fill = guide_legend(ncol = 2), shape = "none")


# Save plot without legend
F1a <- F1 + guides(color = 'none', fill = 'none', shape = 'none')
F1
ggsave("Fig1f.pdf", height = 2.75, width = 2.8, units = 'in', plot = F1a)

