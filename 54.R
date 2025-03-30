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
# Subset Family Data
ps1.com <- carbom.rarefied.sub
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "NA"  # Edit unclassified taxa
ps1.com@phy_tree <- NULL  # Remove phy_tree

# Aggregate by Family level with specific detection and prevalence thresholds
ps1.com.fam <- microbiome::aggregate_rare(ps1.com, "Family", detection = 10/100, prevalence = 55/100)

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
mypal=c("Acetobacteraceae"="#771155",
        "Alcaligenaceae"="#AA4488",
        "Bacillaceae"="#CC99BB",
        "Burkholderiaceae"="#114477",
        "Clostridiaceae"="#4477AA",
        "Comamonadaceae"="#77AADD",
        "Enterobacteriaceae"="#117777",
        "Enterococcaceae"="#44AAAA",
        "Erwiniaceae"="#77CCCC",
        "Flavobacteriaceae"="#117744",
        "Lactobacillaceae"="#44AA77",
        "Leptolyngbyaceae"="#88CCAA",
        "Leuconostocaceae"="#777711",
        "Moraxellaceae"="#AAAA44",
        "Morganellaceae"="#DDDD77",
        "Other"="#774411",
        "Paenibacillaceae"="#AA7744",
        "Pasteurellaceae"="#DDAA77",
        "Pectobacteriaceae"="#B26F00",
        "Pseudomonadaceae"="#FFB232",
        "Rhizobiaceae"="#FFCF7F",
        "Streptococcaceae"="#771122",
        "Streptomycetaceae"="#AA4455",
        "Vibrionaceae"="#DD7788",
        "Xanthomonadaceae"="black",
        "Yersiniaceae"="#3F3F3F",
        'Rhodanobacteraceae'="#ff2e31",
        'Sphingomonadaceae'='#da2c47',
        'Staphylococcaceae'='#dd2c6e',
        'Synechococcaceae'="#e28f9f",
        "Bradyrhizobiaceae"="#fe3b5b",
        "Chitinophagaceae"="#b9db45",
        "Halomonadaceae"="#00c8ec",
        "Microbacteriaceae"="#3e4f94",
        "Mycobacteriaceae"="#5377c1",
        "Oxalobacteraceae"="#fd70a0",
        "Rhodobacteraceae"="#fdbd8c",
        "Rhodospirillaceae"="#465b54",
        "Sphingobacteriaceae"="#80aa9c",
        "Cytophagaceae"="#26ff00",
        "Hymenobacteraceae"="#aa00ff",
        "Hyphomicrobiaceae"="#00aaff",
        "Phyllobacteriaceae"="#e46f12"
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
ggsave("Fig1e.pdf", height = 2.75, width = 2.8, units = 'in', plot = F1a)

