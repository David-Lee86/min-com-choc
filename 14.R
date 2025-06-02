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
library(grid)
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

# Select locations
carbom.sub=subset_samples(carbom.sub, Location %in% c("Huila", "Antioquia"))

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
# Merge at the Family level
carbom.rarefied.sub <- tax_glom(carbom.rarefied.sub, "Family")
carbom.rarefied.sub

# Filter method 1: Remove taxa with uncharacterized (NA) Family
carbom.rarefied.sub = subset_taxa(carbom.rarefied.sub, !is.na(Family) & !Family %in% c("", "uncharacterized"))
carbom.rarefied.sub

# Set ps1 as the filtered data
ps1 = carbom.rarefied.sub

# Edit the unclassified taxa
tax_table(ps1)[tax_table(ps1)[, "Family"] == "", "Family"] <- "NA"

# Remove the phy_tree
ps1@phy_tree <- NULL

# Aggregate rarefied data at the Family level with specific detection and prevalence thresholds
ps1_fam <- microbiome::aggregate_rare(ps1, "Family", detection = 5/100, prevalence = 35/100)
ps1_fam

# -------------------------------------------------------------------------------------------------
# Remove a specific taxa (Unknown Family)
ps1_fam <- subset_taxa(ps1_fam, Family != "Unknown")

# Reorder samples by Location and Time2
ps1_fam@sam_data$Location <- factor(ps1_fam@sam_data$Location, levels = c('Huila', 'Antioquia'))
ps1_fam@sam_data$Time2 <- factor(ps1_fam@sam_data$Time2, levels = c('0hr', '24hr', '48hr', '72hr', '96hr'))

#Rename samples
B <- c('0hr'='0h',
       '24hr'='24h',
       '48hr'='48h',
       '72hr'='72h',
       '96hr'='96h')

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

library(ggh4x)

# Create the plot (bar plot for Family taxa abundance across time)
F1 = ggplot(data = psmelt(ps1_fam), mapping = aes_string(x = "Time2", y = "Abundance")) +
  geom_bar(aes(fill = Family, color = Family), stat = "identity", position = "fill", width = 0.9, size = 0.01) +
  scale_fill_manual(values = mypal) + scale_color_manual(values = mypal) +
  theme_classic() + theme(legend.position = "right") +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.25, size = 8, color = "black"),
    axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 8, color = "black", face = "bold"),
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
    strip.background = element_blank(),
    panel.spacing = unit(0.1, "cm")
  ) +
  labs(x = "Sample", y = "Relative abundance") +
  facet_wrap(~Location, scales = "free_x", nrow = 1) +
  force_panelsizes(cols = c(5, 4)) +
  guides(fill = guide_legend(ncol = 2)) +
  guides(shape = "none") +
  theme(
    plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm"),
    strip.text = element_text(size = 8, face = "bold", margin = margin(t = 1, b = 1)),
    panel.spacing = unit(0.5, "cm"),
    strip.placement = "outside"
  ) +
  scale_x_discrete(expand = c(0, 0), label = B) +
  scale_y_continuous(expand = c(0, 0))
F1

# Remove the legend and adjust the plot
F1a = F1 + guides(color = 'none', fill = 'none', shape = 'none')
F1a

# Save the plot
ggsave("Fig2e.pdf", height = 1.5, width = 3, units = 'in', plot = F1a)


### extract legend option 1
library(ggpubr)
library(egg)
# Extract the legend. Returns a gtable
LEG <- get_legend(F1)
# Convert extracted legend to a ggplot object
LEG = as_ggplot(LEG)
LEG
ggsave("Fig2e_legend.pdf", height=5, width=6.5, units='in', plot = LEG)



