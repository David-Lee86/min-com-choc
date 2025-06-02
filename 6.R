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

carbom.sub <- subset_taxa(carbom.sub, Domain == "Fungi")

# -------------------------------------------------------------------------------------------------
# Data filtering
carbom.sub.prun <- prune_samples(sample_sums(carbom.sub) >= 100, carbom.sub)
carbom.sub.prun <- prune_taxa(taxa_sums(carbom.sub.prun) > 2, carbom.sub.prun)


# -------------------------------------------------------------------------------------------------
# Rarefaction (normalize sequencing depth)
carbom.rarefied <- rarefy_even_depth(
  carbom.sub.prun, 
  rngseed = 1, 
  sample.size = min(sample_sums(carbom.sub.prun)), 
  replace = FALSE
)

# -------------------------------------------------------------------------------------------------
# Alpha Diversity Calculation
rich = estimate_richness(carbom.rarefied, measures="Shannon")

## create a column with sample name
rich <- tibble::rownames_to_column(rich, "Sample_name")

## get metadata
metadata = samples_df %>% as.data.frame

# merge two data frames by ID
mydata <- merge(rich,metadata, by="Sample_name")

#Save reorder Alpha diversity table
write.csv(mydata, file = "alpha_diversity_table.csv")


# -----------------------------------------------------------------------------------
# Clean Environment and Detach Non-Base Packages
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

# -----------------------------------------------------------------------------------
# Load Required Libraries
library(agricolae)
library(tidyverse)
library(BiocParallel)

# -----------------------------------------------------------------------------------
# Set Working Directory (RStudio Only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register Parallel Processing with 6 Cores
register(SnowParam(6))

# -----------------------------------------------------------------------------------
# Import Data
data_set <- read.csv("alpha_diversity_table.csv")

# Subset Data (Select Shannon)
data_set <- data_set %>% select(Sample_name, Sample_description3, Time, Replicate, Shannon) %>%
  rename(value = Shannon)

# Define Factor Levels for Sample Description
data_set$Sample_description3 <- factor(
  data_set$Sample_description3, 
  levels = c('soil', 'leaf', 'pod', 'hand', 'tool', 'fly', 'fermentation_box',
             'fermentation_0hrs', 'fermentation_24hrs', 'fermentation_48hrs', 
             'fermentation_72hrs', 'fermentation_96hrs')
)

# Get Maximum Values for Label Placement
abs_max <- max(data_set$value)

maxs <- data_set %>%
  group_by(Sample_description3) %>%
  summarise(value = min(max(value), 
                        max(value[value <= (quantile(value, 0.75) + 1.5 * 
                                              diff(quantile(value, c(0.25, 0.75))))])) + 0.1 * abs_max)

# -----------------------------------------------------------------------------------
# ANOVA and Tukey HSD Test
Tukey_test <- aov(value ~ Sample_description3 + Replicate, data = data_set) %>%
  HSD.test("Sample_description3", group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "Sample_description3") %>%
  rename(Letters_Tukey = groups) %>%
  select(-value) %>%
  left_join(maxs, by = "Sample_description3")

# Define Colors and Labels
B <- c('fermentation_0hrs'='0',
       'fermentation_24hrs'='24',
       'fermentation_48hrs'='48',
       'fermentation_72hrs'='72',
       'fermentation_96hrs'='96')

mypal <- c('#e3c9a8','#baa692','#4a8182','#71aca8','#d4e4e3')
mypal <- c('#2d4800', '#528501', '#a9d11c', '#eaabcd', '#c0047b')


# -----------------------------------------------------------------------------------
# Plot Shannon Diversity Index
p <- ggplot(data_set, aes(x = Sample_description3, y = value)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2, aes(fill = Sample_description3), 
               colour = "black", alpha = 0.7) +
  geom_text(data = Tukey_test, aes(label = Letters_Tukey), size = 6, colour = "black") +
  scale_fill_manual(values = mypal)

# Extract boxplot statistics
dat <- ggplot_build(p)$data[[1]]

# Final plot with additional formatting
p <- p + 
  geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), 
               colour = "black", size = 0.5) +
  theme_classic() + 
  labs(y = "Shannon Diversity Index", x = "Time (h)") +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", vjust = 1),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    axis.ticks = element_line(colour = "black", size = 0.35),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.35),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_blank()
  ) +
  geom_point(aes(fill = Sample_description3), size = 1.25, shape = 16, alpha = 0.8, 
             colour = "black", stroke = 0.1, position = position_jitterdodge(3)) +
  scale_x_discrete(label = B)
p
# Save figure
ggsave("Fig1d.pdf", height = 3, width = 3.2, units = "in", plot = p)


