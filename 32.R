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
library(knitr)
library(rmdformats)
library(vegan)
library(ape)
library(DESeq2)
library(microbiome)
library(microbiomeutilities)
library(ggpubr)
library(data.table)
library(ggtree)
library(gridExtra)
library(grid)
library(scales)
library(hrbrthemes)
library(tidyverse)
library(BiocParallel)

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))
# -------------------------------------------------------------------------------------------------
# Import data from excel file
otu_mat<- read_excel("SYNCOM_validation_ITS.xlsx", sheet = "OTU_matrix")
tax_mat<- read_excel("SYNCOM_validation_ITS.xlsx", sheet = "Taxonomy_table")
samples_df <- read_excel("SYNCOM_validation_ITS.xlsx", sheet = "Sample_name")
# -------------------------------------------------------------------------------------------------
## format tax_mat object
# Function to process taxonomy names
process_tax_names <- function(tax_name) {
  ifelse(grepl("__", tax_name), sub("^.+__", "", tax_name), NA)
}

# Apply the function to all columns in tax_mat except ASV using dplyr
tax_mat <- tax_mat %>%
  mutate(across(-ASV, ~ process_tax_names(.)))

# -------------------------------------------------------------------------------------------------
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
# Data filtering
carbom.sub.prun = prune_samples(sample_sums(carbom)>=500, carbom)
carbom.sub.prun

#remove a specific taxa
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Kingdom !="Viridiplantae")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Kingdom !="NA")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Kingdom !="Metazoa")
carbom.sub.prun

carbom.sub.prun <- prune_taxa(taxa_sums(carbom.sub.prun) > 2, carbom.sub.prun)
carbom.sub.prun = prune_samples(sample_sums(carbom.sub.prun)>=1000, carbom.sub.prun)
carbom.sub.prun

# -------------------------------------------------------------------------------------------------
# Rarefaction (normalize sequencing depth)
carbom.rarefied <- rarefy_even_depth(
  carbom.sub.prun, 
  rngseed = 1, 
  sample.size = min(sample_sums(carbom.sub.prun)), 
  replace = FALSE
)

# -------------------------------------------------------------------------------------------------
#Subset Family
ps1 = carbom.rarefied
tax_table(ps1)[tax_table(ps1)[, "Family"] == "", "Family"] <- "NA"# edit the unclassified taxa
ps1@phy_tree <- NULL# remove the phy_tree
ps1.fam <- microbiome::aggregate_rare(ps1, "Family", detection = 1/100, prevalence = 15/100)
ps1.fam

#remove a specific taxa
ps1.fam <- subset_taxa(ps1.fam, Family !="Unknown")

mypal <- c(
  "Aspergillaceae" = "#587c7c",
  "Debaryomycetaceae" = "#003f5e",
  "Dipodascaceae" = "#007c84",
  "Glomerellaceae" = "#e8d666",
  "Nectriaceae" = "#9ea615",
  "Other" = "#f4f4f4",
  "Phaffomycetaceae" = "#bbe0ce",
  "Pichiaceae" = "#fddde6",
  "Pyriculariaceae" = "#f7a08c",
  "Saccharomycetaceae" = "#f15640",
  "Trichomonascaceae" = "#1c1c1c",
  "Ustilaginaceae" = "#e5e1e0",
  "Chaetomiaceae" = "#7c5ea8",
  "Mycosphaerellaceae" = "#078eb6",
  "Sordariaceae" = "#83cdbe",
  "Clavicipitaceae" = "#934171",
  "Malasseziaceae" = "#592852",
  "Schizosaccharomycetaceae" = "#927555",
  "Cryptococcaceae" = "#ee217c",
  "Sclerotiniaceae" = "#36804c",
  "Unikaryonidae" = "#aeff00",
  "Pleosporaceae" = "#aa00ff",
  "Rhynchogastremataceae" = "#00aaff",
  "Saccharomycetales_fam_Incertae_sedis" = "#c9de57",
  "Saccharomycodaceae" = "#9675b4",
  "Saccharomycopsidaceae" = "#555555",  # different gray from Trichomonascaceae
  "Ceratocystidaceae" = "#ffb347",  # new warm orange
  "Cladosporiaceae" = "#749c2e",    # olive green
  "Hypocreales_fam_Incertae_sedis" = "#c5c134",  # yellow-olive
  "Marasmiaceae" = "#80ced6",       # teal
  "Mucoraceae" = "grey",         # soft pink
  "Neodevriesiaceae" = "#0d6f76",   # dark teal
  "Neopyrenochaetaceae" = "#fe3b5b"
)

#Reorder samples
ps1.fam@sam_data$Description2<- factor(ps1.fam@sam_data$Description2, levels = c("0hr","NMC_48hr","NMC_96hr",
                                                                                     "MC_48hr", "MC_96hr"))
#Rename samples
B=c("0hr"='0h',
    "NMC_48hr"='No SYNCOM(48h)',
    "NMC_96hr"='No SYNCOM(96h)',
    "MC_48hr"='SYNCOM(48h)',
    "MC_96hr"='SYNCOM(96h)')

# Create a named vector for the legend labels
new_labels <- c(
"Hypocreales_fam_Incertae_sedis"="Hypocreales fam Incertae sedis",
  "Saccharomycetales_fam_Incertae_sedis"="Saccharomycetales fam Incertae sedis"
)

#plot1 with legend
A1 = ggplot(data = psmelt(ps1.fam), mapping = aes_string(x = "Description2",y = "Abundance")) +
  geom_bar(aes(fill=Family, color=Family), stat="identity", position="fill", width = 0.9, size=0.01) +
  theme_minimal() + theme(legend.position = "right") +
  scale_fill_manual(values = mypal, labels = new_labels) +
  scale_color_manual(values = mypal, labels = new_labels) +
  theme(axis.text.y   = element_text(angle=0, hjust = 0, vjust = 0.5, size=8, color = "black"),
        #axis.text.x   = element_text(angle=90, hjust = 1, vjust = 0.5, size=8, color = "black"),
        axis.text.x   = element_text(angle=35, hjust = 1, vjust = 1, size=8, color = "black", 
                                     margin = margin(t = 0, b = 0)), # Added margin
        axis.title.y  = element_text(angle=90, hjust = 0.5, vjust = 0.5, size=8, color = "black", face="bold"),
        axis.title.x  = element_blank(),
        plot.title = element_blank(),
        legend.title = element_text(size = 11, face="bold"),
        legend.key.size = unit(0.2, "cm"),
        legend.text=element_text(size=7, face=c("plain", "italic")),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.line = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(size=8, face="bold"),
        panel.spacing = unit(0.1, "cm")) +
  labs(x="Sample", y="Relative abundance") +
  guides(fill=guide_legend(ncol=2)) +
  guides(shape="none") +
  theme(plot.margin=unit(c(0.1,0,0.1,0.1),"cm")) +
  labs(title="ps1.fam") +
  scale_x_discrete(label = B, expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
A1
A1a = A1 + guides(color = "none", shape="none", fill="none")
A1a
ggsave("Fig5f.pdf", height=2.4, width=2, units='in', plot = A1a)


### extract legend option 1
library(ggpubr)
library(egg)
# Extract the legend. Returns a gtable
A1_leg <- get_legend(A1)
# Convert extracted legend to a ggplot object
LEG = as_ggplot(A1_leg)
LEG
ggsave("Fig5f_legend.pdf", height=5, width=6.5, units='in', plot = LEG)


############################################################################################################################
taxa_list = get_taxa_unique(ps1.fam, "Family")
cat(taxa_list)

