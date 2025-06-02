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

#-------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

# -------------------------------------------------------------------------------------------------
# Import data from excel file
otu_mat<- read_excel("SYNCOM_validation_16S.xlsx", sheet = "OTU_matrix")
tax_mat<- read_excel("SYNCOM_validation_16S.xlsx", sheet = "Taxonomy_table")
samples_df <- read_excel("SYNCOM_validation_16S.xlsx", sheet = "Sample_name")

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
carbom
# -------------------------------------------------------------------------------------------------
# Data filtering
carbom.sub = carbom

#number of reads for each sample
read.count=sample_sums(carbom.sub)
read.count
sample_reads_list = read.count %>% as.data.frame

#remove samples with less than 1000 total reads
carbom.sub.prun = prune_samples(sample_sums(carbom.sub)>=1000, carbom.sub)
carbom.sub.prun

#remove a specific taxa
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Order !="Chloroplast")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Family !="Chloroplast")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Family !="Mitochondria")
carbom.sub.prun

#number of singletons
length(which(taxa_sums(carbom.sub.prun) <= 1))
#Remove singletons (taxa with 0 and 1 abundance)
carbom.sub.prun <- prune_taxa(taxa_sums(carbom.sub.prun) > 2, carbom.sub.prun)
#number of singletons
length(which(taxa_sums(carbom.sub.prun) <= 1))
carbom.sub.prun
read.count=sample_sums(carbom.sub.prun)
read.count

#remove samples with less than 1000 total reads
carbom.sub.prun
carbom.sub.prun = prune_samples(sample_sums(carbom.sub.prun)>=500, carbom.sub.prun)
carbom.sub.prun
read.count=sample_sums(carbom.sub.prun)
read.count

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
ps1.fam <- microbiome::aggregate_rare(ps1, "Family", detection = 1/100, prevalence = 47/100)
ps1.fam

#remove a specific taxa
ps1.fam <- subset_taxa(ps1.fam, Family !="Unknown")

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
        "Ilumatobacteraceae"="#26ff00",
        "Hymenobacteraceae"="#aa00ff",
        "Nocardiaceae"="#aa00ff",
        "Hyphomicrobiaceae"="#00aaff",
        "Propionibacteriaceae"="#00aaff",
        "Phyllobacteriaceae"="#e46f12"
)

#Reorder samples
ps1.fam@sam_data$Treatment<- factor(ps1.fam@sam_data$Treatment, levels = c("NMC","MC-F","MC"))
ps1.fam@sam_data$Description2<- factor(ps1.fam@sam_data$Description2, levels = c("0hr","NMC_48hr","NMC_96hr", "MC_48hr", "MC_96hr"))

#Rename samples
B=c("0hr"='0h',
    "NMC_48hr"='No SYNCOM(48h)',
    "NMC_96hr"='No SYNCOM(96h)',
    "MC_48hr"='SYNCOM(48h)',
    "MC_96hr"='SYNCOM(96h)')

#plot1 with legend
A1 = ggplot(data = psmelt(ps1.fam), mapping = aes_string(x = "Description2",y = "Abundance")) +
  geom_bar(aes(fill=Family, color=Family), stat="identity", position="fill", width = 0.9, size=0.01) +
  theme_minimal() + theme(legend.position = "right") +
  scale_fill_manual(values = mypal) + scale_color_manual(values = mypal) +
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
  #scale_x_discrete(limits = sample.order) +#label = B
  #facet_wrap(~Treatment, scales="free_x", nrow=1) +
  guides(fill=guide_legend(ncol=2)) +
  guides(shape="none") +
  theme(plot.margin=unit(c(0.1,0,0.1,0.1),"cm")) +
  labs(title="ps1.fam") +
scale_x_discrete(label = B, expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
A1
A1a = A1 + guides(color = "none", shape="none", fill="none")
A1a
ggsave("Fig5e.pdf", height=2.4, width=2, units='in', plot = A1a)

### extract legend option 1
library(ggpubr)
library(egg)
# Extract the legend. Returns a gtable
A1_leg <- get_legend(A1)
# Convert extracted legend to a ggplot object
LEG = as_ggplot(A1_leg)
LEG
ggsave("Fig5e_legend.pdf", height=5, width=6.5, units='in', plot = LEG)

############################################################################################################################
taxa_list = get_taxa_unique(ps1.fam, "Family")
cat(taxa_list)
