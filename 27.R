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

#remove samples with less than 1000 total reads
carbom.sub.prun = prune_samples(sample_sums(carbom.sub)>=1000, carbom.sub)

#remove a specific taxa
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Order !="Chloroplast")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Family !="Chloroplast")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Family !="Mitochondria")
carbom.sub.prun

#Remove singletons (taxa with 0 and 1 abundance)
carbom.sub.prun <- prune_taxa(taxa_sums(carbom.sub.prun) > 2, carbom.sub.prun)
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
# Alpha Diversity Calculation
rich = estimate_richness(carbom.rarefied)

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
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))
# -----------------------------------------------------------------------------------
#import data
data_set = read.csv("alpha_diversity_table.csv")

# -----------------------------------------------------------------------------------
# Subset Data (Select Shannon)
data_set <- data_set %>% select(Sample_name, Description2, Rep, Shannon) %>%
  rename(value = Shannon)

# Define Factor Levels for Sample Description
data_set$Description2<- factor(
  data_set$Description2,
  levels = c("0hr","NMC_48hr","NMC_96hr", "MC_48hr", "MC_96hr")
  )

# Get Maximum Values for Label Placement
abs_max <- max(data_set$value)

maxs <- data_set %>%
  group_by(Description2) %>%
  summarise(value = min(max(value), 
                        max(value[value <= (quantile(value, 0.75) + 1.5 * 
                                              diff(quantile(value, c(0.25, 0.75))))])) + 0.1 * abs_max)

# -----------------------------------------------------------------------------------
# ANOVA and Tukey HSD Test
Tukey_test <- aov(value ~ Description2 + Rep, data = data_set) %>%
  HSD.test("Description2", group = TRUE) %>%
  .$groups %>%
  as_tibble(rownames = "Description2") %>%
  rename(Letters_Tukey = groups) %>%
  select(-value) %>%
  left_join(maxs, by = "Description2")

# Define Colors and Labels
B=c("0hr"='0h',
    "NMC_48hr"='No SYNCOM(48h)',
    "NMC_96hr"='No SYNCOM(96h)',
    "MC_48hr"='SYNCOM(48h)',
    "MC_96hr"='SYNCOM(96h)')

mypal = c('#258a98','#c9de57','#7e9103','#d388ae','#8f4694')

# -----------------------------------------------------------------------------------
# plot
p = ggplot(data_set) + aes(x=Description2, y=value) +
  geom_boxplot(outlier.shape=NA, lwd=0.2,aes(fill=Description2), colour = "black", alpha=0.7) +
  geom_text(data=Tukey_test, aes(label=Letters_Tukey), size=4, colour="black") +
  scale_fill_manual(values = mypal)
dat <- ggplot_build(p)$data[[1]]
A1 = p + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), colour="black", size=0.5) +
  theme_classic() + ylab("Shannon Diversity Index") + xlab("Time (hr)") + theme(legend.position="none") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=8, color = "black"),
        axis.text.x   = element_text(angle=35, hjust = 1, vjust = 1, size=7, color = "black"),
        axis.title.y  = element_text(size=8, face="bold"),
        axis.title.x  = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", size = 8),
        axis.ticks = element_line(colour = "black", size = 0.7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.35),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_blank()) +
  geom_point(aes(fill=Description2), size=2, shape=21, alpha=0.8, colour="black", stroke=0.1, position=position_jitterdodge(0.5)) +
  scale_x_discrete(label = B)
A1
ggsave("Fig5a.pdf", height=2.4, width=1.6, units='in')






