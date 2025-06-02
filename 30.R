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
#remove samples with less than 1000 total reads
carbom.sub.prun = prune_samples(sample_sums(carbom)>=500, carbom)

#remove a specific taxa
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Kingdom !="Viridiplantae")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Kingdom !="NA")
carbom.sub.prun <- subset_taxa(carbom.sub.prun, Kingdom !="Metazoa")
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
# PERMANOVA
ps1 = carbom.rarefied

#Calculate Beta diversity
set.seed(1)
bray.pca.dist <- phyloseq::distance(ps1, method = "bray")# Calculate bray curtis distance matrix
metadf <- data.frame(sample_data(ps1))# make a data frame from the sample_data

#PERMANOVA analysis: test main effects with interactions
adonis2(bray.pca.dist ~ Treatment * Time + Rep, data = metadf, permutations = 9999)# Adonis test

# -------------------------------------------------------------------------------------------------
#Plot Variance
LC = read.csv("PERMANOVA_SYNCOM_ITS.csv")
data_melt <- melt(LC)

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)
data_melt$Effect <- factor(data_melt$Effect, levels=unique(data_melt$Effect))

mypal = c("Residual"='white',
          "Treatment"='#c4b496',
          "Time"='#5b7a97',
          "Treatment:Time"='#003f5e'
)

p=ggplot(data=data_melt, aes(x=variable, y=value, fill=Effect)) +
  geom_bar(stat="identity", width=0.5, colour="black",size=0.4)
A1 = p + scale_fill_manual(values = mypal) +
  theme_classic() + guides(color = FALSE, fill=FALSE) +
  theme(axis.text.y   = element_text(angle=0, hjust = 1, vjust = 0.5, size=12, color = "black"),
        axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=12, face="bold"),
        axis.title.x  = element_blank(),
        legend.text=element_text(size=12),
        axis.ticks = element_line(colour = "black", size = 0.4),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "white", size = 1),
        panel.border =element_blank(), plot.title = element_text(hjust = 0.5, face="bold", size = 0),
        axis.line.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill=guide_legend(title="")) +
  labs(x="", y="Variance explained (%)") +
  theme(legend.position="right") +
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0))) +
  guides(fill = guide_legend(override.aes=list(size = 6),title="Location", ncol=1))
A1

# -------------------------------------------------------------------------------------------------
##CAP plot
ps1 = carbom.rarefied

#Effect of Treatment
R_value <- grobTree(textGrob(label = ~italic(R)^{"2"}~{"="}~0.51885, x=0.01, y=0.4, hjust=0, gp=gpar(col="black", fontsize=10)))
P_value <- grobTree(textGrob(label = ~italic(p)~{"<"}~"0.0001", x=0.01, y=0.33, hjust=0, gp=gpar(col="black", fontsize=10)))
#Rename samples
B=c('NMC'="No SYNCOM", "MC"="SYNCOM")
ps1
ps1@sam_data$Treatment<- factor(ps1@sam_data$Treatment, levels = c("NMC","MC"))
ordcap = ordinate(ps1, "CAP", "bray", ~ Treatment + Condition(Time + Rep))
distBC_CAP.plot = plot_ordination(ps1, ordcap, "Sample_name", color="Treatment")
distBC_CAP.plot$layers
distBC_CAP.plot$layers <- distBC_CAP.plot$layers[-1]
mypal = c('#7e9103', '#c9de57')
A2 = distBC_CAP.plot + 
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed") +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed") +
  geom_point(aes(fill=Treatment, color=Treatment), shape = 21, size = 6, colour = "black",stroke=0.6) +
  scale_fill_manual(values = mypal, labels = B) +
  theme_classic() + guides(color="none", shape="none") + guides(fill = guide_legend(title="Time")) + 
  theme(axis.text.y   = element_text(angle=0, hjust = 0, vjust = 0.5, size=12, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=12, color = "black"),
        axis.title.y  = element_text(size=12, face="bold"),
        axis.title.x  = element_text(size=12, face="bold"),
        legend.title = element_blank(),
        legend.text=element_text(size=10),
        legend.key.size = unit(0.07, "cm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.7), plot.title = element_text(hjust = 0.5, face="bold", size=10)) +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=1, override.aes = list(size = 7))) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) +
  annotation_custom(R_value) + annotation_custom(P_value)
A2

##..................................................................................................................
#Effect of Time
R_value <- grobTree(textGrob(label = ~italic(R)^{"2"}~{"="}~0.22409, x=0.1, y=0.18, hjust=0, gp=gpar(col="black", fontsize=10)))
P_value <- grobTree(textGrob(label = ~italic(p)~{"<"}~"0.0001", x=0.1, y=0.1, hjust=0, gp=gpar(col="black", fontsize=10)))
ps1
ps1@sam_data$Time<- factor(ps1@sam_data$Time, levels = c("0hr","48hr","96hr"))
ordcap = ordinate(ps1, "CAP", "bray", ~ Time + Condition(Treatment + Rep))
distBC_CAP.plot = plot_ordination(ps1, ordcap, "Sample_name", color="Time")
distBC_CAP.plot$layers
distBC_CAP.plot$layers <- distBC_CAP.plot$layers[-1]
mypal = c('#e3c9a8','#baa692','#4a8182')
mypal = c('#b33a89','#e56281','#ede1e1')
A3 = distBC_CAP.plot + 
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed") +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed") +
  geom_point(aes(fill=Time, color=Time), shape = 21, size = 6, colour = "black",stroke=0.6) +
  scale_fill_manual(values = mypal) +
  theme_classic() + guides(color="none", shape="none") + guides(fill = guide_legend(title="Time")) + 
  theme(axis.text.y   = element_text(angle=0, hjust = 0, vjust = 0.5, size=12, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=12, color = "black"),
        axis.title.y  = element_text(size=12, face="bold"),
        axis.title.x  = element_text(size=12, face="bold"),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.07, "cm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.7), plot.title = element_text(hjust = 0.5, face="bold", size=10)) +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(ncol=1, override.aes = list(size = 7))) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) +
  annotation_custom(R_value) + annotation_custom(P_value) #+
A3

# -------------------------------------------------------------------------------------------------
##empty plot
E=ggplot() + theme_void()

### final plot merged
A1a = A1 + guides(color = "none", shape="none", fill="none")
A2a = A2 + guides(color = "none", shape="none", fill="none")
A3a = A3 + guides(color = "none", shape="none", fill="none")

Plot5 = egg::ggarrange(A1a,E,A2a,A3a, ncol = 4, widths = c(0.07,0.03,1.2,1.2), heights=c(1))
ggsave("Fig5d.pdf", height=2.7, width=6, units='in', plot = Plot5)

### extract legend option 1
library(ggpubr)
library(egg)
# Extract the legend. Returns a gtable
A2_L <- get_legend(A2)
A1_L <- get_legend(A1)
A3_L <- get_legend(A3)
# Convert extracted legend to a ggplot
A2_L = as_ggplot(A2_L)
A1_L = as_ggplot(A1_L)
A3_L = as_ggplot(A3_L)

##plot
legend = egg::ggarrange(A2_L, A3_L, A1_L,ncol = 3, widths = c(1,1,1), heights=c(1))
ggsave("Fig5d_legend.pdf", height=3.5, width=8.5, units='in', plot = legend)


