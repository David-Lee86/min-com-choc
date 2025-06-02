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
# Set seed for reproducibility
set.seed(1)

# Define dataset
ps1 <- carbom.rarefied.sub
metadf <- data.frame(sample_data(ps1))  # Convert sample data to a data frame

# Calculate Bray-Curtis distance matrix
bray_dist <- phyloseq::distance(ps1, method = "bray")

# PERMANOVA analysis
diversity_result <- adonis2(bray_dist ~ Time2 + Location + Harvest + Replicate, 
                            data = metadf, permutations = 9999)
print(diversity_result)

# -------------------------------------------------------------------------------------------------
# Effect of Time
mypal <- c('#2d4800','#528501','#a9d11c','#eaabcd','#c0047b')
B = c("0hr"="0h",
      "24hr"="24h",
      "48hr"="48h",
      "72hr"="72h",
      "96hr"="96h")
R_value <- grobTree(textGrob(label = ~italic(R)^{"2"}~{"="}~0.49173, x=0.40, y=0.13, hjust=0, gp=gpar(col="black", fontsize=18)))
P_value <- grobTree(textGrob(label = ~italic(p)~{"<"}~"0.0001", x=0.40, y=0.05, hjust=0, gp=gpar(col="black", fontsize=18)))

ps1@sam_data$Time2 <- factor(ps1@sam_data$Time2, levels = c("0hr", "24hr", "48hr", "72hr", "96hr"))
ordcap <- ordinate(ps1, "CAP", "bray", ~ Time2 + Condition(Location + Harvest + Replicate))
distBC_CAP.plot <- plot_ordination(ps1, ordcap, "Sample_name", color = "Time2")
distBC_CAP.plot$layers <- distBC_CAP.plot$layers[-1]
A1 <- distBC_CAP.plot +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_point(aes(fill = Time2), shape = 21, size = 7, color = "black", stroke = 0.6) +
  scale_fill_manual(values = mypal, label = B) +
  theme_classic() +
  theme(
    axis.text = element_text(size=18, color = "black"),
    axis.title = element_text(size=18, face = "bold"),
    legend.text = element_text(size=18),
    legend.key.size = unit(0.07, "cm"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  annotation_custom(R_value) + annotation_custom(P_value) +
  guides(fill = guide_legend(title = "Time", nrow = 2, override.aes = list(size = 7)))
A1

# Effect of Location
mypal <- c("#fc0558", "#9eb744", '#526400')
R_value <- grobTree(textGrob(label = ~italic(R)^{"2"}~{"="}~0.06982, x=0.42, y=0.85, hjust=0, gp=gpar(col="black", fontsize=18)))
P_value <- grobTree(textGrob(label = ~italic(p)~{"="}~"0.0053", x=0.42, y=0.77, hjust=0, gp=gpar(col="black", fontsize=18)))

ps1@sam_data$Location <- factor(ps1@sam_data$Location, levels = c('Huila','Santander','Antioquia'))
ordcap <- ordinate(ps1, "CAP", "bray", ~ Location + Condition(Time2 + Harvest + Replicate))
distBC_CAP.plot <- plot_ordination(ps1, ordcap, "Sample_name", color = "Location")
distBC_CAP.plot$layers <- distBC_CAP.plot$layers[-1]
A2 <- distBC_CAP.plot +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_point(aes(fill = Location), shape = 21, size = 7, color = "black", stroke = 0.6) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  theme(
    axis.text = element_text(size=18, color = "black"),
    axis.title = element_text(size=18, face = "bold"),
    legend.text = element_text(size=18),
    legend.key.size = unit(0.07, "cm"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  annotation_custom(R_value) + annotation_custom(P_value) +
  guides(fill = guide_legend(title = "Location", nrow = 2, override.aes = list(size = 7)))
A2


# Effect of Harvest
mypal <- c("plum4", "darkmagenta")
R_value <- grobTree(textGrob(label = ~italic(R)^{"2"}~{"="}~0.03075, x=0.40, y=0.95, hjust=0, gp=gpar(col="black", fontsize=18)))
P_value <- grobTree(textGrob(label = ~italic(p)~{"="}~"0.0309", x=0.40, y=0.87, hjust=0, gp=gpar(col="black", fontsize=18)))

ps1@sam_data$Harvest <- factor(ps1@sam_data$Harvest, levels = c("H1", "H2"))
ordcap <- ordinate(ps1, "CAP", "bray", ~ Harvest + Condition(Time2 + Location + Replicate))
distBC_CAP.plot <- plot_ordination(ps1, ordcap, "Sample_name", color = "Harvest")
distBC_CAP.plot$layers <- distBC_CAP.plot$layers[-1]
A3 <- distBC_CAP.plot +
  geom_hline(yintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray70", linetype = "dashed") +
  geom_point(aes(fill = Harvest), shape = 21, size = 7, color = "black", stroke = 0.6) +
  scale_fill_manual(values = mypal, labels = c("Harvest 1", "Harvest 2")) +
  theme_classic() +
  theme(
    axis.text = element_text(size=18, color = "black"),
    axis.title = element_text(size=18, face = "bold"),
    legend.text = element_text(size=18),
    legend.key.size = unit(0.07, "cm"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.6),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  annotation_custom(R_value) + annotation_custom(P_value) +
  guides(fill = guide_legend(title = "Harvest", nrow = 1, override.aes = list(size = 7)))
A3


# -------------------------------------------------------------------------------------------------------
# Empty plot
E <- ggplot() + theme_void()

### final plot merged
A1a = A1 + guides(color = "none", shape="none", fill="none")
A2a = A2 + guides(color = "none", shape="none", fill="none")
A3a = A3 + guides(color = "none", shape="none", fill="none")

# Plot figure
Plot2 <- egg::ggarrange(A1a, A2a, A3a, ncol = 3, widths = c(1, 1, 1))
Plot2

# Save the plot
ggsave("Fig2c.pdf", height = 11, width = 27, units = 'cm', plot = Plot2)


# ----------------------------------------------------------------------------------
### extract legend option 1
library(ggpubr)
library(egg)

#Rename labels
# Extract the legend. Returns a gtable
A1_L <- get_legend(A1a)
A2_L <- get_legend(A2a)
A3_L <- get_legend(A3a)
# Convert extracted legend to a ggplot
A1_L = as_ggplot(A1_L)
A2_L = as_ggplot(A2_L)
A3_L = as_ggplot(A3_L)

#plot
LEG <- egg::ggarrange(A1_L, A2_L, A3_L, ncol = 3, widths = c(1, 1, 1))
LEG
ggsave("Fig2c_legend.pdf", height=5, width=3, units='in', plot = LEG)


