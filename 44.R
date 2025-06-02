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
library(RColorBrewer)
library(microbiome)
library(microbiomeutilities)
library(ggpubr)
library(data.table)
library(ggtree)
library(gridExtra)
library(grid)
library(scales)
library(BiocParallel)
library(tidyverse)

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
# Data filtering
carbom2.prun = prune_samples(sample_sums(carbom)>=500, carbom)
carbom2.prun = subset_samples(carbom2.prun, Sample_description3 != "Blank")#blank
carbom2.prun = subset_taxa(carbom2.prun, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
carbom2.prun <- prune_taxa(taxa_sums(carbom2.prun) > 2, carbom2.prun)

carbom2.prun

# -------------------------------------------------------------------------------------------------
# Subset fermentation samples (removing unwanted groups)
unwanted_samples <- c("soil", "leaf", "pod", "hand", "tool", "fly", "fermentation_box")
carbom2.prun <- subset_samples(carbom2.prun, !Sample_description3 %in% unwanted_samples)

# -------------------------------------------------------------------------------------------------
#subset Santander
carbom2.prun.sub=subset_samples(carbom2.prun, Location %in% c("Santander"))

#select a specific taxa
carbom2.prun.sub <- subset_taxa(carbom2.prun.sub, Domain =="Fungi")

# -------------------------------------------------------------------------------------------------
#Merge at the genus level
carbom2.prun.sub <- tax_glom(carbom2.prun.sub, "Genus")
carbom2.prun.sub

#Filter method 1: remove taxa with uncharacterized (NA) Genus
carbom2.prun.sub = subset_taxa(carbom2.prun.sub, !is.na(Genus) & !Genus %in% c("", "uncharacterized"))
carbom2.prun.sub

# -------------------------------------------------------------------------------------------------
#build the model
dds <- phyloseq_to_deseq2(carbom2.prun.sub, ~ Sample_description3)
dds
colData(dds)

# -------------------------------------------------------------------------------------------------
#Option #3
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds), 1, gm_mean)
diagdds = estimateSizeFactors(dds, geoMeans = geoMeans)
dds = DESeq(diagdds, fitType="local")

# -------------------------------------------------------------------------------------------------
#Build contrasts
res1 = results(dds, contrast=c("Sample_description3", "fermentation_24hrs", "fermentation_0hrs")) %>% as.data.frame
res1$Id <- rownames(res1)
res1$Contrast <- rep("res1",nrow(res1))
res1$Facet <- rep("fermentation",nrow(res1))
res1 = cbind(as(res1, "data.frame"), as(tax_table(carbom2.prun.sub)[rownames(res1), ], "matrix"))

res2 = results(dds, contrast=c("Sample_description3", "fermentation_48hrs", "fermentation_0hrs")) %>% as.data.frame
res2$Id <- rownames(res2)
res2$Contrast <- rep("res2",nrow(res2))
res2$Facet <- rep("fermentation",nrow(res2))
res2 = cbind(as(res2, "data.frame"), as(tax_table(carbom2.prun.sub)[rownames(res2), ], "matrix"))

res3 = results(dds, contrast=c("Sample_description3", "fermentation_72hrs", "fermentation_0hrs")) %>% as.data.frame
res3$Id <- rownames(res3)
res3$Contrast <- rep("res3",nrow(res3))
res3$Facet <- rep("fermentation",nrow(res3))
res3 = cbind(as(res3, "data.frame"), as(tax_table(carbom2.prun.sub)[rownames(res3), ], "matrix"))

res4 = results(dds, contrast=c("Sample_description3", "fermentation_96hrs", "fermentation_0hrs")) %>% as.data.frame
res4$Id <- rownames(res4)
res4$Contrast <- rep("res4",nrow(res4))
res4$Facet <- rep("fermentation",nrow(res4))
res4 = cbind(as(res4, "data.frame"), as(tax_table(carbom2.prun.sub)[rownames(res4), ], "matrix"))

# -------------------------------------------------------------------------------------------------
#Number of taxa with significant differential abundance
sum(res1$padj < 0.05, na.rm=TRUE)
sum(res2$padj < 0.05, na.rm=TRUE)
sum(res3$padj < 0.05, na.rm=TRUE)
sum(res4$padj < 0.05, na.rm=TRUE)

# -------------------------------------------------------------------------------------------------
#merge contrast results
Res_fer <- rbind(res1,res2,res3,res4)

# Extract significant genes with padj < 0.05 and absolute log2FoldChange > 2
sig_taxa <- Res_fer %>%
  filter(padj < 0.05, abs(log2FoldChange) > 2) %>%
  select(Id) %>%
  distinct()

# Select only significant genes from Res_fer
Res_fer <- Res_fer %>%
  filter(Id %in% sig_taxa$Id)

# -------------------------------------------------------------------------------------------------
# Import taxID info from excel file
tax_ID <- read_excel("colombia_metagenome.xlsx", sheet = "taxID")

#work with only taxa of interest
tax_ID=tax_ID[tax_ID$ASV %in% sig_taxa$Id,]

#save file
write.csv(tax_ID, file = "tax_ID.csv", row.names = FALSE)

# -------------------------------------------------------------------------------------------------
#import tree
tree <- read.tree("santander_tree_fungi_genus_v2.nw")
tree

#view clustering with ggtree
p=ggtree(tree, color="black", size=0.2, branch.length='none')
p + coord_cartesian(clip = 'off') +
  theme(panel.background = element_blank(),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)

# Get the order of tip labels from the tree (TaxID order)
yaxis_order <- tree %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

# ------------------------------------------plot phylogenetic tree for taxa-------------------------------------------------------
#Import annotation file
metadata = read.csv("tax_ID.csv", header=TRUE)

# Remove multiple columns in one step
metadata <- metadata %>% select(-c(Id, ASV, name))

#add colunm for naming
metadata$taxID2 = metadata$taxID

#determine if all samples are present in the tree file and vice versa
metadata$taxID %in% tree$tip.label
tree$tip.label %in% metadata$taxID

#taxa present in the tree but not metadata file
setdiff(tree$tip.label, metadata$taxID)

#taxa present in metadata file but not in the tree
setdiff(metadata$taxID, tree$tip.label)

# -------------------------------------------------------------------------------------------------
#plot
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)

mypal=c("black",'#a0ae00','#e0e004','#195850','#535121','#e0aa02','#f35828','#de5275','#8e67a0',"#c1204c",'#5e3a49','#d4c29e',"#c1204c","pink",'#6aceee','#117845',"#55003b")
p = ggtree(tree, color="black", size=0.2)  %<+% metadata +
  geom_tiplab(aes(label=taxID2), offset=1.5, hjust = 0, vjust=0.5, size=0, align = TRUE, linetype = "dotted", linesize = 0.2) +
  geom_tippoint(aes(fill=Phylum), size =0.7, alpha=1, shape=21, stroke=0) +
  scale_fill_manual(values = mypal) + guides(color=FALSE)
p + theme(legend.title = element_text(size=6, face="bold"),
          legend.text=element_text(size=4),
          legend.key.size = unit(0.1, "cm"),
          plot.margin=margin(0, 0, 0, 0),
          legend.position = 'right') +
  guides(fill = guide_legend(override.aes = list(shape = 22, size =3)))

# --------------------------------------------reorder x and y axes-----------------------------------------------------
#select columns
columns_to_keep = c("ASV", "taxID")
tax_ID = subset(tax_ID, select=columns_to_keep)

#rename column
colnames(tax_ID)[which(names(tax_ID) == "ASV")] <- "Id"

# merge taxIDs to Res_fer by ID
Res_fer <- merge(Res_fer,tax_ID,by="Id")

#order the x axis
Res_fer$Contrast <- factor(Res_fer$Contrast, levels = c('res1','res2','res3','res4'))

#Reorder y axis taxIDs by phylogenetic tree taxa order
Res_fer$taxID <- factor(Res_fer$taxID, levels = yaxis_order)

# Reorder the 'Genus' column based on taxID order and yaxis_order
Res_fer <- Res_fer %>%
  arrange(factor(taxID, levels = yaxis_order)) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# ----------------------make a dataframe with x and y coordinates of significant taxa to highlight---------------------------------------------------------------------------
#sort heat map matrix based on yaxis_order
Res_fer = Res_fer %>% arrange(factor(taxID, levels = yaxis_order))

#make the list of cells to highlight as significant on y axis (note if there is NA in this column, this would cause problems)
Res_fer$y_position <- factor(Res_fer$taxID)
Res_fer$y_position = as.numeric(Res_fer$y_position)

#make the list of cells to highlight as significant on x axis
Res_fer$x_position <- factor(Res_fer$Contrast)
Res_fer$x_position = as.numeric(Res_fer$x_position)

#select columns
frames = select(Res_fer, c('Id','taxID','Genus','log2FoldChange','padj','x_position','y_position'))

#extract significant results (padj < 0.05)
sig_values = subset(frames, padj < 0.05)

# filter significant results on min log2fold change
sig_values = subset(sig_values, abs(log2FoldChange) > 2)

# -------------------------------------------------------------------------------------------------
#plot heatmap
B=c('0h vs 24h','0h vs 48h','0h vs 72h','0h vs 96h')

p = ggplot(data=Res_fer) +
  geom_raster(aes(x=Contrast, y=Genus, fill=log2FoldChange)) +
  scale_fill_gradient2(low="blue", high="red", na.value="grey", name="log2FoldChange") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.1)) +
  scale_y_discrete(position = "right", expand = c(0,0)) + scale_x_discrete(label = B, expand = c(0,0)) +
  geom_rect(data=sig_values, size=0.2, fill=NA, colour="black", aes(xmin=x_position - 0.4, xmax=x_position + 0.4, ymin=y_position - 0.4, ymax=y_position + 0.4))
p
p + guides(color=FALSE, shape=FALSE, fill=FALSE) +
  theme(axis.text.x=element_text(angle=60, hjust = 1, vjust = 1, size=8, color="black"),
        axis.text.y=element_text(angle=0, hjust = 1, vjust = 0.5, size=5, color="black"),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.title = element_blank(),
        legend.text=element_blank()) +
  labs(x = "Fermentation time (h)") +
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"))

# ------------------------ merge plots-------------------------------------------------------------------------
mypal=c('Dipodascaceae'='#34d6c9',
        'Phaffomycetaceae'='#00a6c0',
        'Debaryomycetaceae'='#055895',
        'Trichomonascaceae'='#9fc2aa',
        'Ustilaginaceae'='#a6e65a',
        'Aspergillaceae'='#027850',
        'Nectriaceae'='#b3a281',
        'Glomerellaceae'='#7d674d',
        'Sordariaceae'='#ff90c5',
        'Chaetomiaceae'='#fe499a',
        'Sclerotiniaceae'='#a51672',
        'Mycosphaerellaceae'='#d83a35',
        'Saccharomycetaceae'='#fe7e01',
        'Pichiaceae'='#fcec65',
        'Cryptococcaceae'='#ad9bc9',
        'Malasseziaceae'='#7b49ac',
        'Pyriculariaceae'='#ad27ac',
        'Clavicipitaceae'='#617b7d')
##build phylo tree
A1 = ggtree(tree, color="black", size=0.1)  %<+% metadata +
  geom_tiplab(aes(label=taxID2), offset=2, hjust = 0, vjust=0.5, size=0, align = TRUE, linetype = "dotted", linesize = 0.05) +
  geom_tippoint(aes(color=Family), size =1.5, alpha=1, shape=16, stroke=0) +
  scale_color_manual(values = mypal) + theme(plot.margin=margin(0, 0, 0, 0))
A1
A1 = A1 + theme(legend.title = element_text(size=5, face="bold"),
                legend.text=element_text(size=5),
                legend.key.size = unit(0.04, "cm"),
                plot.margin=margin(0, 0, 0, 0),
                legend.position = 'left') +
  guides(color = guide_legend(override.aes = list(shape = 16, size =3), ncol=2))
A1

##build heatmap
B=c('0h vs 24h','0h vs 48h','0h vs 72h','0h vs 96h')
A2 = ggplot(data=Res_fer) +
  geom_raster(aes(x=Contrast, y=Genus, fill=log2FoldChange)) +
  scale_fill_gradient2(low="blue", high="red", na.value="grey", name="log2FoldChange") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.4)) +
  scale_y_discrete(position = "right", expand = c(0,0)) + scale_x_discrete(label = B, expand = c(0,0)) +
  geom_rect(data=sig_values, size=0.2, fill=NA, colour="black", aes(xmin=x_position - 0.4, xmax=x_position + 0.4, ymin=y_position - 0.4, ymax=y_position + 0.4)) +
  guides(color=FALSE, shape=FALSE) +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size=4, color="black"),
        axis.text.y=element_text(angle=0, hjust = 1, vjust = 0.5, size=4, color="black"),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(x = "Fermentation time (h)") +
  guides(fill = guide_colourbar(barwidth = 4,
                                barheight = 0.75,
                                direction = "horizontal",
                                ticks = TRUE,
                                title.position = 'top',
                                ticks.colour = "black")) +
  theme(
    legend.title = element_text(size = 5, angle = 0, face="bold"),
    legend.text = element_text(size = 4.5),
    legend.title.align = 0.5,
    legend.direction = "vertical",
    legend.margin=margin(c(0,0,0,0)))
A2

# ------------------------------------------------------------------------------------------------
##empty plot
A4=ggplot() + theme_void()

A1a = A1 + guides(color=FALSE, fill=FALSE)
A2a = A2 + guides(color=FALSE, fill=FALSE)

##plot figure
Plot2 = egg::ggarrange(A4,A1a,A2a, ncol = 3, widths = c(0.4,0.3,2))
Plot2
ggsave("Supplementary_Fig1e.pdf", height=5, width=4.55, units='cm', plot = Plot2)

# --------------------------------------------------------------------------------------------
### extract legend option 1
library(ggpubr)
library(egg)

# Extract the legend. Returns a gtable
A1_leg <- get_legend(A1)
A2_leg <- get_legend(A2)

# Convert extracted legend to a ggplot object
LEG1 = as_ggplot(A1_leg)
LEG2 = as_ggplot(A2_leg)

##plot legends
A3 = egg::ggarrange(LEG1, LEG2, ncol=1, heights=c(1,1))
A3
ggsave("Supplementary_Fig1e_legend.pdf", height=5, width=6.5, units='in', plot = A3)

