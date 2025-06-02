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
#-------------------------------------------------------------------------------------------------
# Load necessary libraries
library(seqinr)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(ggplot2)
library(BiocParallel)
library(readxl)
library(ggtreeExtra)
library(dplyr)
library(reshape2)
library(treeio)
library(ggnewscale)
library(msa)
library(phangorn)
library(microseq)

#-------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

#-------------------------------------------------------------------------------------------------
# load the sequences from the file
# change "DNA" to "RNA" or "AA" if necessary
seqs <- readAAStringSet("MAGs_55_concatenated_gene_sequences.fasta", format = "fasta")

##look at some of the sequences (optional)
seqs

##get fasta header
headers <- names(seqs)
headers

#-----------------------multiple sequence alignment--------------------------------------------------------------------------
#aligned<-msa(seqs,"ClustalW")

# Option 1: write the MAS alignment to a new FASTA file
#writeXStringSet(as(unmasked(aligned), "XStringSet"), file="55MAG_aligned_clustalW.fasta")

#------------------------trim-------------------------------------------------------------------------
##read alignment FASTA file for microseq
df <- readFasta("55MAG_aligned_clustalW.fasta")

## trim aligment with microseq
dft = msaTrim(df, gap.end = 0.5, gap.mid = 0.75)#gap.mid = 0.9

# save trim alignment to FASTA file
writeFasta(dft, "55MAG_aligned_clustalW_trim.fasta", width=80)

#---------------------distance matrix----------------------------------------------------------------------------
# read in the aligned data
dna <- read.alignment("55MAG_aligned_clustalW_trim.fasta", format = "fasta")

# create a distance matrix for the alignment  with seqinr
D <- dist.alignment(dna, matrix = "similarity")

#-----------------------plot--------------------------------------------------------------------------
#create tree object with ape
tree <- njs(D)
class(tree) #all trees created using {ape} package will be of class phylo
tree <- ladderize(tree)

#check the formatting of the tip.labels in the tree
tree$tip.label

#plot tree
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5, options(ignore.negative.edge=TRUE))
p

#-----------------------get metadata--------------------------------------------------------------------------
# Import metadata
metadata <- read.csv("MAG_55_metadata.csv")

# inspect the first 6 Sample_IDs
head(metadata$Sample_ID) 

#determine if all samples are present in the tree file and vice versa
metadata$Sample_ID %in% tree$tip.label
tree$tip.label %in% metadata$Sample_ID

#show any sample IDs that are not on the tree
metadata$Sample_ID[!tree$tip.label %in% metadata$Sample_ID]#character(0) = none

#--------------------------plot-----------------------------------------------------------------------
##plot
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5, options(ignore.negative.edge=TRUE))
p
p1 = rotate_tree(p, -90)# rotate tree
p1

#-------------------------get nodes------------------------------------------------------------------------
#Get internal node numbers
ggtree(tree) + geom_text(aes(label=node), hjust=-.5,  size=2)
ggsave("tree_circular.jpeg", width = 8, height = 8,  units='in', dpi=600)

ggtree(tree, layout="equal_angle", open.angle=10, size=0.1) + geom_text(aes(label=node), hjust=-.3,  size=2)
ggsave("tree_node_number.pdf", width = 8, height = 8,  units='in')

#------------------------annotate-------------------------------------------------------------------------
mypal=c("Ascomycota" = "#f14a98",
        "Bacteroidetes" = "#a6dbdf",
        "Firmicutes" = "#9ACD32",
        "Proteobacteria" = "#BC80BD"
)

A1 = ggtree(tree, color="black", layout="equal_angle", size=0.7, options(ignore.negative.edge=TRUE)) %<+% metadata +
  #geom_tippoint(aes(color=Phylum), size = 2, alpha=1, shape=16) +
  #geom_tiplab(aes(label=Name), offset=0.01, size=1, color="black") +##,    align = FALSE, linetype = "blank", linesize = 0.5, hjust = 0, vjust=0, 
  geom_hilight(node=61, fill="#8e67a0", alpha = 0.6, extendto=0.1, color="#8e67a0", size=0.05) +
  geom_hilight(node=64, fill="#8e67a0", alpha = 0.6, extendto=0.1, color="#8e67a0", size=0.05) +
  geom_hilight(node=65, fill="#9ACD32", alpha = 0.6, extendto=0.1, color="#9ACD32", size=0.05) +
  geom_hilight(node=67, fill="#a6dbdf", alpha = 0.6, extendto=0.1, color="#a6dbdf", size=0.05) +
  geom_hilight(node=74, fill="#f14a98", alpha = 0.6, extendto=0.1, color="#f14a98", size=0.05) +
  geom_hilight(node=79, fill="#f14a98", alpha = 0.6, extendto=0.1, color="#f14a98", size=0.05) +
  geom_hilight(node=92, fill="#f14a98", alpha = 0.6, extendto=0.1, color="#f14a98", size=0.05) +
  scale_color_manual(values = mypal) + #guides(color=FALSE) +
  guides(color=guide_legend(nrow=5, override.aes = list(size = 5))) +
  theme(legend.title = element_text(size=6, face="bold"),
        legend.text=element_text(size=6),
        legend.position="right",
        legend.key.size = unit(0.4, "cm")) +
  labs(fill = "Phylum")
A1
ggsave("Fig4a.pdf", height=7, width=9, units='in')


# Create a data frame for legend only
legend_df <- data.frame(Phylum = names(mypal))
# Plot legend
legend_plot <- ggplot(legend_df, aes(x = Phylum, y = 1, fill = Phylum)) +
  geom_point(aes(fill=Phylum), shape = 21, size = 4, colour = "black",stroke=0.3) +
  scale_fill_manual(values = mypal) +
  theme_void() +
  theme(legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.key.size = unit(0, "cm"),
        legend.spacing = unit(0, "cm"),  # Space between legend items
        legend.margin = margin(0, 0, 0, 0), # Margins around legend box
        ) +
  guides(fill = guide_legend(title = "Phylum", nrow = 4)) +
  guides(color = "none", shape="none") #+
#scale_x_discrete(expand = c(0.1, 0.1)) + scale_y_discrete(expand = c(0.1, 0.1)) +
# Show the legend
legend_plot
ggsave("Fig4a_legend.pdf", height=7, width=7, units='in', plot = legend_plot)

