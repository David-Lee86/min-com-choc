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
library(adegenet)
library(ape)
library(BiocParallel)
library(DECIPHER)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(ggnewscale)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(grid)
library(gridExtra)
library(hrbrthemes)
library(knitr)
library(MASS)
library(RColorBrewer)
library(readr)
library(readxl)
library(reshape2)
library(rmdformats)
library(scales)
library(seqinr)
library(treeio)
library(viridis)
library (ComplexHeatmap)
library(circlize)
library(magick)
library(tidyverse)
#-------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

######################################################### 44 MAGs ################################
# Import data metabolome data matrix
AC <- read.table("MAGs_44_metabolites.tsv", sep = "\t", header = TRUE)

AC <- AC %>% rename("MAG" = "X")

# sort data
AC <- AC[order(AC$MAG),]

# Replace 'bin_001A' with 'bin_001' in the AC$MAG column
AC$MAG <- gsub("bin_001A", "bin_001", AC$MAG)

#replace name of values in column
AC$MAG <- str_replace(AC$MAG, "bin_", "m")

# Define the row names from the sample_name column
AC = AC %>% remove_rownames %>% column_to_rownames(var="MAG")

#-------------------------------------------------------------------------------------------------
# create distance matrix
D <- dist(AC, method = "euclidean")

#Classical Multidimensional_scaling: project the distance matrix to a lower dimension k,
#where desired k = 2 (x y coordinates), while trying to preserve the distances between data points
Tmds <- cmdscale(D, k = 2)  %>% as.data.frame

#create a column with the row names
Tmds <- tibble::rownames_to_column(Tmds, "Name")

#--------------------------------Add taxonomy information-----------------------------------------------------------------
# Import metadata
metadata <- read.csv("MAG_55_metadata.csv")

# merge two data frames by ID
AC3 <- merge(Tmds, metadata, by="Name")

##Define the row names from the Id column
rownames(AC3) = AC3$Name

#rename axis
colnames(AC3)[which(names(AC3) == "V1")] <- "xvalue"
colnames(AC3)[which(names(AC3) == "V2")] <- "yvalue"

#-------------------------------------------------------------------------------------------------
# Extract unique values from the 'Family' column of AC3
A <- AC3 %>% pull(Family) %>% unique() %>% sort()
# Print each unique value on a separate line
cat(A, sep = "\n")
output <- paste0('"', paste(A, collapse = '", "'), '"')
cat(output)

list2 =c('#96ae02','#e5b4d4','#a797ca','#a7d8ec','#aadbd7','#659e98','#e6efd0',"#f27a63", "#be057c", "#528502")

result <- paste("\"", A, "\"=\"", list2, "\",", sep = "")
cat(result, sep = "\n")

#-------------------------------------------------------------------------------------------------
#plot
mypal=c(
"Acetobacteraceae"="#96ae02",
"Erwiniaceae"="#e5b4d4",
"Lactobacillaceae"="#a797ca",
"Leuconostocaceae"="#a7d8ec",
"Pectobacteriaceae"="#a59579",
"Pichiaceae"="#659e98",
"Rhodanobacteraceae"="#e0e004",
"Saccharomycetaceae"="#f27a63",
"Saccharomycodaceae"="#be057c",
"unclassified"="#e6efd0",
"Bacillaceae"="#176c91"
)

p <-ggplot(AC3, aes(xvalue, yvalue, fill=Name))
A1 = p + 
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed", size = 0.5) +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed", size = 0.5) +
  #geom_point(aes(fill=Family, color=Family), shape = 1, size = 1.5, stroke=1) +
  geom_point(aes(fill=Family, color=Family), shape = 21, size = 3.25, colour = "black",stroke=0.3) +
  scale_fill_manual(values = mypal) +
  scale_color_manual(values = mypal) +
  theme_classic() +
  theme(axis.text.y   = element_text(angle=0, hjust = 1, vjust = 0.5, size=8, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=8, color = "black"),
        axis.title.y  = element_text(size=9, face="bold"),
        axis.title.x  = element_text(size=9, face="bold"),
        legend.title = element_text(size=9, face="bold"),
        legend.text=element_text(size=6),
        legend.key.size = unit(0.1, "cm"),
        plot.title = element_text(size=8, hjust = 0.5, vjust = 0, color = "black", face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3), axis.ticks = element_line(size=0.3),
        panel.border = element_blank()
  ) +
  ylab("Coord.2") + xlab("Coord.1") + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) + 
  labs(title = "Full community") +
  #guides(color="none", shape="none", fill="none") +# + guides(fill = guide_legend(title="Family")) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1)) +
  guides(fill = guide_legend(ncol = 4)) 
A1


###################################################### 10 MAGs #########################################
# Import data metabolome data matrix
LC <- read.table("MAGs_10_metabolites.tsv", sep = "\t", header = TRUE)

LC <- LC %>% rename("MAG" = "X")

# sort data
LC <- LC[order(LC$MAG),]

# Replace 'bin_001A' with 'bin_001' in the LC$MAG column
LC$MAG <- gsub("bin_001A", "bin_001", LC$MAG)

#replace name of values in column
LC$MAG <- str_replace(LC$MAG, "bin_", "m")

# Define the row names from the sample_name column
LC = LC %>% remove_rownames %>% column_to_rownames(var="MAG")

#-------------------------------------------------------------------------------------------------
# create distance matrix
D <- dist(LC, method = "euclidean")

#Classical Multidimensional_scaling: project the distance matrix to a lower dimension k,
#where desired k = 2 (x y coordinates), while trying to preserve the distances between data points
mds <- cmdscale(D, k = 2)  %>% as.data.frame

#create a column with the row names
mds <- tibble::rownames_to_column(mds, "Name")

#-----------------------------------Add taxonomy information-------------------------------------
# Import metadata
metadata <- read.csv("MAG_55_metadata.csv")

# merge two data frames by ID
LC3 <- merge(mds, metadata, by="Name")

##Define the row names from the Id column
rownames(LC3) = LC3$Name

#rename axis
colnames(LC3)[which(names(LC3) == "V1")] <- "xvalue"
colnames(LC3)[which(names(LC3) == "V2")] <- "yvalue"

#-------------------------------------------------------------------------------------------------
#plot
p <-ggplot(LC3, aes(xvalue, yvalue, fill=Name))
A2 = p + 
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed", size = 0.5) +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed", size = 0.5) +
  #geom_point(aes(fill=Family, color=Family), shape = 1, size = 4, stroke=1) +
  geom_point(aes(fill=Family, color=Family), shape = 21, size = 3.25, colour = "black",stroke=0.3) +
  scale_fill_manual(values = mypal) +
  scale_color_manual(values = mypal) +
  theme_classic() +
  theme(axis.text.y   = element_text(angle=0, hjust = 1, vjust = 0.5, size=8, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=8, color = "black"),
        axis.title.y  = element_text(size=9, face="bold"),
        axis.title.x  = element_text(size=9, face="bold"),
        legend.title = element_text(size=9, face="bold"),
        legend.text=element_text(size=10),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size=8, hjust = 0.5, vjust = 0, color = "black", face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3), axis.ticks = element_line(size=0.3),
        panel.border = element_blank()
  ) +
  ylab("Coord.2") + xlab("Coord.1") + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) + 
  labs(title = "Defined community") +
  #guides(color="none", shape="none", fill="none") +# + guides(fill = guide_legend(title="Family")) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1))
A2

############################################### SYNCOM #############################################################
# Import data metabolome data matrix
LC = read.csv("SYNCOM_metabolites.csv")

# sort data
LC <- LC[order(LC$Genome),]

#library(tidyverse)
#replace name of values in column
LC$Genomes <- str_replace(LC$Genomes, "C", "D")

# Define the row names from the sample_name column
LC = LC %>% remove_rownames %>% column_to_rownames(var="Genomes")

#-------------------------------------------------------------------------------------------------
# create distance matrix
D <- dist(LC, method = "euclidean")

#Classical Multidimensional_scaling: project the distance matrix to a lower dimension k,
#where desired k = 2 (x y coordinates), while trying to preserve the distances between data points
mds <- cmdscale(D, k = 2)  %>% as.data.frame

#create a column with the row names
mds <- tibble::rownames_to_column(mds, "Sample_ID2")

#-------------------------------------------------------------------------------------------------
metadata = read.csv("SYNCOM_metadata.csv")

# merge two data frames by ID
LC3 <- merge(mds, metadata, by="Sample_ID2")

##Define the row names
rownames(LC3) = LC3$Sample_ID2

#rename axis
colnames(LC3)[which(names(LC3) == "V1")] <- "xvalue"
colnames(LC3)[which(names(LC3) == "V2")] <- "yvalue"

##plot
p <-ggplot(LC3, aes(xvalue, yvalue, fill=Sample_ID))
A3 = p + 
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed", size = 0.5) +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed", size = 0.5) +
  #geom_point(aes(fill=Family, color=Family), shape = 1, size = 4, stroke=1) +
  geom_point(aes(fill=Family, color=Family), shape = 21, size = 3.25, colour = "black",stroke=0.3) +
  scale_fill_manual(values = mypal) +
  scale_color_manual(values = mypal) +
  theme_classic() +
  theme(axis.text.y   = element_text(angle=0, hjust = 1, vjust = 0.5, size=8, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=8, color = "black"),
        axis.title.y  = element_text(size=9, face="bold"),
        axis.title.x  = element_text(size=9, face="bold"),
        legend.title = element_text(size=9, face="bold"),
        legend.text=element_text(size=10),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size=8, hjust = 0.5, vjust = 0, color = "black", face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.3), axis.ticks = element_line(size=0.3),
        panel.border = element_blank()
  ) +
  ylab("Coord.2") + xlab("Coord.1") + theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) + 
  labs(title = "SYNCOM") +
  #guides(color="none", shape="none", fill="none") +# + guides(fill = guide_legend(title="Family")) +
  scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.1, 0.1))
A3


###################################################################################################################
##empty plot
A4=ggplot() + theme_void()

### final plot merged
A1a = A1 + guides(color = "none", shape="none", fill="none")
A2a = A2 + guides(color = "none", shape="none", fill="none")
A3a = A3 + guides(color = "none", shape="none", fill="none")

##plot figure
Plot2 = egg::ggarrange(A1a,A2a,A3a, nrow = 1, heights = c(1), ncol = 3, widths = c(1,1,1))
Plot2
ggsave("Fig4b.pdf", height=4.5, width=12, units='cm', plot = Plot2)


# Create a data frame for legend only
legend_df <- data.frame(Family = names(mypal))
# Plot legend
legend_plot <- ggplot(legend_df, aes(x = Family, y = 1, fill = Family)) +
  geom_point(aes(fill=Family), shape = 21, size = 3.25, colour = "black",stroke=0.3) +
  scale_fill_manual(values = mypal) +
  theme_void() +
  theme(legend.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=10),
        legend.key.size = unit(0.1, "cm"),
        legend.margin = margin(10, 10, 10, 10),
        plot.margin = margin(1, 1, 1, 1)) +
  guides(fill = guide_legend(title = "Family", ncol = 2)) +#nrow = 4
  guides(color = "none", shape="none") #+
  #scale_x_discrete(expand = c(0.1, 0.1)) + scale_y_discrete(expand = c(0.1, 0.1)) +
# Show the legend
legend_plot
ggsave("Fig4b_legend.pdf", height=5, width=10, units='in', plot = legend_plot)
