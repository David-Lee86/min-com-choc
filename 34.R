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
library(DESeq2)
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
# Import data from excel file
otu_mat<- read_excel("SYNCOM_GCMS.xlsx", sheet = "OTU_matrix")
tax_mat<- read_excel("SYNCOM_GCMS.xlsx", sheet = "Taxonomy_table")
samples_df <- read_excel("SYNCOM_GCMS.xlsx", sheet = "Sample_name")
# -------------------------------------------------------------------------------------------------
# Phyloseq objects need to have row.names. Define the row names from the ASV column
otu_mat = otu_mat %>% remove_rownames %>% column_to_rownames(var="ASV")##adds row names and remove sample_name column

# Define the row names from the tax_mat data frame 
tax_mat = tax_mat %>% remove_rownames %>% column_to_rownames(var="ASV")##adds row names and remove sample_name column

# Define the row names from the samples_df data frame 
row.names(samples_df) <- samples_df$Name

# Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform data to phyloseq objects
ASV = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAM = sample_data(samples_df)
rownames(SAM) <-SAM$Name

sample_names(ASV)
rownames(SAM)

#Merge phyloseq objects
carbom = phyloseq(ASV, TAX, SAM)
carbom

#----------------------------------------------------------------------------------------------------
carbom2.prun = carbom

#trouble shooting dataframes
LC <- as.data.frame(carbom2.prun@otu_table)

##add a pseudo-count value of '1' to all values in datafrme
#carbom2.prun@otu_table = carbom2.prun@otu_table +1

# Check if all rows contain at least one 0
all_rows_contain_zero <- all(apply(LC[, -1] == 0, 1, any))
cat("All rows contain at least one 0:", all_rows_contain_zero, "\n")

#----------------------------------------------------------------------------------------------------
#build the model
dds <- phyloseq_to_deseq2(carbom2.prun, ~ Description4)
dds
colData(dds)

#----------------------------------------------------------------------------------------------------
#Option #1 Differential expression analysis
dds <- DESeq(dds, test="Wald", fitType="parametric")

#----------------------------------------------------------------------------------------------------
#Build contrasts
res1 = results(dds, contrast=c("Description4", "SYNCOM", "No_SYNCOM")) %>% as.data.frame
res1$Id <- rownames(res1)
res1$Contrast <- rep("res1",nrow(res1))
res1$Facet <- rep("fermentation",nrow(res1))
res1 = cbind(as(res1, "data.frame"), as(tax_table(carbom2.prun)[rownames(res1), ], "matrix"))

res2 = results(dds, contrast=c("Description4", "Santander", "No_SYNCOM")) %>% as.data.frame
res2$Id <- rownames(res2)
res2$Contrast <- rep("res2",nrow(res2))
res2$Facet <- rep("fermentation",nrow(res2))
res2 = cbind(as(res2, "data.frame"), as(tax_table(carbom2.prun)[rownames(res2), ], "matrix"))

res3 = results(dds, contrast=c("Description4", "Huila", "No_SYNCOM")) %>% as.data.frame
res3$Id <- rownames(res3)
res3$Contrast <- rep("res3",nrow(res3))
res3$Facet <- rep("fermentation",nrow(res3))
res3 = cbind(as(res3, "data.frame"), as(tax_table(carbom2.prun)[rownames(res3), ], "matrix"))

res4 = results(dds, contrast=c("Description4", "Antioquia", "No_SYNCOM")) %>% as.data.frame
res4$Id <- rownames(res4)
res4$Contrast <- rep("res4",nrow(res4))
res4$Facet <- rep("fermentation",nrow(res4))
res4 = cbind(as(res4, "data.frame"), as(tax_table(carbom2.prun)[rownames(res4), ], "matrix"))

#----------------------------------------------------------------------------------------------------
#Number of DEGs genes
k = 2

##Number of DEGs genes (padj < 0.05, & log2FoldChange > 1)
sum(res1$padj < 0.05 & abs(res1$log2FoldChange) > k, na.rm = TRUE)
sum(res2$padj < 0.05 & abs(res2$log2FoldChange) > k, na.rm = TRUE)
sum(res3$padj < 0.05 & abs(res3$log2FoldChange) > k, na.rm = TRUE)
sum(res4$padj < 0.05 & abs(res4$log2FoldChange) > k, na.rm = TRUE)

#---------------------------merge and select only significant taxa----------------------------------
#merge contrast results
Res_fer <- rbind(res1,res2,res3,res4)

#log2FoldChange > k
k = 2

# Extract list of significant genes with padj < 0.05 and |log2FoldChange| > 2
sig_genes <- Res_fer %>%
  filter(padj < 0.05, abs(log2FoldChange) > k) %>%
  distinct(Id)  # Get unique gene IDs

# Filter the original Res_fer for only the significant genes
Res_fer <- Res_fer %>%
  filter(Id %in% sig_genes$Id)

#-------------------------prepare heatmap data matrix---------------------------------------------
columns_to_keep = c("Id", "Contrast", "log2FoldChange")
heatmap.dataset2 = subset(Res_fer, select=columns_to_keep)

# Convert data to short format (unmelt)
heatmap.dataset2 = reshape2::dcast( heatmap.dataset2 , Id~Contrast, value.var="log2FoldChange") %>% as.data.frame

# Set up row names
row.names(heatmap.dataset2)=heatmap.dataset2$Id

#Format the dataframe for clustering
tree.matrix2 <- heatmap.dataset2 %>% select (-Id)

#-----------------------------Run clustering for y axis-----------------------------------------
#create distance matrix
dd <- dist(scale(tree.matrix2), method = "euclidean")

#Hierarchical Clustering
tree.y <- hclust(dd, method = "ward.D2")

#view clustering with ggtree
p=ggtree(tree.y, color="black", size=0.2)
p
A1=p + coord_cartesian(clip = 'off') + scale_x_reverse() + coord_flip() +
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.2, "cm"),
        legend.text=element_text(size=2),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)
A1

y_axis <- tree.y %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

#--------------------------Run clustering for x axis------------------------------------------
##transpose data matrix
tree.matrix3 = t(tree.matrix2) %>% as.data.frame

#create distance matrix
dd <- dist(scale(tree.matrix3), method = "euclidean")

#Hierarchical Clustering
tree.x <- hclust(dd, method = "ward.D2")

#view clustering with ggtree
p=ggtree(tree.x, color="black", size=0.2)
p
A3=p + coord_cartesian(clip = 'off') + #scale_x_reverse() + coord_flip() +
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.2, "cm"),
        legend.text=element_text(size=2)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE) +
  theme(plot.margin = unit(c(0,0,0,10), "cm"))
A3

x_axis <- tree.x %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

#----------------------------------order x and y axes---------------------------------------------------
#order the x axis
#Res_fer$Contrast <- factor(Res_fer$Contrast, levels = x_axis)
Res_fer$Contrast <- factor(Res_fer$Contrast, levels = c('res4','res3','res2',"res1"))

#Reorder y axis IDs by clustering order (enter the list using the y_axis vector)
Res_fer$Id <- factor(Res_fer$Id, levels = y_axis)

#-----------------------make a dataframe with x and y coordinates of significant taxa to highlight-------
#sort heat map matrix based on the y_axis
Res_fer = Res_fer %>% arrange(factor(Id, levels = y_axis))

#make the list of cells to highlight as significant on y axis (note if there is NA in this column, this would cause problems)
Res_fer$x_position <- factor(Res_fer$Id)
Res_fer$x_position = as.numeric(Res_fer$x_position)

#make the list of cells to highlight as significant on x axis
Res_fer$y_position <- factor(Res_fer$Contrast)
Res_fer$y_position = as.numeric(Res_fer$y_position)

#make the list of cells to highlight as significant
frames = select(Res_fer, c('Id','Component_Name3','log2FoldChange','padj','x_position','y_position'))

# Extract significant results (padj < 0.05)
sig_values = subset(frames, padj < 0.05)

# filter significant results on min log2fold change (log2FoldChange > 2)
sig_values = subset(sig_values, abs(log2FoldChange) > 2)

#----------------------------------------------------------------------------------------------------
# plot heatmap
B=c("res1" = "SYNCOM vs No SYNCOM",
    "res2" = "Santander vs No SYNCOM",
    "res3" = "Huila vs No SYNCOM",
    "res4" = "Antioquia vs No SYNCOM"
)

##make y axis labels
##select columns
columns_to_keep = c("Id", "Component_Name3")
y_relabel = subset(sig_values, select=columns_to_keep)
y_relabel = y_relabel %>% distinct(Id, .keep_all = TRUE)
# Generating the desired output
output <- paste0('"', y_relabel$Id, '"="', y_relabel$Component_Name3, '"', collapse = ",\n")
# Print the output
cat(output)

C=c(
  "ASV12"="propan-2-yl acetate (ethereal)",
  "ASV33"="ethyl acetate (ethereal)",
  "ASV70"="butan-2-one (ethereal, acetone)",
  "ASV5"="3-methylbutyl propanoate (sweet, ripe tropical fruit)",
  "ASV62"="propan-2-ol (alcoholic)",
  "ASV21"="2-methylbut-3-en-2-ol (herbal)",
  "ASV6"="3,7-dimethyloctan-1-ol (floral, waxy, citrus green)",
  "ASV68"="2-methylpropan-1-ol (ethereal)",
  "ASV58"="ethanol (alcoholic)",
  "ASV64"="propan-1-ol (alcoholic)",
  "ASV65"="butan-1-ol (fermented)",
  "ASV14"="2-methylpyrazine (nutty)",
  "ASV16"="2-methylpropyl acetate (fruity)",
  "ASV66"="pentan-1-ol (fermented, fusel)",
  "ASV1"="ethenylbenzene (balsamic, sweet)",
  "ASV61"="hexanal (green)",
  "ASV28"="3-methylbutan-1-ol (fermented)",
  "ASV31"="nonanal (aldehydic)",
  "ASV76"="oxolan-2-one (creamy)",
  "ASV37"="Pentanoic acid, 2-methyl-, methyl ester (fruity)",
  "ASV73"="nonan-2-one (fruity)",
  "ASV20"="2,3,5,6-tetramethylpyrazine (nutty)",
  "ASV34"="2,3,5-trimethylpyrazine (nutty)",
  "ASV55"="pentan-2-yl acetate (herbal)",
  "ASV7"="ethyl octanoate (fruity, wine, waxy, sweet)",
  "ASV23"="3-methylbut-2-enyl acetate (fruity)",
  "ASV29"="ethyl hexanoate (fruity)",
  "ASV50"="2,3-dimethylpyrazine (nutty)",
  "ASV51"="heptan-2-yl acetate (brown)",
  "ASV3"="ethyl 2-phenylacetate (floral)",
  "ASV60"="4-methylpentanoic acid (cheesy)",
  "ASV44"="butane-2,3-diol (creamy, fruity)",
  "ASV4"="2-phenylethyl acetate (floral)",
  "ASV45"="3-hydroxybutan-2-one (buttery)",
  "ASV30"="3-methylbutyl acetate (fruity)",
  "ASV46"="2-amino-3-methylbutanoic acid",
  "ASV9"="butanoic acid (cheesy)",
  "ASV71"="propanoic acid (acidic)",
  "ASV8"="Propargyl alcohol (floral, geranium)",
  "ASV26"="2,5-dimethylpyrazine (chocolate, nutty)",
  "ASV43"="3-methylbutanoic acid (cheesy)",
  "ASV54"="5-methylfuran-2-carbaldehyde (caramellic)",
  "ASV19"="octan-1-ol (waxy)",
  "ASV10"="1-(1H-pyrrol-2-yl)ethanone (musty)",
  "ASV18"="2-(2-aminoethylamino)ethanol (mild ammonia)",
  "ASV42"="3-oxobutan-2-yl acetate (fruity, pungent)",
  "ASV41"="2-phenylbut-2-enal (sweet, honey, cocoa, nutty)",
  "ASV52"="5989-33-3 (earthy)"
)

#-----------------------------transform the data for plotting-------------------------------------
Res_fer2=Res_fer

# Set the Winsorizing percentiles (e.g., 1% and 99%)
winsorize_percentile_low <- 0.05
winsorize_percentile_high <- 0.95
# Calculate the Winsorized values
lower_limit <- quantile(Res_fer2$log2FoldChange, winsorize_percentile_low)
upper_limit <- quantile(Res_fer2$log2FoldChange, winsorize_percentile_high)
Res_fer2$log2FoldChange[Res_fer2$log2FoldChange < lower_limit] <- lower_limit
Res_fer2$log2FoldChange[Res_fer2$log2FoldChange > upper_limit] <- upper_limit

#----------------------------------------------------------------------------------------------------
p <- ggplot(data = Res_fer2) +
  geom_raster(aes(y = Contrast, x=Id, fill=log2FoldChange)) +
  scale_fill_gradient2(low="#02a0c5", high="#cc0001", na.value="grey", name="log2FoldChange") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.2)) +
  scale_x_discrete(label = C, expand = c(0,0)) + scale_y_discrete(label = B, expand = c(0,0), position = "right") +#label = B,
  geom_rect(data=sig_values, size=0.2, fill=NA, colour="black", aes(xmin=x_position - 0.45, xmax=x_position + 0.45, ymin=y_position - 0.45, ymax=y_position + 0.45))
p
A2 = p + guides(color="none", shape="none") +
  theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size=6, color="black"),
        axis.text.y=element_text(angle=0, hjust = 1, vjust = 0.5, size=6, color="black"),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.title = element_text(size=8, face="bold"),
        legend.text=element_text(size=4),
        legend.key.size = unit(0.3, "cm"),
        legend.margin=margin(c(0,0,0,0)),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(x = "Fermentation time (hr)") +
  guides(fill = guide_colourbar(title.position ='top',
                                frame.linewidth = 0.2,
                                ticks.colour = "black",
                                frame.colour ='black',
                                barwidth = 3,
                                barheight = 0.5,
                                #barheight = NULL,
                                title="Enrichment",
                                title.hjust =0.5,
                                color='red',
                                legend.position="bottom")) +
  theme(legend.title = element_text(size = 5, angle = 0, face="bold")) +
  theme(legend.margin=margin(c(0,0,0,0))) +
  theme(legend.position="bottom")
A2

#---------------------------------------------------------------------------------------------------------
##empty plot
E=ggplot() + theme_void()

## final plot merged
A1a = A1 + guides(color = "none", shape="none", fill="none")
A2a = A2 + guides(color = "none", shape="none", fill="none")
A3a = A3 + guides(color = "none", shape="none", fill="none")

##plot figure
Plot2 = egg::ggarrange(A1a,A2a, nrow = 2, heights = c(0.5,2), ncol = 1, widths = c(1))
Plot2
ggsave("Fig5h.pdf", height=6, width=13, units='cm', plot = Plot2)

## extract legend option 1
library(ggpubr)
library(egg)
# Extract the legend. Returns a gtable
A2_L <- get_legend(A2)
# Convert extracted legend to a ggplot
A2_L = as_ggplot(A2_L)
##plot
legend = egg::ggarrange(A2_L,
                        ncol = 1, widths = c(1), heights=c(1))
ggsave("Fig5h_legend.pdf", height=2, width=6, units='in', plot=legend)
