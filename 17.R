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
library(ggplot2)
library(reshape2)
library(scales)
library(readxl)
library(forcats)
library(RColorBrewer)
library(ggrepel)
library(ggdendro)
library(ggtree)
library(dplyr)
library(BiocParallel)
library(tidyverse)
# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

# -------------------------------------------------------------------------------------------------
#Import data as .csv
LC<- read_excel("Colombia_flavor_heatmap.xlsx", sheet = "heatmap")

# -------------------------------------------------------------------------------------------------
names(LC)

#attribute sums
colSums(LC %>% select(where(is.numeric), -1))

#subset if needed
LC.sub = LC

LC.sub = LC.sub %>% remove_rownames %>% column_to_rownames(var="Name")

# Set as a data matric
LC.sub=data.matrix(LC.sub)

# -------------------------------------------------------------------------------------------------
# Rescale data per column and reshape for heatmap
heatmap.dataset2 <- as.data.frame(LC.sub) %>%
  rownames_to_column("X1") %>%
  pivot_longer(-X1, names_to = "X2", values_to = "value") %>%
  group_by(X2) %>%
  mutate(rescale = rescale(value, to = c(0, 6))) %>%
  ungroup() %>%
  select(X1, X2, rescale) %>%
  pivot_wider(names_from = X2, values_from = rescale)

# -------------------------------------------------------------------------------------------------
#Format the dataframe for clustering
heatmap.dataset2 = heatmap.dataset2 %>% remove_rownames %>% column_to_rownames(var="X1")

#Format the dataframe for clustering 
tree.matrix2 <- heatmap.dataset2

# -------------------------------------------------------------------------------------------------
#Run clustering for samples
# create distance matrix
dd <- dist(scale(tree.matrix2), method = "euclidean")

# run Hierarchical Clustering of distance matrix
tree <- hclust(dd, method = "ward.D2")

#plot clustering with ggtree
p=ggtree(tree, color="black", size=0.3)
p
A1=p +
  coord_cartesian(clip = 'off') +
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.2, "cm"),
        legend.text=element_text(size=2)) +
    guides(color=FALSE, shape=FALSE, fill=FALSE) +
  theme(plot.margin = unit(c(0,0,0,10), "cm"))
A1

# Get the order of tip labels from the tree
y_axis <- tree %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

# -------------------------------------------------------------------------------------------------
#Run clustering for flavours
tree.matrix3=t(tree.matrix2)

# create distance matrix
dd <- dist(scale(tree.matrix3), method = "euclidean")

# run Hierarchical Clustering of distance matrix
tree <- hclust(dd, method = "ward.D2")

#plot clstering with ggtree
p=ggtree(tree, color="black", size=0.3)
p
A3=p +
  coord_cartesian(clip = 'off') + scale_x_reverse() + coord_flip() +
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.2, "cm"),
        legend.text=element_text(size=2),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)
A3

# Get the order of tip labels from the tree
x_axis <- tree %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

# -------------------------------------------------------------------------------------------------
## get clustering groups:
gr <- cutree(tree, k = 2) %>% as.data.frame
colnames(gr) <- c("cls")
gr <- tibble::rownames_to_column(gr, "variable")#create a column with the row names, call it "Id"
gr = mutate(gr, grp1 = case_when(cls==1~"C1",cls==2~"C2",cls==3~"C3",cls==4~"C4",cls==5~"C5",
                                 cls==6~"C6",cls==7~"C7",cls==8~"C8",cls==9~"C9",cls==10~"C10")) %>% select (-cls)
# -------------------------------------------------------------------------------------------------
#colour attributes by groups
metadata <- read_excel("Colombia_flavor_heatmap.xlsx", sheet = "groupings")

#Reorder flavor_order by flavor clustering (x axis)
metadata$Attribute <- factor(metadata$Attribute, levels = x_axis)

mypal=c('#eb1a3c','#e84c8b','#f69997','#f59535','#ffda69','#e0e14f','#a3c049',
        '#85ccbc','#2a9cd1','#bc3e85','#52534e','#dddad5','#7b49ac')

A5 <- ggplot(metadata, aes(Attribute, y=1, fill=Flavour_group)) + geom_tile() +
  #pc <- ggplot(MAA, aes(Attribute, y=anno, fill=type)) + geom_tile() + 
  scale_fill_manual(values = mypal, labels = c("astringency",
                                               "bitterness",
                                               "browned fruit",
                                               "cocoa",
                                               "floral",
                                               "fresh fruit",
                                               "nutty",
                                               "off-flavours",
                                               "roast",
                                               "spice",
                                               "sweet",
                                               "total acid",
                                               "woody")) +
  coord_cartesian(clip = 'off') + #scale_x_reverse() + coord_flip() +
  theme(axis.text.x=element_blank(),
        #axis.text.x=element_text(angle=45, hjust = 1, vjust = 1, size=7, color="black"),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(size = 10,colour='black'),
        legend.position="bottom",
        legend.key.size = unit(0.2, "cm")) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"),
        legend.margin=margin(c(0,0,0,0))) +
  guides(fill=guide_legend(ncol=2, override.aes = list(size = 4)))

A5
# -------------------------------------------------------------------------------------------------
# Finalize dataframe for plotting
heatmap.dataset2 <- tibble::rownames_to_column(heatmap.dataset2, "Id")
rownames(heatmap.dataset2) = heatmap.dataset2$Id

#Format the dataframe for ggplot
colnames(heatmap.dataset2)
Res_fer <- melt(heatmap.dataset2)

# add the clustering of variables
Res_fer <- merge(Res_fer, gr, by="variable")

# -------------------------------------------------------------------------------------------------
#order data
#Reorder flavor_order by flavor clustering (x axis)
Res_fer$variable <- factor(Res_fer$variable, levels = x_axis)

#Reorder y_axis by sample clustering (y axis)
Res_fer$Id <- factor(Res_fer$Id, levels = y_axis)

# -------------------------------------------------------------------------------------------------
### rename x axis
x_axis
output <- paste0('"', paste(x_axis, collapse = '", "'), '"')
cat(output)
#remove replace values in string
output=gsub("_", " ", output)
cat(output)
B=c("EMMW", "dark wood", "roast", "cocoa", "tobacco", "resin", "putrid/manure", "overripe",
    "spices", "OFRF", "smoky", "lactic", "dirty/dusty", "musty", "fruit", "brown", "dried",
    "bitterness", "astringency", "nut flesh", "orange blossom", "savory/umami", "GGVH", "dark",
    "nut skins", "berry", "flowers", "caramel/panela", "YOWF", "light wood", "citrus", "acetic", "tropical")

### rename y axis
y_axis
output <- paste0('"', paste(y_axis, collapse = '", "'), '"')
cat(output)
#remove replace values in string
output=gsub("_", " ", output)
cat(output)
C=c("Ghana", "Ivory Coast", "Antioquia harvest 1", "Antioquia harvest 2", "Santander harvest 1",
    "Santander harvest 2", "Madagascar", "Huila harvest 1", "Huila harvest 2")
# -------------------------------------------------------------------------------------------------
library("hrbrthemes")
p = ggplot(data=Res_fer) +
  geom_raster(aes(x=variable, y=Id, fill=value)) +
  scale_fill_distiller(palette = "RdPu",direction = 1) + #theme_ipsum() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_y_discrete(label = C, position = "right", expand = c(0,0)) +
  scale_x_discrete(label = B, expand = c(0,0))
  #scale_x_discrete(label = function(x) stringr::str_trunc(x, 20))
p
A2 = p + theme(axis.text.x=element_text(angle=90, hjust = 1, vjust = 0.5, size=10, color="black"),
               axis.text.y=element_text(angle=0, hjust = 1, vjust = 0.5, size=10, color="black"),
               panel.background = element_blank(),
               axis.ticks = element_blank(),
               plot.title = element_blank(),
               axis.title.x  = element_blank(),
               axis.title.y  = element_blank(),
               legend.title = element_text(size = 10,colour="black"),
               legend.text= element_text(size = 10,colour='black'),
               legend.key.size = unit(0.5, "cm")) +
  theme(plot.margin = unit(c(0,0,0,4), "cm")) +
  guides(fill = guide_colourbar(title.position ='top',
                                frame.linewidth = 0.3,
                                ticks.colour = "black",
                                frame.colour ='black',
                                barwidth = 6,
                                #barheight = NULL,
                                title="Attribute",
                                title.hjust =0.5,
                                color='red',
                                legend.position="bottom")) +
  theme(legend.margin=margin(c(0,0,0,0))) +
  #guides(color="none", shape="none", fill="none") +
  theme(legend.position="bottom")
  #facet_grid(~grp1, scales="free_x")
  #+ facet_wrap(~grp1, scales="free_x")
A2

# -------------------------------------------------------------------------------------------------
##empty plot
E=ggplot() + theme_void()

### final plot merged
A1a = A1 + guides(color = "none", shape="none", fill="none")
A2a = A2 + guides(color = "none", shape="none", fill="none")
A3a = A3 + guides(color = "none", shape="none", fill="none")
A5a = A5 + guides(color = "none", shape="none", fill="none")

#plot using aplot
library(aplot)
Plot2 = A2a %>%
  insert_left(A1a, width=0.04) %>%
  insert_top(A5a, height=0.08) %>%
  insert_top(A3a, height=0.08) %>%
  insert_top(E, height=0.03) %>%
  insert_left(E, width=0.2)
Plot2
ggsave("Fig3b.pdf", height=3, width=8, units='in', plot = Plot2)

### extract legend option 1
library(ggpubr)
library(egg)

# Extract the legend. Returns a gtable
A2_L <- get_legend(A2)
A5_L <- get_legend(A5)

# Convert extracted legend to a ggplot
A2_L = as_ggplot(A2_L)
A5_L = as_ggplot(A5_L)
##plot
legend = egg::ggarrange(A2_L,A5_L,
                       ncol = 2, widths = c(1,1), heights=c(1))
ggsave("Fig3b_legend.pdf", height=4, width=8, units='in', plot = legend)

