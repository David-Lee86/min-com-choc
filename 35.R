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
library(readxl)
library(dplyr)
library(plyr)
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
#Import data
LC <- read.csv("SYNCOM_sensory.csv", check.names = FALSE)

# -------------------Rescaling---------------------------------------------------------------------
# Reorder data -convert to long format
data_melt <- melt(LC)

#change column names of data_melt
colnames(data_melt) <- c("X1", "X2", "value")

# Rescaling
data_melt <- ddply(data_melt, .(X2), transform, rescale = rescale(value,to=c(0,6)))

# Reorder data -convert back to short format (unmelt)
heatmap.dataset2 = dcast(data_melt, X1~X2, value.var="rescale") %>% as.data.frame

# -------------------------------------------------------------------------------------------------
#Format the dataframe for clustering
library(tidyverse)
heatmap.dataset2 = heatmap.dataset2 %>% remove_rownames %>% column_to_rownames(var="X1")

#Format the dataframe for clustering 
tree.matrix2 <- heatmap.dataset2

# -----------------Run clustering for y_axis--------------------------------------------------------
# create distance matrix
dd <- dist(scale(tree.matrix2), method = "euclidean")

# run Hierarchical Clustering of distance matrix
tree.y <- hclust(dd, method = "single")

#plot clustering with ggtree
p=ggtree(tree.y, color="black", size=0.3)
A1=p +
  coord_cartesian(clip = 'off') + #scale_x_reverse() + coord_flip() +
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.2, "cm"),
        legend.text=element_text(size=2)) +
    guides(color=FALSE, shape=FALSE, fill=FALSE) +
  theme(plot.margin = unit(c(0,0,0,10), "cm"))
  #geom_tiplab()
A1

y_axis <- tree.y %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

# --------------------Run clustering for x_axis-----------------------------------------------------------------------------
#transpose the dataframe
tree.matrix3=t(tree.matrix2)

# create distance matrix
dd <- dist(scale(tree.matrix3), method = "euclidean")

# run Hierarchical Clustering of distance matrix
tree.x <- hclust(dd, method = "single")

#plot clstering with ggtree
p=ggtree(tree.x, color="black", size=0.3)
A3=p +
  coord_cartesian(clip = 'off') + scale_x_reverse() + coord_flip() +
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 2),
        legend.key.size = unit(0.2, "cm"),
        legend.text=element_text(size=2),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)
A3

x_axis <- tree.x %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

# --------------------create annotation bar---------------------------------------------------------
#Import annotations
metadata <- read_excel("chocolate_flavour_groupings.xlsx", sheet = "flavour_groups_SYNCOM")

#Reorder x_axis by flavor clustering (x axis)
metadata$Attribute <- factor(metadata$Attribute, levels = x_axis)

mypal=c('#eb1a3c','#e84c8b','#f69997','#f59535','#ffda69','#e0e14f','#a3c049',
        '#85ccbc','#2a9cd1','#bc3e85','#52534e','#dddad5','#7b49ac')

A5 <- ggplot(metadata, aes(Attribute, y=1, fill=Flavour_group)) + geom_tile() +
  scale_fill_manual(values = mypal, labels = c("total_acid" = "total acid",
                                               "astringency" = "astringency",
                                               "bitterness" = "bitterness",
                                               "browned_fruit" = "browned fruit",
                                               "cocoa" = "cocoa",
                                               "floral" = "floral",
                                               "fresh_fruit" = "fresh fruit",
                                               "nutty" = "nutty",
                                               "off-flavours" = "off-flavours",
                                               "roast" = "roast",
                                               "spice" = "spice",
                                               "sweet" = "sweet",
                                               "woody" = "woody"
)) +
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

#Reorder x_axis by flavor clustering
Res_fer$variable <- factor(Res_fer$variable, levels = x_axis)

#Reorder y_axis by sample clustering
Res_fer$Id <- factor(Res_fer$Id, levels = y_axis)

# -------------------------------------------------------------------------------------------------
### rename y values
y_axis
output <- paste0('"', y_axis, '" = "', gsub("_", " ", y_axis), '"', collapse = ",\n")
cat(output)
C=c(
  "Antioquia" = "Antioquia",
  "No_SYNCOM" = "No SYNCOM",
  "Ghana" = "Ghana",
  "Ivory_Coast" = "Ivory Coast",
  "SYNCOM" = "SYNCOM",
  "Santander" = "Santander",
  "Huila" = "Huila",
  "Madagascar" = "Madagascar"
)

### rename x values
x_axis
output <- paste0('"', x_axis, '" = "', gsub("_", " ", x_axis), '"', collapse = ",\n")
cat(output)
B=c(
  "floral_earthy,_mushroom,_moss,_woodsy" = "EMMW",
  "spice_tobacco" = "tobacco",
  "wood_darkwood" = "dark wood",
  "cocoa" = "cocoa",
  "roastdegree" = "roast",
  "astringency" = "astringency",
  "bitterness" = "bitterness",
  "floral_grassy,_greenvegetal,_herbal" = "GGVH",
  "spice_savory,_umami" = "savory/umami",
  "nutty_nutflesh" = "nut flesh",
  "spice_spices" = "spices",
  "fruit_overripe" = "overripe fruit",
  "over_fermented,_rottenfruit" = "OFRF",
  "acid_lactic" = "lactic",
  "putrid,_manure" = "putrid/manure",
  "dirty,_dusty" = "dirty/dusty",
  "smoky" = "smoky",
  "wood_resin" = "resin",
  "musty" = "musty",
  "wood_lightwood" = "ligh twood",
  "floral_orangeblossom" = "orange blossom",
  "fruit_dried" = "dried fruit",
  "nutty_nutskins" = "nut skin",
  "caramel_panela" = "caramel/panela",
  "acid_fruit" = "fruit acid",
  "fruit_tropical" = "tropical friut",
  "fruit_brown" = "brown fruit",
  "acid_acetic" = "acetic",
  "fruit_dark" = "dark fruit",
  "fruit_yellow,_orange_whiteflesh" = "YOWF",
  "floral_flowers" = "flowers",
  "fruit_berry" = "berry",
  "fruit_citrus" = "citrus"
)

# -------------------------------------------------------------------------------------------------
library("hrbrthemes")
p = ggplot(data=Res_fer) +
  geom_raster(aes(x=variable, y=Id, fill=value)) +
  scale_fill_distiller(palette = "RdPu",direction = 1) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_y_discrete(position = "right", label = C, expand = c(0,0)) + scale_x_discrete(label = B, expand = c(0,0))
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
  theme(legend.position="bottom")
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
ggsave("Fig5i.pdf", height=2.5, width=6, units='in')


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
ggsave("Fig5i_legend.pdf", height=2, width=6, units='in', plot=legend)
