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
library(ggplot2)
library(readxl)
library(dplyr)
library(BiocParallel)
library(tidyverse)
library(ggtree)
library(gridExtra)
library(grid)
library(scales)
library(hrbrthemes)
#-------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))
#-------------------------------------------------------------------------------------------------
# Import data
heatmap.dataset2 = read.csv("Colombia_RF_Inflection_point.csv")

#-------------------------------------------------------------------------------------------------
# Set up row names to X
row.names(heatmap.dataset2)=heatmap.dataset2$Feature

#Format the dataframe for clustering
tree.matrix2 <- heatmap.dataset2 %>% select (-Feature)

# ------------------------y axis-------------------------------------------------------------------------
#create distance matrix
dd <- dist(scale(tree.matrix2), method = "euclidean")

#Hierarchical Clustering
tree.y <- hclust(dd, method = "ward.D2")

#view clustering with ggtree
p=ggtree(tree.y, color="black", size=0.2)
p + coord_cartesian(clip = 'off') +
  theme(panel.background = element_blank(),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)

# Get the order of tip labels from the tree
y_axis <- tree.y %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

# -----------------------x axis--------------------------------------------------------------------------
## transpose
tree.matrix3 = t(tree.matrix2)

#create distance matrix
dd <- dist(scale(tree.matrix3), method = "euclidean")

#Hierarchical Clustering
tree.x <- hclust(dd, method = "ward.D2")

#view clustering with ggtree
p=ggtree(tree.x, color="black", size=0.2)
p + coord_cartesian(clip = 'off') +
  theme(panel.background = element_blank(),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)

# Get the order of tip labels from the tree
x_axis <- tree.x %>%
  fortify() %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label) %>%
  unique()

#-------------------------------------------------------------------------------------------------
library(reshape)
#Reorder data -convert to long format
Res_fer <- melt(heatmap.dataset2)
colnames(Res_fer) <- c("Id", "variable","value")

#Reorder y axis IDs by clustering order (enter the list using the y_axis vector)
Res_fer$Id <- factor(Res_fer$Id, levels = y_axis)

#Reorder x axis IDs by clustering order (enter the list using the x_axis vector)
Res_fer$variable <- factor(Res_fer$variable, levels = x_axis)

#-------------------------------------------------------------------------------------------------
### rename y values
y_axis
#C = output=gsub("_", " ", y_axis)
#cat(C)

### rename x values
x_axis
B = output=gsub("_", " ", x_axis)
cat(B)

### rename x values extra
flavour_names = x_axis %>% as.data.frame
colnames(flavour_names) <- c("First_label")
AC<- read_excel("chocolate_flavour_groupings.xlsx", sheet = "flavour_groups")
flavour_names <- merge(flavour_names,AC,by="First_label")
flavour_names = flavour_names %>% arrange(factor(First_label, levels = x_axis))
AC2 = flavour_names %>% pull(Attribute)
B = output=gsub("_", " ", AC2)
cat(B)

# -------------------------------------------------------------------------------------------------
## highlight taxa
taxa = heatmap.dataset2 %>% select(c('Feature'))
taxa = mutate(taxa, col2 = case_when(Feature=='Saccharomyces'~"#fc0558", Feature=='Torulaspora'~"#fc0558"))
taxa[is.na(taxa)] <- "black"
taxa = taxa %>% arrange(factor(Feature, levels = y_axis))
taxa_colour = taxa$col2
taxa_colour

# -------------------------------------------------------------------------------------------------
## plot
p = ggplot(data=Res_fer) +
  geom_raster(aes(x=variable, y=Id, fill=value)) +
  scale_fill_gradient2(low="#ff6b12", high="#0d7595", na.value="grey", name="Feature importance") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_y_discrete(position = "left", expand = c(0,0)) + scale_x_discrete(label = B, expand = c(0,0)) #+ #label = B #labels=c("Saccharomyces"=expression(bold(Saccharomyces)), parse=TRUE)
A2 = p +
  theme(axis.text.x=element_text(angle=35, hjust = 1, vjust = 1, size=5.2, color="black", margin=margin(-1,0,0,0)),
        axis.text.y=element_text(angle=0, hjust = 1, vjust = 0.5, size=5.2, color=taxa_colour, margin=margin(0,-1,0,0)),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.title.x  = element_blank(),
        axis.title.y  = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  labs(x = "Fermentation time (hr)") +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 3.5,
                                direction = "vertical",
                                frame.linewidth = 0.2,
                                ticks = TRUE,
                                title.position = 'left',
                                frame.colour = "black",
                                ticks.colour = "black")) +
  theme(legend.margin=margin(c(0,-5,0,-10)),
        legend.title = element_text(size=5.2, face="bold", hjust = 0.5, angle=90),
        legend.text=element_text(size=5.2),
        #legend.title.align = 90,
        legend.direction = "vertical")
A2

#view clustering with ggtree
p=ggtree(tree.y, color="black", size=0.12)
A1 = p + coord_cartesian(clip = 'off') +
  theme(panel.background = element_blank(),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)
A1

#plot clstering for flavor
p=ggtree(tree.x, color="black", size=0.12)
A3= p + coord_cartesian(clip = 'off') + coord_flip() +
  theme(panel.background = element_blank(),
        plot.margin=margin(0, 0, 0, 0)) +
  guides(color=FALSE, shape=FALSE, fill=FALSE)
A3 = A3 + scale_x_reverse() 
A3

#colour attributes by flavour groups
x_axis
#Reorder x axis factors by x_axis
flavour_names$First_label <- factor(flavour_names$First_label, levels = x_axis)

#replace name of values in column
library(tidyverse)
flavour_names$Flavour_group <- str_replace(flavour_names$Flavour_group, "_", " ")

mypal=c("astringency"='#eb1a3c',
        "bitterness"='#e84c8b',
        "browned fruit"='#f69997',
        "cocoa"='#f59535',
        "floral"='#ffda69',
        "fresh fruit"='#e0e14f',
        "nutty"='#a3c049',
        "off-flavours"='#85ccbc',
        "roast"='#2a9cd1',
        "spice"='#bc3e85',
        "sweet"='#52534e',
        "total acid"='#dddad5',
        "woody"='#7b49ac')

A5 <- ggplot(flavour_names, aes(First_label, y=1, fill=Flavour_group)) + geom_tile() +
  scale_fill_manual(values = mypal) +
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
  guides(fill=guide_legend(ncol=2, override.aes = list(size = 4))) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
A5

#colour attributes by quality groups
x_axis
#Reorder x axis factors by x_axis
flavour_names$First_label <- factor(flavour_names$First_label, levels = x_axis)

mypal2=c('bulk'='#a59579',
         'fine/flavour'='#f6f0e7',
         'undesirable'='#52534e')

A6 <- ggplot(flavour_names, aes(First_label, y=1, fill=Quality_group)) + geom_tile() +
  scale_fill_manual(values = mypal2) +
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
  guides(fill=guide_legend(ncol=2, override.aes = list(size = 4))) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0))
A6

#######################################################################################################
##empty plot
E=ggplot() + theme_void()

### final plot merged
A1a = A1 + guides(color = "none", shape="none", fill="none")
A2a = A2 #+ guides(color = "none", shape="none", fill="none")
A3a = A3 + guides(color = "none", shape="none", fill="none")
A5a = A5 + guides(color = "none", shape="none", fill="none")
A6a = A6 + guides(color = "none", shape="none", fill="none")

#plot using aplot
library(aplot)
Plot2 = A2a %>%
  insert_left(A1a, width=0.02) %>%
  insert_top(A5a, height=0.03) %>%
  insert_top(A6a, height=0.03) %>%
  insert_top(A3a, height=0.04) %>%
  insert_top(E, height=0.03) %>%
  insert_left(E, width=0.2)
Plot2
ggsave("Fig3d_iv.pdf", height=1.8, width=4.2, units='in', plot = Plot2)



### extract legend option 1
library(ggpubr)
library(egg)

# Extract the legend. Returns a gtable
A2_L <- get_legend(A2)
A5_L <- get_legend(A5)
A6_L <- get_legend(A6)

# Convert extracted legend to a ggplot
A2_L = as_ggplot(A2_L)
A5_L = as_ggplot(A5_L)
A6_L = as_ggplot(A6_L)
##plot
legend = egg::ggarrange(A2_L,A5_L,A6_L,
                        ncol = 3, widths = c(1,1,1), heights=c(1))
ggsave("Fig3d_iv_legend.pdf", height=4, width=8, units='in', plot = legend)

