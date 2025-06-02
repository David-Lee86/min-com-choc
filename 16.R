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
library(readxl)
library(dplyr)
library(data.table)
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
#Plot Variance
LC = read.csv("Colombia_flavour_PERMANOVA.csv")
data_melt <- melt(LC)

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)
data_melt$Effect <- factor(data_melt$Effect, levels=unique(data_melt$Effect))

mypal = c('white','#b33a89','#e56281','#ffc25d','#d2df6d','#ede1e1','#53c0cc','#1f4f58')

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
# Import the saved data
df <- read.csv("Colombia_flavour_ordination_data.csv")

#Reorder samples
df$Location<- factor(df$Location,levels = c('Santander','Huila', 'Antioquia'))

mypal = c("Santander" = "#ffaabb",
          "Huila" = "#bbcc33",
          "Antioquia" = "#44bb99")

# Rebuild your custom plot
A2 <- ggplot(df, aes(x = CAP1, y = CAP2)) + 
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed") +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed") +
  geom_point(aes(fill=Location, shape=Location), size = 7, shape = 21, colour = "black", stroke=0.7) +
  scale_fill_manual(values = mypal) +
  theme_classic() +
  guides(color = FALSE, shape = FALSE, fill = FALSE) +
  theme(
    axis.text.y   = element_text(angle=0, hjust = 0, vjust = 0.5, size=14, color = "black"),
    axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=14, color = "black"),
    axis.title.y  = element_text(size=14, face="bold"),
    axis.title.x  = element_text(size=14, face="bold"),
    axis.ticks = element_line(colour = "black", size = 0.35),
    legend.text = element_text(size=12),
    legend.key.size = unit(0.4, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.1),
    plot.margin = margin(0, 0, 0, 0),
    panel.border = element_rect(colour = "black", fill=NA, size = 0.35),
    plot.title = element_text(hjust = 0.5, face="bold")
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 6), title = "Location", nrow = 3)) +
  annotate("text", size=4, label = expression(italic(R)^2 == 0.28342), x = -0.7, y = 1.2) +
  annotate("text", size=4, label = expression(italic(p) < 0.0001), x = -0.7, y = 1) +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank(), legend.margin = margin(c(0, 0, 0, 0)))
A2
# -------------------------------------------------------------------------------------------------
##empty plot
A4=ggplot() + theme_void()

### final plot merged
A1a = A1 + guides(color = "none", shape="none", fill="none")
A2a = A2 + guides(color = "none", shape="none", fill="none")

Plot5 = egg::ggarrange(A1a,A4,A2a,
                       ncol = 3, widths = c(0.1,0.05,1.2), heights=c(1))
ggsave("Fig3a.pdf", height=3.5, width=5, units='in', plot = Plot5)


### extract legend option 1
library(ggpubr)
library(egg)

# Extract the legend. Returns a gtable
A2_L <- get_legend(A2)
A1_L <- get_legend(A1)

# Convert extracted legend to a ggplot
A2_L = as_ggplot(A2_L)
A1_L = as_ggplot(A1_L)
##plot
legend = egg::ggarrange(A2_L,A1_L,
                        ncol = 2, widths = c(1,1), heights=c(1))
ggsave("Fig3a_legend.pdf", height=3.5, width=8, units='in', plot = legend)



