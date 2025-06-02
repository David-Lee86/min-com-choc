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
library(ggpubr)
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
##Import data
total = read.csv("SYNCOM_GCMS_PCA.csv")

B=c('Antioquia' = 'Antioquia',
    'Huila' = 'Huila',
    "Santander" = "Santander",
    "Ghana" = "Ghana",
    "SYNCOM" = "SYNCOM",
    "No_SYNCOM" = "No SYNCOM"
    )

#adjust sample order
total$Description4<- factor(total$Description4, levels = c("No_SYNCOM","SYNCOM","Santander",'Huila','Antioquia'))

mypal = c('#b579ad','#c8b2e0',"#f1daec",'#97a341','#ccd278')

# Manually enter variance explained values
var_exp_PC1 <- 0.56
var_exp_PC2 <- 0.20

# plot
A1 = ggplot(total, aes(x = PC1, y = PC2)) +
  geom_point() + 
  labs(
    x = paste0("Axis 1 (", scales::percent(var_exp_PC1), ")"),
    y = paste0("Axis 2 (", scales::percent(var_exp_PC2), ")")
  ) +
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed") +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed") +
  geom_point(aes(fill=Description4, color=Description4), shape = 21, size = 6, colour = "black",stroke=0.6) +
  scale_fill_manual(labels = B, values = mypal) +#labels = B
  theme_classic() + guides(color=FALSE, shape=FALSE) + guides(fill = guide_legend(title="Sample")) +
  theme(axis.text.y = element_text(angle=0, hjust = 0, vjust = 0.5, size=12, color = "black"),
        axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 1, size=12, color = "black"),
        axis.title.y = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        legend.key.size = unit(0.07, "cm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.35), plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(legend.position="bottom") +
  theme(plot.margin=unit(c(0.1,1,0.1,0.1),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0))) +
  guides(fill=guide_legend(nrow=5, override.aes = list(size = 6))) +
  scale_x_discrete(expand = c(0.1, 0.1)) + scale_y_discrete(expand = c(0.1, 0.1))
A1

##plot
A1a = A1 + guides(color = "none", shape="none", fill="none")
ggsave("Fig5g.pdf", height=2, width=2.75, units='in', plot = A1a)


### extract legend option 1
library(ggpubr)
library(egg)
# Extract the legend. Returns a gtable
A1_L <- get_legend(A1)
# Convert extracted legend to a ggplot
A1_L = as_ggplot(A1_L)
A1_L
ggsave("Fig5g_legend.pdf", height=3, width=3, units='in', plot = A1_L)
