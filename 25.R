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
library(knitr)
library(rmdformats)
library(ggpubr)
library(data.table)
library(scales)
library(BiocParallel)
library(ggvenn)

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))
# -------------------------------------------------------------------------------------------------

##import data
LC = read.csv("venn_diagram_community_metabolites.csv")

# -------------------------------------------------------------------------------------------------
##format data for plotting
# Convert dataframe to a list
LC <- as.list(LC)

# Remove empty strings from all elements in the list
LC <- lapply(LC, function(x) x[nzchar(x)])

# Count the number of strings from each element in the list
lengths(LC)

# -------------------------------------------------------------------------------------------------
##plot
ggvenn(LC)

# Change the fill color
names(LC) <- c("Full community","Defined community","SYNCOM")
A1 = ggvenn(LC, text_size=2.5, text_color = "black",
  #fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  fill_color = c("#0073C2FF", "#f19d06", "#7a972b"),
  fill_alpha = 0.7,
  stroke_size = 0.5, set_name_size = 3.7
)
A1
ggsave("Fig4c.pdf", height=3, width=4, units='in', plot = A1)



