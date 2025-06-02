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
library("BiocParallel")
library("readxl")
library("reshape")
library(dplyr)

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

# -------------------------------------------------------------------------------------------------
# Import metadata from excel file- Linux
LC = read.csv("colombia_structure_result.csv")
class(LC)

# -------------------------------------------------------------------------------------------------
# Reorder Data - Convert to Long Format
data_melt <- melt(LC)

# Reorder x-axis to match the order in the data
data_melt$Sample_name <- as.character(data_melt$Sample_name)
data_melt$Sample_name <- factor(data_melt$Sample_name, levels=unique(data_melt$Sample_name))

# Define Color Palette for Cacao Accessions
mypal = c(
  'Amelonado' = "#55003b",
  'Contamana' = '#c1204c',
  'Criollo' = '#a0ae00',
  'Curaray' = '#e0e004',
  'Guiana' = '#508104',
  'Iquitos' = '#535121',
  'Marañon' = '#e0aa02',
  'Nacional' = '#f35828',
  'Nanay' = '#8e67a0',
  'Purús' = "#cfc097",
  'Amelonado-Criollo' = '#fc0558',
  'Santander' = '#a2d79f',
  'Huila' = "#7aa4a2",
  'Antioquia' = "#05788a"
)

# -------------------------------------------------------------------------------------------------
# Reorder y-axis and Farms
data_melt$variable <- factor(data_melt$variable, levels = c(
  'Amelonado', 'Contamana', 'Criollo', 'Curaray', 'Guiana', 'Iquitos', 'Marañon', 
  'Nacional', 'Nanay', 'Purús'
))

data_melt$Population <- factor(data_melt$Population, levels = c("Reference", "Amelonado-Criollo", "Santander", "Huila", "Antioquia"))

# New facet label names
B <- list(
  'Reference' = "Reference",
  'Amelonado-Criollo' = "Amelonado\nCriollo",
  'Santander' = "Santander",
  'Huila' = "Huila",
  'Antioquia' = "Antioquia"
)

# Create labeller function for facets
hospital_labeller <- function(variable, value) {
  return(B[value])
}

library(ggh4x)

# Plot
A1 = ggplot(data_melt, aes(fill = variable, y = value, x = Sample_name)) +
  geom_bar(aes(fill = variable, color = variable), position = "fill", stat = "identity", width = 0.99, size = 0) +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(color = "none", shape = "none") +
  scale_fill_manual(values = mypal) + scale_color_manual(values = mypal) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 7, color = "black"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 7, face = "bold"),
    axis.title.x = element_text(size = 7, face = "bold"),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 7, face = c("plain", "italic")),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.line = element_blank(),
    panel.border = element_blank(),
    panel.spacing = unit(0.2, "cm")
  ) +
  labs(x = "Cacao accessions", y = "Membership (K=10)") +
  facet_wrap(~Population, scales = "free_x", labeller = hospital_labeller, nrow = 1) +
  force_panelsizes(cols = c(0.7, 0.3, 0.5, 0.2, 0.29)) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 4))) +
  theme(
    strip.text = element_text(size = 7, face = "bold", margin = margin(t = 0.5, b = 0.5)),
    panel.spacing = unit(0.2, "cm"),
    strip.placement = "outside"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

A1

# -------------------------------------------------------------------------------------------------
# Final Plot Merged
A1a = A1 + guides(color = "none", shape = "none", fill = "none")
A1a

# Save the plot
ggsave("Fig2a.pdf", height = 6, width = 9, units = 'cm', plot = A1a)



