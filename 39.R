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
# Load required libraries
library("ggplot2")
library("readxl")
library("dplyr")
library(data.table)
library(tidyverse)
library("BiocParallel")

# -------------------------------------------------------------------------------------------------
# Set working directory automatically (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# -------------------------------------------------------------------------------------------------
# Import data
rl_filtered = read.csv("colombia_read_length_data.csv")

# -------------------------------------------------------------------------------------------------
colnames(rl_filtered)[which(names(rl_filtered) == "Length.bp.")] <- "Length(bp)"

# Calculate row sums excluding the 'Length(bp)' column
rl_filtered$Frequency <- rowSums(rl_filtered[, !names(rl_filtered) %in% "Length(bp)"], na.rm = TRUE)

##select columns
columns_to_keep = c("Length(bp)", "Frequency")
rl_sub = subset(rl_filtered, select=columns_to_keep)

# -------------------------------------------------------------------------------------------------
# Calculate the weighted sum and total frequency
weighted_sum <- sum(rl_sub$`Length(bp)` * rl_sub$Frequency)
total_frequency <- sum(rl_sub$Frequency)

# Calculate the average read length
average_length <- weighted_sum / total_frequency
average_length

# -------------------------------------------------------------------------------------------------
p <- ggplot(rl_sub, aes(x = `Length(bp)`, y = Frequency)) +
  geom_line(color = "black", linewidth = 0.2) +
  geom_vline(xintercept = average_length, color = "red", linetype = "dashed", linewidth = 0.2) +
  labs(
    x = "Read length (bp)", 
    y = "Frequency", 
    title = "Read length (bp) frequency distribution "
  )
p
A1 <- p + theme_classic() +
  theme(
    axis.text.y   = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 9, color = "black"),
    axis.text.x   = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 9, color = "black"),
    axis.title.y  = element_text(size = 9, face = "bold", margin = margin(r = 1)),
    axis.title.x  = element_text(size = 9, face = "bold", margin = margin(t = 1)),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size = 0.01),
    axis.ticks = element_line(colour = "black", size = 0.3),
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9, face = "plain"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.title = element_blank(),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    legend.margin = margin(c(0, 0, 0, 0))
  ) +
  theme(legend.position = "none")
A1
ggsave("Extended_Data_Fig1h.pdf", height=4.5, width=5.25, units='cm', plot = A1)


