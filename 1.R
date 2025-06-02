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
library(ggplot2)
library(readxl)
library(dplyr)
library(knitr)
library(rmdformats)
library(RColorBrewer)
library(ggpubr)
library(data.table)
library(grid)
library(scales)
library(BiocParallel)

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

# -------------------------------------------------------------------------------------------------
# Import data
pH = read.csv("colombia_pH_data.csv")
tem = read.csv("colombia_temperature_data.csv")


# -------------------------------------------------------------------------------------------------
# Filter the Santander data
data_set <- subset(tem, Location == 'Santander')

# Create the temperature plot (A1)
p = ggplot(data_set, aes(x = Time, y = Temperature, group = Depth, color = Depth)) +
  geom_vline(xintercept = c(48, 96), linetype = 2, color = 2, size = 0.2) +
  geom_point(aes(fill = Depth), col = 'black', size = 0.5, stroke = 0, position = position_jitterdodge(5)) +
  geom_smooth(alpha = .15, aes(fill = Depth), size = 0.5) +
  scale_x_continuous(breaks = seq(0, 170, 24)) +
  scale_color_manual(values = c("#1b676b", "#88c425")) +
  scale_fill_manual(values = c("grey", "#bef202")) +
  theme_classic() +
  ylab("Temperature (°C)") + xlab("Time (h)") + ggtitle("Santander") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 9, color = "black"),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_text(size = 9, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.35),
        axis.ticks = element_line(colour = "black", size = 0.35),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.margin = margin(-8, 0, 0, 0),
        legend.spacing.y = unit(0, "cm"))

# Display plot (A1)
A1 = p

# Create the pH plot (A2)
data_set <- subset(pH, Location == 'Santander')
B = c("cotyledon" = "cotyledon", "testa" = "testa/pulp")

p = ggplot(data_set, aes(x = Time, y = pH, group = Fraction, color = Fraction)) +
  geom_vline(xintercept = c(48, 96), linetype = 2, color = 2, size = 0.2) +
  geom_point(aes(fill = Fraction), col = 'black', size = 0.5, stroke = 0, position = position_jitterdodge(5)) +
  geom_smooth(alpha = .05, aes(fill = Fraction), size = 0.5) +
  scale_x_continuous(breaks = seq(0, 170, 24)) +
  scale_color_manual(values = c("Red", "Purple"), labels = B) +
  scale_fill_manual(values = c("Red", "Purple"), labels = B) +
  theme_classic() +
  ylab("pH") + xlab("Time (h)") + ggtitle("Santander") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 9, color = "black"),
        axis.title.y = element_text(size = 9, face = "bold"),
        axis.title.x = element_text(size = 9, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.35),
        axis.ticks = element_line(colour = "black", size = 0.35),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 9),
        panel.border = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.margin = margin(-8, 0, 0, 0),
        legend.spacing.y = unit(0, "cm"))

# Display plot (A2)
A2 = p

##########################################################################################################################
# Create an empty plot (E)
E = ggplot() + theme_void()

# Combine plots A1, E, and A2 into a single layout
Plot2 = egg::ggarrange(A1, E, A2,
                       nrow = 1,
                       widths = c(1, 0.05, 1),
                       heights = c(1))

# Display and save the final plot
Plot2
ggsave("Fig1a.pdf", height = 5.5, width = 11.5, units = 'cm', plot = Plot2)

