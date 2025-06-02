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


##############################################################################################################################################
library(ggplot2)
library("ggplot2")
library("readxl")
library("dplyr")
library(knitr)
library(rmdformats)
library(RColorBrewer)
library(ggpubr) # publication quality figures, based on ggplot2
library(data.table) # alternative to data.frame
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

# ------------------------------------------------------------------------------------------------
# Huila and Antioquia Temperature and pH plots

# Define custom colors
mypal_temp <- c("#1b676b", "#88c425")
mypal_temp2 <- c("grey", "#bef202")
mypal_pH <- c("red", "purple")
mypal_pH2 <- c("red", "purple")

# Function to customize plot theme
custom_theme <- function() {
  theme_classic() + 
    theme(
      axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 7, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 7, color = "black"),
      axis.title.y = element_text(size = 8, face = "bold"),
      axis.title.x = element_text(size = 8, face = "bold"),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 0.35),
      axis.ticks = element_line(colour = "black", size = 0.15),
      legend.title = element_blank(),
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 8, face = "plain"),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8, margin = margin(b = -5)),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
}

# ------------------------------------------------------------------------------------------------
# Temperature plots for Huila and Antioquia

# Huila Temperature
data_set <- subset(tem, Location == 'Huila')
p_Huila <- ggplot(data_set, aes(x = Time, y = Temperature, group = Depth, color = Depth)) +
  geom_vline(xintercept = c(48, 96), linetype = 2, color = 2, size = 0.2) +
  geom_point(aes(fill = Depth), col = 'black', size = 0.5, stroke = 0, position = position_jitterdodge(5)) +
  geom_smooth(alpha = .15, aes(fill = Depth), size = 0.5) + 
  scale_x_continuous(breaks = seq(0, 170, 24)) +
  custom_theme() +
  scale_color_manual(values = mypal_temp) + 
  scale_fill_manual(values = mypal_temp2) +
  ggtitle("Huila") +
  ylab("Temperature (°C)") + 
  xlab("Time (h)")
p_Huila

# Antioquia Temperature
data_set <- subset(tem, Location == 'Antioquia')
p_Antioquia <- ggplot(data_set, aes(x = Time, y = Temperature, group = Depth, color = Depth)) +
  geom_vline(xintercept = c(48, 96), linetype = 2, color = 2, size = 0.2) +
  geom_point(aes(fill = Depth), col = 'black', size = 0.5, stroke = 0, position = position_jitterdodge(5)) +
  geom_smooth(alpha = .15, aes(fill = Depth), size = 0.5) + 
  scale_x_continuous(breaks = seq(0, 170, 24)) +
  custom_theme() +
  scale_color_manual(values = mypal_temp) + 
  scale_fill_manual(values = mypal_temp2) +
  ggtitle("Antioquia") +
  ylab("Temperature (°C)") + 
  xlab("Time (h)")
p_Antioquia

# ------------------------------------------------------------------------------------------------
# pH plots for Huila and Antioquia

# Rename pH labels
B <- c("cotyledon" = "cotyledon", "testa" = "testa/pulp")

# Huila pH
data_set <- subset(pH, Location == 'Huila')
p_Huila_pH <- ggplot(data_set, aes(x = Time, y = pH, group = Fraction, color = Fraction)) +
  geom_vline(xintercept = c(48, 96), linetype = 2, color = 2, size = 0.2) +
  geom_point(aes(fill = Fraction), col = 'black', size = 0.5, stroke = 0, position = position_jitterdodge(5)) +
  geom_smooth(alpha = .05, aes(fill = Fraction), size = 0.5) + 
  scale_x_continuous(breaks = seq(0, 170, 24)) +
  custom_theme() +
  scale_color_manual(values = mypal_pH, labels = B) + 
  scale_fill_manual(values = mypal_pH2, labels = B) +
  ggtitle("Huila") +
  ylab("pH") + 
  xlab("Time (h)")
p_Huila_pH

# Antioquia pH
data_set <- subset(pH, Location == 'Antioquia')
p_Antioquia_pH <- ggplot(data_set, aes(x = Time, y = pH, group = Fraction, color = Fraction)) +
  geom_vline(xintercept = c(48, 96), linetype = 2, color = 2, size = 0.2) +
  geom_point(aes(fill = Fraction), col = 'black', size = 0.5, stroke = 0, position = position_jitterdodge(5)) +
  geom_smooth(alpha = .05, aes(fill = Fraction), size = 0.5) + 
  scale_x_continuous(breaks = seq(0, 170, 24)) +
  custom_theme() +
  scale_color_manual(values = mypal_pH, labels = B) + 
  scale_fill_manual(values = mypal_pH2, labels = B) +
  ggtitle("Antioquia") +
  ylab("pH") + 
  xlab("Time (h)")
p_Antioquia_pH

# ------------------------------------------------------------------------------------------------
# Empty plot for layout alignment
empty_plot <- ggplot() + theme_void()

# ------------------------------------------------------------------------------------------------
# Merge final plots
p_Huila_final <- p_Huila + guides(color = "none", shape = "none", fill = "none")
p_Antioquia_final <- p_Antioquia + guides(color = "none", shape = "none", fill = "none")
p_Huila_pH_final <- p_Huila_pH + guides(color = "none", shape = "none", fill = "none")
p_Antioquia_pH_final <- p_Antioquia_pH + guides(color = "none", shape = "none", fill = "none")

# Arrange plots into a grid
final_plot <- egg::ggarrange(p_Huila_final, empty_plot, p_Antioquia_final, empty_plot, empty_plot, empty_plot, p_Huila_pH_final, empty_plot, p_Antioquia_pH_final,
                             ncol = 3, widths = c(0.9, 0.05, 1), heights = c(1, 0.05, 1))

# Save final plot
ggsave("Fig2b.pdf", height = 6.5, width = 8, units = 'cm', plot = final_plot)

# ------------------------------------------------------------------------------------------------
# Extract and save legends for both Temperature and pH plots
A1_L <- get_legend(p_Huila)
A3_L <- get_legend(p_Huila_pH)

# Convert extracted legends to ggplot
A1_L <- as_ggplot(A1_L)
A3_L <- as_ggplot(A3_L)

# Arrange legends into a single plot
legend_plot <- egg::ggarrange(A1_L, A3_L, ncol = 1, widths = c(1), heights = c(1, 1))
legend_plot

# Save legends
ggsave("Fig2b_legend.pdf", height = 6.5, width = 8, units = 'cm', plot = legend_plot)

