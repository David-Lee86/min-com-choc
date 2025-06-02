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
library("ggplot2")
library("readxl")
library("dplyr")
library(knitr)
library(rmdformats)
library(ggpubr) # publication quality figures, based on ggplot2
library(data.table) # alternative to data.frame
library(scales)
library("BiocParallel")

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

# -------------------------------------------------------------------------------------------------
# Import metadata from excel file- Linux
df = read.csv("santander_bean_colour_analysis.csv")

# -------------------------------------------------------------------------------------------------
# Define the row names from the sample_name column
library(tidyverse)
df2 = df %>% remove_rownames %>% column_to_rownames(var="Id")
df2 <- df2 %>% select (-c(Time, Rep))# remove the column variety 

# Perform PCA using prcomp
pca_result <- prcomp(df2, center = TRUE, scale. = TRUE)
pca_result
summary(pca_result)

### scree plot: see how much variance our variables are able to explain
plot(pca_result, type = "l",main="")

# Extract proportion of variance explained
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Create a data frame with PCA results
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df <- tibble::rownames_to_column(pca_df, "Id")
# merge two data frames by ID
total <- merge(df,pca_df,by="Id")

# -------------------------------------------------------------------------------------------------
#adjust sample order
total$Time<- factor(total$Time, levels = c('0hr',"24hr","48hr","72hr","96hr","120hr","144hr","168hr"))

#Rename samples
B=c(
  "0hr"="0h",
  "24hr"="24h",
  "48hr"="48h",
  "72hr"="72h",
  "96hr"="96h",
  "120hr"="120h",
  "144hr"="144h",
  "168hr"="168h"
    )

mypal = c("#437076","#6c9c9e","#92c6c1","#e8e2d2","#f3cec8","#fa98ad","#dd7959",'#8e5e4d')

# plot
A1 = ggplot(total, aes(x = PC1, y = PC2)) +
  geom_point() + 
  labs(x = gsub("\\(\\s*", "(", gsub("\\s*\\)", ")", paste("Axis 1 (", scales::percent(variance_explained[1]), ")"))),
       y = gsub("\\(\\s*", "(", gsub("\\s*\\)", ")", paste("Axis 2 (", scales::percent(variance_explained[2]), ")")))) +
  geom_hline(yintercept = 0, colour="gray70", linetype="dashed") +
  geom_vline(xintercept = 0, colour="gray70", linetype="dashed") +
  geom_point(aes(fill=Time, color=Time), shape = 21, size = 3, colour = "black",stroke=0.3) +
  scale_fill_manual(values = mypal, labels = B) +
  theme_classic() + guides(color=FALSE, shape=FALSE) + guides(fill = guide_legend(title="Sample")) +
  theme(axis.text.y = element_text(angle=0, hjust = 0, vjust = 0.5, size=9, color = "black"),
        axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 1, size=9, color = "black"),
        axis.title.y = element_text(size=9, face="bold"),
        axis.title.x = element_text(size=9, face="bold"),
        legend.title = element_blank(),
        legend.text=element_text(size=9),
        legend.key.size = unit(0.07, "cm"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0.05),
        panel.border = element_rect(colour = "black", fill=NA, size=0.35), plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(legend.position="bottom") +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0))) +
  guides(fill=guide_legend(nrow=2, override.aes = list(size = 5))) +
  scale_y_continuous(expand = c(0.05,0.05)) + scale_x_continuous(expand = c(0.05,0.05))
A1
A1a = A1 + guides(color = "none", shape="none", fill="none")
ggsave("Extended_Data_Fig1f.pdf", height=4.5, width=4.5, units='cm', plot=A1a)



### extract legend option 1
library(ggpubr)
library(egg)
A1L = A1 + guides(fill=guide_legend(ncol=1, override.aes = list(size = 4, stroke=0.6)))
# Extract the legend. Returns a gtable
A1L <- get_legend(A1L)
# Convert extracted legend to a ggplot
A1L = as_ggplot(A1L)
A1L
ggsave("Extended_Data_Fig1f_legend.pdf", height=3, width=3, units='in', plot = A1L)



