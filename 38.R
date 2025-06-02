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
library(ggplot2)
library(dplyr)
library(grid)
library(ggpubr)
library(egg)
library(BiocParallel)

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# Register parallel processing with 6 cores
register(SnowParam(6))

# -------------------------------------------------------------------------------------------------
# Import data
dataset = read.csv("santander_pH_tem_colour_data.csv")

# ----------------------------------------- Temperature --------------------------------------------------------
p = ggplot(dataset, aes(x=Grayscale, y=Temperature)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A1 = p + theme_classic() + ylab("Temperature (°C)") + xlab("Grayscale") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A1
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Grayscale, dataset$Temperature, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
# Convert r_value to character to ensure it's treated as text
r_value_text <- as.character(round(r_value, 5))
# Format p-value in scientific notation
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
# Create grob for R-value with bold text
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
# Create grob for p-value with scientific notation and bold text
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
# Add the annotations to the plot
A1 <- A1 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A1


p = ggplot(dataset, aes(x=Luminance, y=Temperature)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A2 = p + theme_classic() + ylab("Temperature (°C)") + xlab("Luminance") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A2
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Luminance, dataset$Temperature, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A2 <- A2 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A2


p = ggplot(dataset, aes(x=Blue, y=Temperature)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A3 = p + theme_classic() + ylab("Temperature (°C)") + xlab("Blue") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A3
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Blue, dataset$Temperature, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A3 <- A3 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A3


p = ggplot(dataset, aes(x=Green, y=Temperature)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A4 = p + theme_classic() + ylab("Temperature (°C)") + xlab("Green") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A4
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Green, dataset$Temperature, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A4 <- A4 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A4


p = ggplot(dataset, aes(x=Red, y=Temperature)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A5 = p + theme_classic() + ylab("Temperature (°C)") + xlab("Red") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A5
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Red, dataset$Temperature, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A5 <- A5 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A5


figure = egg::ggarrange(A3,A4,A5,A1,A2, ncol = 5, widths = c(1,1,1,1,1), heights = c(1))
F1= annotate_figure(figure,
                    left = text_grob("Temperature (°C)", color = "black", rot = 90, size = 4, face="bold")
)


# -------------------------------------- Cotyledon_pH --------------------------------------------------------
p = ggplot(dataset, aes(x=Grayscale, y=Cotyledon_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A1 = p + theme_classic() + ylab("Cotyledon pH") + xlab("Grayscale") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A1
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Grayscale, dataset$Cotyledon_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A1 <- A1 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A1


p = ggplot(dataset, aes(x=Luminance, y=Cotyledon_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A2 = p + theme_classic() + ylab("Cotyledon pH") + xlab("Luminance") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A2
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Luminance, dataset$Cotyledon_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A2 <- A2 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A2


p = ggplot(dataset, aes(x=Blue, y=Cotyledon_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A3 = p + theme_classic() + ylab("Cotyledon pH") + xlab("Blue") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A3
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Blue, dataset$Cotyledon_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A3 <- A3 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A3


p = ggplot(dataset, aes(x=Green, y=Cotyledon_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A4 = p + theme_classic() + ylab("Cotyledon pH") + xlab("Green") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A4
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Green, dataset$Cotyledon_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A4 <- A4 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A4


p = ggplot(dataset, aes(x=Red, y=Cotyledon_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A5 = p + theme_classic() + ylab("Cotyledon pH") + xlab("Red") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A5
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Red, dataset$Cotyledon_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A5 <- A5 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A5


figure = egg::ggarrange(A3,A4,A5,A1,A2, ncol = 5, widths = c(1,1,1,1,1), heights = c(1))
F2= annotate_figure(figure,
                    left = text_grob("Cotyledon pH", color = "black", rot = 90, size = 4, face="bold"),
)


# -------------------------------------- Testa/plup_pH --------------------------------------------------------
p = ggplot(dataset, aes(x=Grayscale, y=Testae_plup_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A1 = p + theme_classic() + ylab("Testa/pulp pH") + xlab("Grayscale") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A1
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Grayscale, dataset$Testae_plup_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A1 <- A1 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A1


p = ggplot(dataset, aes(x=Luminance, y=Testae_plup_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A2 = p + theme_classic() + ylab("Testa/pulp pH") + xlab("Luminance") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A2
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Luminance, dataset$Testae_plup_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A2 <- A2 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A2


p = ggplot(dataset, aes(x=Blue, y=Testae_plup_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A3 = p + theme_classic() + ylab("Testa/pulp pH") + xlab("Blue") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A3
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Blue, dataset$Testae_plup_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A3 <- A3 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A3


p = ggplot(dataset, aes(x=Green, y=Testae_plup_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A4 = p + theme_classic() + ylab("Testa/pulp pH") + xlab("Green") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A4
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Green, dataset$Testae_plup_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A4 <- A4 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A4


p = ggplot(dataset, aes(x=Red, y=Testae_plup_pH)) + 
  geom_point(shape = 3, size = 3, alpha = 0.3, stroke = 0.1, color = "black") +
  geom_smooth(method=lm, colour="red", fill = "grey", size=0.4)
A5 = p + theme_classic() + ylab("Testa/pulp pH") + xlab("Red") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=4.2, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=4.2, color = "black"),
        axis.title.y  = element_blank(),
        axis.title.x  = element_text(size=4.2, face="bold", margin = margin(b = 0.1)),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.15),
        axis.ticks = element_line(colour = "black", size = 0.1), axis.ticks.length = unit(0.1, "cm"),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=4.2, face=c("plain")),
        panel.border = element_rect(colour = "black", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.05,0.01,0.01,0.01),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A5
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Red, dataset$Testae_plup_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.9, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.7, hjust = 0, gp = gpar(col = "black", fontsize = 6, fontface = "bold")))
A5 <- A5 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A5


figure = egg::ggarrange(A3,A4,A5,A1,A2, ncol = 5, widths = c(1,1,1,1,1), heights = c(1))
F3= annotate_figure(figure,
                    left = text_grob("Testa/plup pH", color = "black", rot = 90, size = 4, face="bold"),
)


# -------------------------------------- Testa/plup_pH --------------------------------------------------------
## combine all three
#plot using egg
Plot = egg::ggarrange(F1,F2,F3, nrow = 3)
ggsave("Extended_Data_Fig1g.pdf", plot = Plot, height=6.5, width=10.5, units='cm')



