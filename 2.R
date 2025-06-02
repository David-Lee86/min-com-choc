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
library(dplyr)
library(grid)
library(ggpubr)
library(egg)
library(BiocParallel)

# -------------------------------------------------------------------------------------------------
# Automatically set working directory to script location (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# -------------------------------------------------------------------------------------------------
# Import data
dataset = read.csv("santander_temperature_and_pH_data.csv")

# -------------------------------------------------------------------------------------------------
p = ggplot(dataset, aes(x=testa_pH, y=Temperature)) + 
  geom_point(shape = 21, fill = '#5acb2f', size = 2, alpha = 0.6, stroke = 0.3, color = "black") +
  geom_smooth(method=lm, size=0.6, fill = "#4cbefa")
p
A2 = p + theme_classic() + ylab("Temperature (°C)") + xlab("Testa/pulp pH") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=9, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=9, color = "black"),
        axis.title.y  = element_text(size=9, face="bold"),
        axis.title.x  = element_text(size=9, face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.35),
        axis.ticks = element_line(colour = "black", size = 0.35),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=9, face=c("plain")),
        panel.border = element_rect(colour = "white", fill=NA, size=0.4),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A2
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$testa_pH, dataset$Temperature, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.94, hjust = 0, gp = gpar(col = "black", fontsize = 9, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.8, hjust = 0, gp = gpar(col = "black", fontsize = 9, fontface = "bold")))
A2 <- A2 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A2
ggsave("Fig1b.pdf", plot = A2, height=5.5, width=5.7, units='cm')


# -------------------------------------------------------------------------------------------------
p = ggplot(dataset, aes(x=testa_pH, y=Cotyledon_pH)) + 
  geom_point(shape = 21, fill = '#5acb2f', size = 3, alpha = 0.6, stroke = 0.3, color = "black") +
  geom_smooth(method=lm, size=0.6, fill = "#4cbefa")
p
A1 = p + theme_classic() + ylab("Cotyledon pH") + xlab("Testa/pulp pH") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=9, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=9, color = "black"),
        axis.title.y  = element_text(size=9, face="bold"),
        axis.title.x  = element_text(size=9, face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.35),
        axis.ticks = element_line(colour = "black", size = 0.15),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=9, face=c("plain")),
        panel.border = element_blank(),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A1
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$testa_pH, dataset$Cotyledon_pH, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.3, hjust = 0, gp = gpar(col = "black", fontsize = 9, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.1, hjust = 0, gp = gpar(col = "black", fontsize = 9, fontface = "bold")))
A1 <- A1 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A1


p = ggplot(dataset, aes(x=Cotyledon_pH, y=Temperature)) + 
  geom_point(shape = 21, fill = '#5acb2f', size = 3, alpha = 0.6, stroke = 0.3, color = "black") +
  geom_smooth(method=lm, size=0.6, fill = "#4cbefa")
p
A3 = p + theme_classic() + ylab("Temperature (°C)") + xlab("Cotyledon pH") + ggtitle("Santander") + theme(legend.position="bottom") +
  theme(axis.text.y   = element_text(angle=0, hjust = 0.5, vjust = 0.5, size=9, color = "black"),
        axis.text.x   = element_text(angle=0, hjust = 0.5, vjust = 1, size=9, color = "black"),
        axis.title.y  = element_text(size=9, face="bold"),
        axis.title.x  = element_text(size=9, face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size=0.35),
        axis.ticks = element_line(colour = "black", size = 0.15),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text=element_text(size=9, face=c("plain")),
        panel.border = element_blank(),
        plot.title = element_blank()) +
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm")) +
  theme(legend.margin=margin(c(0,0,0,0)))
A3
# Calculate Pearson correlation and p-value manually
cor_test_result <- cor.test(dataset$Cotyledon_pH, dataset$Temperature, method = "pearson")
r_value <- cor_test_result$estimate
p_value <- cor_test_result$p.value
r_value_text <- as.character(round(r_value, 5))
formatted_p_value <- format(p_value, scientific = TRUE, digits = 5)
R_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(r) == bold(.(r_value_text)))),
                                  x = 0.05, y = 0.3, hjust = 0, gp = gpar(col = "black", fontsize = 9, fontface = "bold")))
P_value_grob <- grobTree(textGrob(label = as.expression(bquote(bolditalic(p) == bold(.(formatted_p_value)))),
                                  x = 0.05, y = 0.1, hjust = 0, gp = gpar(col = "black", fontsize = 9, fontface = "bold")))
A3 <- A3 + annotation_custom(R_value_grob) + annotation_custom(P_value_grob)
A3


#plot using egg
Plot = egg::ggarrange(A3,A1, ncol = 2, widths = c(1,1), heights = c(1))
ggsave("Extended_Data_Fig1d.pdf", plot = Plot, height=2, width=4, units='in')

