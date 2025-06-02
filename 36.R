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
library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(egg)
library(ggpubr)

# -------------------------------------------------------------------------------------------------
# Set working directory automatically (RStudio only)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
cat("Current working directory:", getwd(), "\n")

# -------------------------------------------------------------------------------------------------
# Set a theme for all plots
theme_set(theme_classic() +
            theme(axis.text.y = element_text(angle=0, hjust = 0, vjust = 0.5, size=7, color = "black"),
               axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 1, size=7, color = "black"),
               axis.title.y = element_text(size=7, face="bold"),
               axis.title.x = element_text(size=7, face="bold"),
               legend.title = element_blank(),
               legend.text=element_text(size=10),
               legend.key.size = unit(0.3, "cm"),
               panel.background = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black", size = 0),
               panel.border = element_rect(colour = "black", fill=NA, size=0.35),
               plot.title = element_text(hjust = 0.5, face="bold", size=0)) +
  theme(legend.position="bottom") +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")))


# Define colors
mypal <- c('cot_pH.rescaled'="#4c2c51",
           'Cot_pH'="#4c2c51",
           'Testa_pH'="#c1204c",
           'testa_pH.rescaled'="#c1204c",
           'Temp_mid_box'='#e5e1e0',
           'Temp_top'="#7dabcf",
           'tem.rescaled'="#a1af01"
)

B <- c('cot_pH.rescaled'="Cotyledon pH",
       'Cot_pH'="Cotyledon pH",
       'Testa_pH'="Testae/pulp pH ",
       'testa_pH.rescaled'="Testae/pulp pH ",
       'Temp_mid_box'='Bean temp. (mid-box)',
       'Temp_top'="Bean temp. at 7cm depth ",
       'tem.rescaled'="Bean temperature"
)

# -------------------------------------------------------------------------------------------------
# pH Distribution Analysis

# Import and prepare data
tem.pH.data <- read.csv("santander_temperature_pH_distribution.csv") %>%
  select(-Temp_mid_box, -Temp_top)

# Kolmogorov-Smirnov test
ks_result_pH <- ks.test(tem.pH.data$Testa_pH, tem.pH.data$Cot_pH)
p_value_sci <- format(ks_result_pH$p.value, scientific = TRUE, digits = 4)
p_value_sci <- gsub("e", "e", p_value_sci)
p_value_sci

# Prepare data for plotting
data_melt <- melt(tem.pH.data)

# Create plot
pH_plot <- ggplot(data_melt, aes(x=value, fill=variable)) +
  geom_density(alpha=0.5, size=0.25) +
  scale_fill_manual(values = mypal, 
                    labels = B) +
  labs(x = "pH", y = "Density") +
  xlim(3, 7) +
  annotate("text", size=2, 
           label = paste0("italic(p)==", '"', p_value_sci, '"'),  
           x = 5.7, y = 1, parse = TRUE)
pH_plot

# -------------------------------------------------------------------------------------------------
# Temperature Distribution Analysis

# Import and prepare data
tem.pH.data <- read.csv("santander_temperature_pH_distribution.csv") %>%
  select(-Testa_pH, -Cot_pH)

# Kolmogorov-Smirnov test
ks_result_temp <- ks.test(tem.pH.data$Temp_mid_box, tem.pH.data$Temp_top)
p_value_sci <- format(ks_result_temp$p.value, scientific = TRUE, digits = 4)
p_value_sci <- gsub("e", "e", p_value_sci)
p_value_sci

# Prepare data for plotting
data_melt <- melt(tem.pH.data)

# Create plot
temp_plot <- ggplot(data_melt, aes(x=value, fill=variable)) +
  geom_density(alpha=0.5, size=0.25) +
  scale_fill_manual(values = mypal, 
                    labels = B) +
  labs(x = "Temperature", y = "Density") +
  xlim(20, 55) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  annotate("text", size=2, 
           label = paste0("italic(p)==", '"', p_value_sci, '"'), 
           x = 40, y = 0.04, parse = TRUE)
temp_plot

# -------------------------------------------------------------------------------------------------
# Temperature vs Cotyledon pH Analysis

# Import and prepare data
tem.pH.data <- read.csv("santander_temperature_pH_distribution.csv") %>%
  mutate(Temp_mean = rowMeans(select(., Temp_mid_box, Temp_top), na.rm = TRUE)) %>%
  select(-Testa_pH, -Temp_mid_box, -Temp_top) %>%
  mutate(tem.rescaled = rescale(Temp_mean, to = c(0, 1)),
         cot_pH.rescaled = rescale(Cot_pH, to = c(0, 1))) %>%
  select(-Cot_pH, -Temp_mean)

# Kolmogorov-Smirnov test
ks_result_temp_cot <- ks.test(tem.pH.data$tem.rescaled, tem.pH.data$cot_pH.rescaled)
p_value_sci <- format(ks_result_temp_cot$p.value, scientific = TRUE, digits = 4)
p_value_sci <- gsub("e", "e", p_value_sci)
p_value_sci

# Prepare data for plotting
data_melt <- melt(tem.pH.data)

# Create plot
temp_cot_plot <- ggplot(data_melt, aes(x=value, fill=variable)) +
  geom_density(alpha=0.5, size=0.25) +
  scale_fill_manual(values = mypal, 
                    labels = B) +
  labs(x = "Temperature and pH", y = "Density") +
  xlim(0,1) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  annotate("text", size=2, 
           label = paste0("italic(p)==", '"', p_value_sci, '"'),  
           x = 0.3, y = 1.3, parse = TRUE)
temp_cot_plot

# -------------------------------------------------------------------------------------------------
# Temperature vs Testa pH Analysis

# Import and prepare data
tem.pH.data <- read.csv("santander_temperature_pH_distribution.csv") %>%
  mutate(Temp_mean = rowMeans(select(., Temp_mid_box, Temp_top), na.rm = TRUE)) %>%
  select(-Cot_pH, -Temp_mid_box, -Temp_top) %>%
  mutate(tem.rescaled = rescale(Temp_mean, to = c(0, 1)),
         testa_pH.rescaled = rescale(Testa_pH, to = c(0, 1))) %>% 
  select(-Testa_pH, -Temp_mean)
         
# Kolmogorov-Smirnov test
ks_result_temp_testa <- ks.test(tem.pH.data$tem.rescaled, tem.pH.data$testa_pH.rescaled)
p_value_sci <- format(ks_result_temp_testa$p.value, scientific = TRUE, digits = 4)
p_value_sci <- gsub("e", "e", p_value_sci)
p_value_sci

# Prepare data for plotting
data_melt <- melt(tem.pH.data)

# Create plot
temp_testa_plot <- ggplot(data_melt, aes(x=value, fill=variable)) +
  geom_density(alpha=0.5, size=0.25) +
  scale_fill_manual(values = mypal, 
  labels = B) +
  labs(x = "Temperature and pH", y = "Density") +
  xlim(0,1) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  annotate("text", size=2, 
           label = paste0("italic(p)==", '"', p_value_sci, '"'), 
  x = 0.3, y = 0.3, parse = TRUE)
temp_testa_plot


# -------------------------------------------------------------------------------------------------
# Combine and Save Plots

# Create combined plot without legends
combined_plot <- egg::ggarrange(
  temp_testa_plot + guides(fill = "none"),
  temp_cot_plot + guides(fill = "none"),
  temp_plot + guides(fill = "none"),
  pH_plot + guides(fill = "none"),
  ncol = 2
)
combined_plot
# Save combined plot
ggsave("Extended_Data_Fig1c.pdf", plot = combined_plot, height = 6.5, width = 7, units = "cm")


### extract legend option 1
library(ggpubr)
library(egg)
A1L = temp_testa_plot + guides(fill=guide_legend(nrow=2, override.aes = list(size = 4, stroke=0.6)))
A2L = temp_cot_plot + guides(fill=guide_legend(nrow=2, override.aes = list(size = 4, stroke=0.6)))
A3L = temp_plot + guides(fill=guide_legend(nrow=2, override.aes = list(size = 4, stroke=0.6)))
A4L = pH_plot + guides(fill=guide_legend(nrow=2, override.aes = list(size = 4, stroke=0.6)))
# Extract the legend. Returns a gtable
A1L <- get_legend(A1L)
A2L <- get_legend(A2L)
A3L <- get_legend(A3L)
A4L <- get_legend(A4L)
# Convert extracted legend to a ggplot
A1L = as_ggplot(A1L)
A2L = as_ggplot(A2L)
A3L = as_ggplot(A3L)
A4L = as_ggplot(A4L)
##plot
legend = egg::ggarrange(A1L,
                        A2L,
                        A3L,
                        A4L,
                        #ncol = 4,
                        nrow = 4,
                        widths = c(1), heights=c(1,1,1,1))
legend
ggsave("Extended_Data_Fig1c_legend.pdf", height=3, width=3, units='in', plot = legend)

