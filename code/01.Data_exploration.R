
# Required packages and my functions
require(tidyverse)
source("code/functions.R")

# Read in the data
data <- read.csv("data/data.csv")
glimpse(data)

### Step 1: Investigate outliers ###############################################

# Predictor variables outliers (extreme values)
data %>%
  # Remove unwanted or select desired variables
  select(-c(plot, block, treatment)) %>%
  pivot_longer(-sample_id) %>%
  # Define plot order: arrange by ascending values and group by variable
  arrange(value) %>%
  arrange(name) %>%
  ggplot(aes(x = 1:nrow(.), y = value)) +
  geom_point() +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  labs(x = 'Sample ordered by incresing value', y = 'Value of the variable') +
  theme(axis.text.x=element_blank()) +
  MyTheme()
ggsave('output/outliers_extreme_values.jpg', width = 17, unit = 'cm')

# Correlation outliers: Bacterial abundance
data %>%
  # Remove unwanted or select desired variables
  select(-c(sample_id, plot, block, treatment, fungal_abundance, total_abundance)) %>%
  pivot_longer(-bacterial_abundance) %>%
  ggplot(aes(value, bacterial_abundance, group = name)) +
  geom_point() +
  geom_smooth(method = 'loess', level = 0, colour = 'red', linewidth = 0.5) +
  geom_smooth(method = 'glm', level = 0.95, linewidth = 1) +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  labs(x = 'Predictor variable', y = 'Response variable') +
  MyTheme()

# Correlation outliers: Bacterial abundance
data %>%
  # Remove unwanted or select desired variables
  select(-c(sample_id, plot, block, treatment)) %>%
  pivot_longer(-fungal_abundance) %>%
  ggplot(aes(value, fungal_abundance, group = name)) +
  geom_point() +
  geom_smooth(method = 'loess', level = 0, colour = 'red', linewidth = 0.5) +
  geom_smooth(method = 'glm', level = 0.95, linewidth = 1) +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  labs(x = 'Predictor variable', y = 'Response variable') +
  MyTheme()

# Correlation outliers: Total abundance
data %>%
  # Remove unwanted or select desired variables
  select(-c(sample_id, plot, block, treatment)) %>%
  pivot_longer(-total_abundance) %>%
  ggplot(aes(value, total_abundance, group = name)) +
  geom_point() +
  geom_smooth(method = 'loess', level = 0, colour = 'red', linewidth = 0.5) +
  geom_smooth(method = 'glm', level = 0.95, linewidth = 1) +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  labs(x = 'Predictor variable', y = 'Response variable') +
  MyTheme()
ggsave('output/outliers_correlation.jpg', width = 15, unit = 'cm')

# I need to consider removing the nitrogen as well as the fungal abundance
# outliers

### Step 2: Investigate collinearity ###########################################

# Variation inflation factor

data %>%
  select(-c(sample_id, plot, block, treatment, fungal_abundance, 
            bacterial_abundance, ratio, total_abundance)) %>%
  VIF()
# As expected, carbon, nitrogen and C:N ratio co-vary. Despite nitrogen having
# the highest VIF, I will remove C:N as I'd rather conserve nitrogen and carbon
# considering they had notable inverse effects in the correlation plots, whereas
# C:N ratio had no effect (i.e., carbon and nitrogen effects effectively,
# neutralised each other)

data %>%
  select(-c(sample_id, plot, block, treatment, fungal_abundance, 
            bacterial_abundance, ratio, total_abundance, C.N)) %>%
  VIF()

# Correlation matrix
data %>%
  # Remove unwanted or select desired variables
  select(-c(sample_id, plot, block, total_abundance)) %>%
  GGally::ggpairs()

# PCA
pca <- data %>%
  # Remove unwanted or select desired variables
  select(-c(sample_id, plot, block, treatment, total_abundance)) %>%
  # Compute PCA
  prcomp(., scale. = TRUE)
# Plot PCA
AMR::ggplot_pca(pca) +
  MyTheme() +
  scale_y_continuous(limits = c(-2.5, 3)) +
  scale_x_continuous(limits = c(-3.5, 2)) + coord_fixed()
ggsave('output/correlation_PCA.jpg', width = 15, height = 15, unit = 'cm')
