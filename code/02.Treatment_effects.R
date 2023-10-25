
# Required packages and my functions
require(glmmTMB)
library(lme4)
require(MASS)
require(emmeans)
require(performance)
require(DHARMa)
require(tidyverse)
source("code/functions.R")

# Read in the data
data <- read.csv("data/data.csv")
glimpse(data)

### (1) The effect of treatment on fungal to bacterial ratio ##################

m1_ratio <- glmmTMB(ratio ~ treatment + (1 | block),
                    data = data,
                    family = gaussian())

m2_ratio <- glmmTMB(ratio ~ treatment + (1 | block/sample_id),
                   data = data,
                   family = gaussian())

# Compare model performance
compare_performance(m1_ratio, m2_ratio)

# Model validation
res_ratio <- simulateResiduals(m1_ratio)
plot(res_ratio)

# Results
summary(m1_ratio)

# Estimated marginal means
emmeans_ratio <- emmeans(m1_ratio, list(pairwise ~ treatment),
                         adjust ="tukey")

### (2) The effect of treatment on fungal abundance ############################

m1_fungal_abundance <- glmer(fungal_abundance/10000 ~ treatment + (1 | block),
                             data = data,
                             family = poisson)

m2_fungal_abundance <- glmer(fungal_abundance/10000 ~ treatment + 
                               (1 | block/sample_id),
                             data = data,
                             family = poisson)

# Compare model performance
compare_performance(m1_fungal_abundance, m2_fungal_abundance)

# Model validation
res_fungal_abundance <- simulateResiduals(m2_fungal_abundance)
plot(res_fungal_abundance)

# Results
summary(m2_fungal_abundance)

# Estimated marginal means
emmeans_fungal_abundance <- emmeans(m2_fungal_abundance,
                                    list(pairwise ~ treatment),
                                    adjust ="tukey")
summary(m2_fungal_abundance)
### (3) The effect of treatment on bacterial abundance #########################

m1_bacterial_abundance <- glmer(bacterial_abundance/10000 ~ treatment + (1 | block),
                                data = data,
                                family = poisson)

m2_bacterial_abundance <- glmer(bacterial_abundance/10000 ~ treatment +
                                  (1 | block/sample_id),
                                data = data,
                                family = poisson)

# Compare model performance
compare_performance(m1_bacterial_abundance, m2_bacterial_abundance)

# Model validation
res_bacterial_abundance <- simulateResiduals(m2_bacterial_abundance)
plot(res_bacterial_abundance)

# Results
summary(m2_bacterial_abundance)

# Estimated marginal means
emmeans_bacterial_abundance <- emmeans(m2_bacterial_abundance,
                                       list(pairwise ~ treatment),
                                       adjust ="tukey")

### (4) The effect of treatment on carbon ######################################

m1_carbon <- glmmTMB(carbon ~ treatment + (1 | block),
                     data = data,
                     family = gaussian())

m2_carbon <- glmmTMB(carbon ~ treatment + (1 | block/sample_id),
                     data = data,
                     family = gaussian())

# Compare model performance
compare_performance(m1_carbon, m2_carbon)

# Model validation
res_carbon <- simulateResiduals(m1_carbon)
plot(res_carbon)

# Results
summary(m1_carbon)

# Estimated marginal means
emmeans_carbon <- emmeans(m1_carbon, list(pairwise ~ treatment),
                          adjust ="tukey")

# Pairwise differences
pairwise_diff_carbon <- pairs(emmeans_carbon) %>%
  as_tibble()
write_csv(pairwise_diff_carbon, "output/pairwise_diff_carbon.csv")

# Estimated marginal means
emmeans_data_carbon <- emmeans_carbon %>%
  as.data.frame() %>%
  slice(1:6) %>%
  select(treatment, emmean, SE, df, lower.CL, upper.CL)

### (5) The effect of treatment on nitrogen ####################################

m1_nitrogen <- glmmTMB(nitrogen ~ treatment + (1 | block),
                       data = data,
                       family = gaussian())

m2_nitrogen <- glmmTMB(nitrogen ~ treatment + (1 | block/sample_id),
                       data = data,
                       family = gaussian())

# Compare model performance
compare_performance(m1_nitrogen, m2_nitrogen)

# Model validation
res_nitrogen <- simulateResiduals(m2_nitrogen)
plot(res_nitrogen)

# Lets see what happens when we remove that outlier we detected in the data 
# exploration
outliers(res_nitrogen)
data_no_nitrogen_outlier <- data %>%
  slice(-1)

m3_nitrogen <- glmmTMB(nitrogen ~ treatment + (1 | block),
                       data = data_no_nitrogen_outlier,
                       family = gaussian())

m4_nitrogen <- glmmTMB(nitrogen ~ treatment + (1 | block/sample_id),
                       data = data_no_nitrogen_outlier,
                       family = gaussian())

# Compare model performance
compare_performance(m3_nitrogen, m4_nitrogen)

# Model validation
res_nitrogen <- simulateResiduals(m4_nitrogen)
plot(res_nitrogen)
# The distribution is still off
# I'll see what happens with a log transformation to squeeze down the upper 
# extreme values that were detected in the data exploration

m5_nitrogen <- glmmTMB(log(nitrogen) ~ treatment + (1 | block),
                       data = data,
                       family = gaussian())

m6_nitrogen <- glmmTMB(log(nitrogen) ~ treatment + (1 | block/sample_id),
                       data = data,
                       family = gaussian())

m7_nitrogen <- glmmTMB(log(nitrogen) ~ treatment + (1 | block),
                       data = data_no_nitrogen_outlier,
                       family = gaussian())

m8_nitrogen <- glmmTMB(log(nitrogen) ~ treatment + (1 | block/sample_id),
                       data = data_no_nitrogen_outlier,
                       family = gaussian())

# Compare model performance
compare_performance(m5_nitrogen, m6_nitrogen, m7_nitrogen, m8_nitrogen)

# Model validation
res_nitrogen <- simulateResiduals(m8_nitrogen)
plot(res_nitrogen)
# Now the homogeneity of variances is off. I'll try the model with the outlier.
res_nitrogen <- simulateResiduals(m5_nitrogen)
plot(res_nitrogen)
# That looks like the best fit

# Results
summary(m5_nitrogen)

# Estimated marginal means
emmeans_nitrogen <- emmeans(m5_nitrogen, list(pairwise ~ treatment),
                            adjust ="tukey")

### (6) The effect of treatment on carbon to nitrogen ratio ####################

m1_CN <- glmmTMB(C.N ~ treatment + (1 | block),
                 data = data,
                 family = gaussian())

m2_CN <- glmmTMB(C.N ~ treatment + (1 | block/sample_id),
                 data = data,
                 family = gaussian())

# Compare model performance
compare_performance(m1_CN, m2_CN)

# Model validation
res_CN <- simulateResiduals(m1_CN)
plot(res_CN)

# Results
summary(m1_CN)

# Estimated marginal means
emmeans_CN <- emmeans(m1_CN, list(pairwise ~ treatment),
                      adjust ="tukey")

# Pairwise differences
pairwise_diff_CN <- pairs(emmeans_CN) %>%
  as_tibble()


# Estimated marginal means
emmeans_data_CN <- emmeans_CN %>%
  as.data.frame() %>%
  slice(1:6) %>%
  select(treatment, emmean, SE, df, lower.CL, upper.CL)

### (7) The effect of treatment on pH ##########################################

m1_pH <- glmmTMB(pH ~ treatment + (1 | block),
                 data = data,
                 family = gaussian())

m2_pH <- glmmTMB(pH ~ treatment + (1 | block/sample_id),
                 data = data,
                 family = gaussian())

# Compare model performance
compare_performance(m1_pH, m2_pH)

# Model validation
res_pH <- simulateResiduals(m2_pH)
plot(res_pH)

# Results
summary(m2_pH)

# Estimated marginal means
emmeans_pH <- emmeans(m2_pH, list(pairwise ~ treatment),
                      adjust ="tukey")

# Pairwise differences
pairwise_diff_pH <- pairs(emmeans_pH) %>%
  as_tibble()

# Estimated marginal means
emmeans_data_pH <- emmeans_pH %>%
  as.data.frame() %>%
  slice(1:6) %>%
  select(treatment, emmean, SE, df, lower.CL, upper.CL)

### (8) Effect plots ###########################################################

# Carbon
jitter_carbon <- data %>%
  select(treatment, emmean = carbon) %>%
  mutate(lower.CL = emmean) %>%
  mutate(upper.CL = emmean) %>%
  mutate(SE = emmean)

ggplot_carbon <- ggplot(data = emmeans_data_carbon,
                        aes(x = treatment, 
                            ymin = lower.CL,
                            lower = emmean - SE,
                            middle = emmean,
                            upper = emmean + SE,
                            ymax = upper.CL,
                            fill = treatment))  +
  geom_boxplot(stat = 'identity') +
  geom_point(jitter_carbon,
             mapping = aes(x = treatment, y = emmean, fill = treatment),
             position = position_jitter(seed = 16, width = 0.25),
             shape = 21,
             size = 2,
             alpha = 0.66
  ) +
  scale_fill_manual(
    values = c("#D55E00", "#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#F0E442"),
    limits = c('E1', 'E2', 'E3', 'E5', 'L2', 'U')) +
  theme(legend.position= "none") +
  scale_y_continuous(limits = c(0, 3)) +
  xlab("Treatment") +
  ylab("Percent soil organic carbon") +
  MyTheme()

ggplot_carbon2 <- ggplot(data = emmeans_data_carbon,
                         aes(x = treatment, y = emmean)) + 
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.4,
    linewidth = 0.5) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0,
    linewidth = 0.7,
    linetype = 'dotted') +
  geom_point(size = 5, shape = 23, fill = "#36454f") +
  scale_y_continuous(limits = c(0, 2.2),
                     breaks = c(0, 1, 2)) +
  xlab("Treatment") +
  ylab("Soil total carbon (mg/kg)") +
  MyTheme()
