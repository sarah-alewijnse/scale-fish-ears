#### Temp-corrected O2 consumption vs. body mass ####

library(tidyverse)
library(MCMCglmm)
library(here)
library(viridis)

# Load data --------------------------------------------------------------------

data <- read.csv(here("Outputs",
                      "07_C_resp_Models",
                      "Edited_Data",
                      "Scaling_Datasets",
                      "Scaling_Dataset.csv"))

# Load model -------------------------------------------------------------------

o2_mod <- readRDS(here("Outputs",
                       "07_C_resp_Models",
                       "Ch_6_Scaling",
                       "paper_models",
                       "mass_specific.rds"))

# Extract parameters
o2_mod_sum <- summary(o2_mod$Sol)

# Get log o2 consumption
data$ln_O2_mg_kg <- log(data$O2_mg_kg)

# Get model parameters
a <- o2_mod_sum$statistics[1]
b_ln_bm <- o2_mod_sum$statistics[2]
b_inv_temp <- o2_mod_sum$statistics[3]

# Normalise to common temperature ----------------------------------------------

# Fit O2 consumption based on temperature
data$fitted_o2_mg_kg <- a + b_inv_temp * data$Z_inv_temp
plot(data$Z_ln_BM, data$fitted_o2_mg_kg)
plot(data$mean_Temp_K, data$fitted_o2_mg_kg)

# Get residuals
data$resid_o2_ind <- data$ln_O2_mg_kg - data$fitted_o2_mg_kg
plot(data$fitted_o2_mg_kg, data$resid_o2_ind)

plot(data$ln_BM, data$resid_o2_ind)

# Get mean and sd temp
mean_temp <- mean(data$inv_Temp)
sd_temp <- sd(data$inv_Temp)

# Convert 10 C to inverse temp
boltz <- 8.62*10^-5 # Boltzmann constant
standard_inv_Temp <- 1/((273.15+10) * boltz) # inv_temp = 1/kT (T in kelvin)
standard_Z_inv_Temp <- (standard_inv_Temp - mean_temp) / sd_temp # Standardize

# Calculate normalised FMR to 10C
data$FMR10C <- a + b_inv_temp * standard_Z_inv_Temp + data$resid_o2_ind
plot(data$ln_BM, data$FMR10C)
range(data$FMR10C)
range(data$ln_O2_mg_kg)

# Get ci -----------------------------------------------------------------------

# Load function
load(here("Functions",
          "LM_CI.Rdata"))

# Calclate
bm_ci <- lm.ci(data$Z_ln_BM,
                 a = a, b = b_ln_bm,
                 a_low_ci = o2_mod_sum$quantiles[1, 1], 
                 b_low_ci = o2_mod_sum$quantiles[2, 1],
                 a_up_ci = o2_mod_sum$quantiles[1, 5], 
                 b_up_ci = o2_mod_sum$quantiles[2, 5])

# Assign required values to objects --------------------------------------------
# Makes things tidier

min_ln_BM <- min(data$Z_ln_BM)
max_ln_BM <- max(data$Z_ln_BM)
mean_ln_BM <- mean(data$ln_BM)
sd_ln_BM <- sd(data$ln_BM)

# Get species averages ---------------------------------------------------------

# Get means, sd and n

spp_means <- aggregate(data[, c("FMR10C", "Z_ln_BM")], list(data$sciname), mean, na.rm = TRUE)
spp_sds <- aggregate(data[, c("FMR10C", "Z_ln_BM")], list(data$sciname), sd, na.rm = TRUE)
spp_temp <- aggregate(data[, c("mean_Temp")], list(data$sciname), mean, na.rm = TRUE)
spp_n <- aggregate(data[, "n"], list(data$sciname), length)

# Form into nice table 

spp_avgs <- data.frame(sciname = spp_means$Group.1,
                       n = spp_n$x,
                       fmr_mean = spp_means$FMR10C,
                       fmr_sd = spp_sds$FMR10C,
                       mass_mean = spp_means$Z_ln_BM,
                       mass_sd = spp_sds$Z_ln_BM,
                       temp_mean = spp_temp$x)
glimpse(spp_avgs)

# Averages with size = n -------------------------------------------------------

bm_avgs_n <- ggplot(data = spp_avgs, 
                      aes(x = mass_mean, y = fmr_mean)) +
  geom_polygon(data = bm_ci, 
               aes(x_ribbon, y_ribbon), 
               alpha = 0.3) +
  geom_point(aes(size = n, col = temp_mean), 
             alpha = 0.7) +
  geom_segment(aes(x = min_ln_BM,
                   xend = max_ln_BM,
                   y = a + b_ln_bm * min_ln_BM,
                   yend = a + b_ln_bm * max_ln_BM), lwd = 0.5) +
  scale_color_viridis(option = "inferno", 
                      name = "Experienced \ntemperature (°C)",
                      begin = 0.2,
                      end = 0.9) +
  # Customise the theme
  ylab(expression("log oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
  scale_x_continuous(labels = c(0.0, 2.5, 5.0, 7.5, 10.0, 12.5),
                     name = "log body mass (g)",
                     c((0 - mean_ln_BM) / sd_ln_BM,
                       (2.5 - mean_ln_BM) / sd_ln_BM,
                       (5 - mean_ln_BM) / sd_ln_BM,
                       (7.5 - mean_ln_BM) / sd_ln_BM,
                       (10 - mean_ln_BM) / sd_ln_BM,
                       (12.5 - mean_ln_BM) / sd_ln_BM)) +
  scale_y_continuous(breaks = c(3, 4, 5, 6, 7, 8),
                     limits = c(2, 8.2)) +
  theme(panel.background = element_blank(), # Keep the background blank,
        text = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "right",
        legend.title = element_text(),
        legend.key = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
bm_avgs_n

svg(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "mass_specific_comb_spp_avgs_n_bm.svg"), width = 7, height = 5)
bm_avgs_n
dev.off()

# Averages with size = sd -------------------------------------------------------

spp_avgs[is.na(spp_avgs)] <- 0

bm_avgs_sd <- ggplot(data = spp_avgs, 
                       aes(x = mass_mean, y = fmr_mean)) +
  geom_polygon(data = bm_ci, 
               aes(x_ribbon, y_ribbon), 
               alpha = 0.3) +
  geom_point(aes(size = fmr_sd, col = temp_mean), 
             alpha = 0.7) +
  geom_segment(aes(x = min_ln_BM,
                   xend = max_ln_BM,
                   y = a + b_ln_bm * min_ln_BM,
                   yend = a + b_ln_bm * max_ln_BM), lwd = 0.5) +
  scale_color_viridis(option = "inferno", 
                      name = "Experienced \ntemperature (°C)",
                      begin = 0.2,
                      end = 0.9) +
  ## Customise the theme
  ylab(expression("log oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
  scale_x_continuous(labels = c(0.0, 2.5, 5.0, 7.5, 10.0, 12.5),
                     name = "log body mass (g)",
                     c((0 - mean_ln_BM) / sd_ln_BM,
                       (2.5 - mean_ln_BM) / sd_ln_BM,
                       (5 - mean_ln_BM) / sd_ln_BM,
                       (7.5 - mean_ln_BM) / sd_ln_BM,
                       (10 - mean_ln_BM) / sd_ln_BM,
                       (12.5 - mean_ln_BM) / sd_ln_BM)) +
  scale_y_continuous(breaks = c(3, 4, 5, 6, 7, 8),
                     limits = c(2, 8.2)) +
  labs(size = "Standard deviation \nof oxygen consumption") +
  theme(panel.background = element_blank(), # Keep the background blank,
        text = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "right",
        legend.title = element_text(),
        legend.key = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
bm_avgs_sd

svg(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "mass_specific_comb_spp_avgs_sd_bm.svg"), width = 7, height = 5)
bm_avgs_sd
dev.off()

# Error bars -------------------------------------------------------------------

bm_avgs_error <- ggplot(data = spp_avgs, 
                        aes(x = mass_mean, y = fmr_mean)) +
  geom_polygon(data = bm_ci, 
               aes(x_ribbon, y_ribbon), 
               alpha = 0.3) +
  geom_point(aes(col = temp_mean), 
             alpha = 1, size = 3) +
  geom_errorbar(aes(ymin = fmr_mean - fmr_sd,
                    ymax = fmr_mean + fmr_sd,
                    col = temp_mean),
                alpha = 0.4) +
  geom_errorbarh(aes(xmin = mass_mean - mass_sd,
                    xmax = mass_mean + mass_sd,
                    col = temp_mean),
                alpha = 0.4) +
  geom_segment(aes(x = min_ln_BM,
                   xend = max_ln_BM,
                   y = a + b_ln_bm * min_ln_BM,
                   yend = a + b_ln_bm * max_ln_BM), lwd = 0.5) +
  scale_color_viridis(option = "inferno", 
                      name = "Experienced \ntemperature (°C)",
                      begin = 0.2,
                      end = 0.9) +
  ## Customise the theme
  ylab(expression("log"[10]~"oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
  scale_x_continuous(labels = c(0.0, 2.5, 5.0, 7.5, 10.0, 12.5),
                     name = expression("log"[10]~"body mass (g)"),
                     c((0 - mean_ln_BM) / sd_ln_BM,
                       (2.5 - mean_ln_BM) / sd_ln_BM,
                       (5 - mean_ln_BM) / sd_ln_BM,
                       (7.5 - mean_ln_BM) / sd_ln_BM,
                       (10 - mean_ln_BM) / sd_ln_BM,
                       (12.5 - mean_ln_BM) / sd_ln_BM)) +
  scale_y_continuous(breaks = c(3, 4, 5, 6, 7),
                     limits = c(2, 7.2)) +
  labs(size = "Standard deviation \nof oxygen consumption") +
  theme(panel.background = element_blank(), # Keep the background blank,
        text = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "right",
        legend.title = element_text(),
        legend.key = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.text = element_text(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
bm_avgs_error

png(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "mass_specific_comb_spp_avgs_error_bars_bm.png"), width = 7, height = 5,
    res = 250, units = "in")
bm_avgs_error
dev.off()
