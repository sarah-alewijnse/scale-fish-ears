#### Temp-corrected O2 consumption vs. body mass ####

library(tidyverse)
library(MCMCglmm)
library(here)
library(viridis)

# Load data --------------------------------------------------------------------

source(here("code",
            "00_pre_processing",
            "03_oxygen_consumption",
            "02_scaling_dataset_main.R"))
# Ignore duplicated data warning

# Load model -------------------------------------------------------------------

o2_mod <- readRDS(here("outputs",
                       "models",
                       "main",
                       "mass_specific.rds"))

# Extract parameters
o2_mod_sum <- summary(o2_mod$Sol)

# Get model parameters
a <- o2_mod_sum$statistics[1]
b_ln_bm <- o2_mod_sum$statistics[2]
b_inv_temp <- o2_mod_sum$statistics[3]

# Normalise to common temperature ----------------------------------------------

# Fit O2 consumption based on temperature
new_fish_data$fitted_o2_mg_kg <- a + b_inv_temp * new_fish_data$Z_inv_temp
plot(new_fish_data$Z_ln_body_mass, new_fish_data$fitted_o2_mg_kg)
plot(new_fish_data$temp_kelvin_mean, new_fish_data$fitted_o2_mg_kg)

# Get residuals
new_fish_data$resid_o2_ind <- new_fish_data$log_o2_mg_kg - new_fish_data$fitted_o2_mg_kg
plot(new_fish_data$fitted_o2_mg_kg, new_fish_data$resid_o2_ind)

plot(new_fish_data$ln_body_mass, new_fish_data$resid_o2_ind)

# Get mean and sd temp
mean_temp <- mean(new_fish_data$inv_temp_mean)
sd_temp <- sd(new_fish_data$inv_temp_mean)

# Convert 10 C to inverse temp
boltz <- 8.62*10^-5 # Boltzmann constant
standard_inv_Temp <- 1/((273.15+10) * boltz) # inv_temp = 1/kT (T in kelvin)
standard_Z_inv_Temp <- (standard_inv_Temp - mean_temp) / sd_temp # Standardize

# Calculate normalised FMR to 10C
new_fish_data$FMR10C <- a + b_inv_temp * standard_Z_inv_Temp + new_fish_data$resid_o2_ind
plot(new_fish_data$ln_body_mass, new_fish_data$FMR10C)
range(new_fish_data$FMR10C)
range(new_fish_data$ln_o2_mg_kg)

# Get ci -----------------------------------------------------------------------

# Load function
load(here("functions",
          "LM_CI.Rdata"))

# Calclate
bm_ci <- lm.ci(new_fish_data$Z_ln_body_mass,
                 a = a, b = b_ln_bm,
                 a_low_ci = o2_mod_sum$quantiles[1, 1], 
                 b_low_ci = o2_mod_sum$quantiles[2, 1],
                 a_up_ci = o2_mod_sum$quantiles[1, 5], 
                 b_up_ci = o2_mod_sum$quantiles[2, 5])

# Assign required values to objects --------------------------------------------
# Makes things tidier

min_ln_BM <- min(new_fish_data$Z_ln_body_mass)
max_ln_BM <- max(new_fish_data$Z_ln_body_mass)
mean_ln_BM <- mean(new_fish_data$ln_body_mass)
sd_ln_BM <- sd(new_fish_data$ln_body_mass)

# Get species averages ---------------------------------------------------------

# Get means, sd and n

spp_means <- aggregate(new_fish_data[, c("FMR10C", "Z_ln_body_mass")], list(new_fish_data$species), mean, na.rm = TRUE)
spp_sds <- aggregate(new_fish_data[, c("FMR10C", "Z_ln_body_mass")], list(new_fish_data$species), sd, na.rm = TRUE)
spp_temp <- aggregate(new_fish_data[, c("temp_mean")], list(new_fish_data$species), mean, na.rm = TRUE)
spp_n <- aggregate(new_fish_data[, "n"], list(new_fish_data$species), length)

# Form into nice table 

spp_avgs <- data.frame(sciname = spp_means$Group.1,
                       n = spp_n$x,
                       fmr_mean = spp_means$FMR10C,
                       fmr_sd = spp_sds$FMR10C,
                       mass_mean = spp_means$Z_ln_body_mass,
                       mass_sd = spp_sds$Z_ln_body_mass,
                       temp_mean = spp_temp$x)
glimpse(spp_avgs)

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
  ylab(expression("ln oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
  scale_x_continuous(labels = c(0.0, 2.5, 5.0, 7.5, 10.0, 12.5),
                     name = expression("ln body mass (g)"),
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

png(here("plots",
         "mass_specific_body_mass_error_bars.png"), width = 7, height = 5,
    res = 250, units = "in")
bm_avgs_error
dev.off()
