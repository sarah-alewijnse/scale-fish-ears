#### Mass corrected O2 consumption vs. temperature ####

library(tidyverse)
library(MCMCglmm)
library(here)
library(viridis)

# Load new_fish_data -----------------------------------------------------------

source(here("code",
            "00_pre_processing",
            "03_oxygen_consumption",
            "02_scaling_dataset_main.R"))
# Ignore duplicated new_fish_data warning

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

# Normalise to common body mass ------------------------------------------------

# Fit O2 consumption based on body mass alone
new_fish_data$fitted_o2_mg_kg <- a + b_ln_bm * new_fish_data$Z_ln_body_mass
plot(new_fish_data$Z_inv_temp, new_fish_data$fitted_o2_mg_kg)
plot(new_fish_data$ln_body_mass, new_fish_data$fitted_o2_mg_kg)

# Get residuals
new_fish_data$resid_o2_ind <- new_fish_data$ln_o2_mg_kg - new_fish_data$fitted_o2_mg_kg
plot(new_fish_data$fitted_o2_mg_kg, new_fish_data$resid_o2_ind)

# Get mean and sd body mass
mean_bm <- mean(new_fish_data$ln_body_mass)
sd_bm <- sd(new_fish_data$ln_body_mass)

# Convert 300g to z-score
standard_Z_ln_bm <- (log(300) - mean_bm) / sd_bm # Standardize

# Calculate normalised FMR to 300g
new_fish_data$FMR300g <- a + b_ln_bm * standard_Z_ln_bm + new_fish_data$resid_o2_ind
plot(new_fish_data$inv_temp_mean, new_fish_data$FMR300g)
range(new_fish_data$FMR300g)
range(new_fish_data$ln_o2_mg_kg)

# Get ci -----------------------------------------------------------------------

# Load function
load(here("functions",
          "LM_CI.Rdata"))

# Calclate
temp_ci <- lm.ci(new_fish_data$Z_inv_temp_mean,
               a = a, b = b_inv_temp,
               a_low_ci = o2_mod_sum$quantiles[1, 1], 
               b_low_ci = o2_mod_sum$quantiles[3, 1],
               a_up_ci = o2_mod_sum$quantiles[1, 5], 
               b_up_ci = o2_mod_sum$quantiles[3, 5])

# Assign required values to objects --------------------------------------------
# Makes things tidier

min_inv_temp <- min(new_fish_data$Z_inv_temp_mean)
max_inv_temp <- max(new_fish_data$Z_inv_temp_mean)
mean_inv_temp <- mean(new_fish_data$inv_temp_mean)
sd_inv_temp <- sd(new_fish_data$inv_temp_mean)

b_UTD <- -0.65*sd_inv_temp # Times by sd to convert to z-scored scale

# Get species averages ---------------------------------------------------------

# Get means, sd and n

spp_means <- aggregate(new_fish_data[, c("FMR300g", "Z_inv_temp_mean")], list(new_fish_data$species), mean, na.rm = TRUE)
spp_sds <- aggregate(new_fish_data[, c("FMR300g", "Z_inv_temp_mean")], list(new_fish_data$species), sd, na.rm = TRUE)
spp_mass <- aggregate(new_fish_data[, c("ln_body_mass")], list(new_fish_data$species), mean, na.rm = TRUE)
spp_n <- aggregate(new_fish_data[, "n"], list(new_fish_data$species), length)

# Form into nice table 

spp_avgs <- data.frame(sciname = spp_means$Group.1,
                       n = spp_n$x,
                       fmr_mean = spp_means$FMR300g,
                       fmr_sd = spp_sds$FMR300g,
                       temp_mean = spp_means$Z_inv_temp,
                       temp_sd = spp_sds$Z_inv_temp,
                       mass_mean = spp_mass$x)
glimpse(spp_avgs)

# Error bars -------------------------------------------------------------------

temp_avgs_error <- ggplot(data = spp_avgs, 
                       aes(x = temp_mean, y = fmr_mean)) +
  geom_polygon(data = temp_ci, 
               aes(x_ribbon, y_ribbon), 
               alpha = 0.3) +
  geom_point(aes(col = mass_mean), 
             alpha = 1, size = 3) +
  geom_errorbar(aes(ymin = fmr_mean - fmr_sd,
                    ymax = fmr_mean + fmr_sd,
                    col = mass_mean),
                alpha = 0.4) +
  geom_errorbarh(aes(xmin = temp_mean - temp_sd,
                     xmax = temp_mean + temp_sd,
                     col = mass_mean),
                 alpha = 0.4) +
  geom_segment(aes(x = min_inv_temp,
                   xend = max_inv_temp,
                   y = a + b_inv_temp * min_inv_temp,
                   yend = a + b_inv_temp * max_inv_temp), lwd = 0.5) + # IDK why this is the wrong way round it just is
  # geom_segment(aes(x = min_inv_temp, xend = max_inv_temp, 
  #                  y = a + b_UTD * min_inv_temp, 
  #                  yend = a + b_UTD * max_inv_temp,
  #                  linetype = "solid"), lwd = 0.5) +
  # scale_linetype_manual(name = "Alternative \nactivation energies",
  #                       values = c("solid", "dotted"),
  # labels = c("0.26", "0.65")) +
  scale_color_viridis(option = "viridis", 
                      name = expression("ln body mass (g)")) +
  # Customise the theme
  ylab(expression("ln oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
  scale_x_continuous(labels = c(38, 39, 40, 41, 42, 43),
                     name = "Inverse temperature (1/kT)",
                     c((38 - mean_inv_temp) / sd_inv_temp,
                       (39 - mean_inv_temp) / sd_inv_temp,
                       (40 - mean_inv_temp) / sd_inv_temp,
                       (41 - mean_inv_temp) / sd_inv_temp,
                       (42 - mean_inv_temp) / sd_inv_temp,
                       (43 - mean_inv_temp) / sd_inv_temp)) +
  scale_y_continuous(breaks = c(2, 3, 4, 5, 6, 7),
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
temp_avgs_error

png(here("plots",
         "mass_specific_temp_error_bars.png"), width = 7, height = 5,
    res = 250, units = "in")
temp_avgs_error
dev.off()
