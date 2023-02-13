#### Mass corrected O2 consumption vs. temperature ####

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

# Normalise to common body mass ------------------------------------------------

# Fit O2 consumption based on body mass alone
data$fitted_o2_mg_kg <- a + b_ln_bm * data$Z_ln_BM
plot(data$Z_inv_temp, data$fitted_o2_mg_kg)
plot(data$ln_BM, data$fitted_o2_mg_kg)

# Get residuals
data$resid_o2_ind <- data$ln_O2_mg_kg - data$fitted_o2_mg_kg
plot(data$fitted_o2_mg_kg, data$resid_o2_ind)

# Get mean and sd body mass
mean_bm <- mean(data$ln_BM)
sd_bm <- sd(data$ln_BM)

# Convert 300g to z-score
standard_Z_ln_bm <- (log(300) - mean_bm) / sd_bm # Standardize

# Calculate normalised FMR to 300g
data$FMR300g <- a + b_ln_bm * standard_Z_ln_bm + data$resid_o2_ind
plot(data$inv_Temp, data$FMR300g)
range(data$FMR300g)
range(data$ln_O2_mg_kg)

# Get ci -----------------------------------------------------------------------

# Load function
load(here("Functions",
          "LM_CI.Rdata"))

# Calclate
temp_ci <- lm.ci(data$Z_inv_temp,
               a = a, b = b_inv_temp,
               a_low_ci = o2_mod_sum$quantiles[1, 1], 
               b_low_ci = o2_mod_sum$quantiles[3, 1],
               a_up_ci = o2_mod_sum$quantiles[1, 5], 
               b_up_ci = o2_mod_sum$quantiles[3, 5])

# Assign required values to objects --------------------------------------------
# Makes things tidier

min_inv_temp <- min(data$Z_inv_temp)
max_inv_temp <- max(data$Z_inv_temp)
mean_inv_temp <- mean(data$inv_Temp)
sd_inv_temp <- sd(data$inv_Temp)

b_UTD <- -0.65*sd_inv_temp # Times by sd to convert to z-scored scale

# Get species averages ---------------------------------------------------------

# Get means, sd and n

spp_means <- aggregate(data[, c("FMR300g", "Z_inv_temp")], list(data$sciname), mean, na.rm = TRUE)
spp_sds <- aggregate(data[, c("FMR300g", "Z_inv_temp")], list(data$sciname), sd, na.rm = TRUE)
spp_mass <- aggregate(data[, c("ln_BM")], list(data$sciname), mean, na.rm = TRUE)
spp_n <- aggregate(data[, "n"], list(data$sciname), length)

# Form into nice table 

spp_avgs <- data.frame(sciname = spp_means$Group.1,
                       n = spp_n$x,
                       fmr_mean = spp_means$FMR300g,
                       fmr_sd = spp_sds$FMR300g,
                       temp_mean = spp_means$Z_inv_temp,
                       temp_sd = spp_sds$Z_inv_temp,
                       mass_mean = spp_mass$x)
glimpse(spp_avgs)

# Averages with size = n -------------------------------------------------------

temp_avgs_n <- ggplot(data = spp_avgs, 
                    aes(x = temp_mean, y = fmr_mean)) +
  geom_polygon(data = temp_ci, 
               aes(x_ribbon, y_ribbon), 
               alpha = 0.3) +
  geom_point(aes(size = n, col = mass_mean), 
             alpha = 0.7) +
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
  scale_color_viridis(option = "mako", 
                      name = "log body mass (g)",
                      begin = 0.2,
                      end = 0.9) +
  # Customise the theme
  ylab(expression("log oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
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
temp_avgs_n

svg(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "mass_specific_comb_spp_avgs_n_temp.svg"), width = 7, height = 5)
temp_avgs_n
dev.off()

# Averages with size = sd -------------------------------------------------------

spp_avgs[is.na(spp_avgs)] <- 0

temp_avgs_sd <- ggplot(data = spp_avgs, 
                      aes(x = temp_mean, y = fmr_mean)) +
  geom_polygon(data = temp_ci, 
               aes(x_ribbon, y_ribbon), 
               alpha = 0.3) +
  geom_point(aes(size = fmr_sd, col = mass_mean), 
             alpha = 0.7) +
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
  scale_color_viridis(option = "mako", 
                      name = "log body mass (g)",
                      begin = 0.2,
                      end = 0.9) +
  # Customise the theme
  ylab(expression("log oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
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
temp_avgs_sd

svg(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "mass_specific_comb_spp_avgs_sd_temp.svg"), width = 7, height = 5)
temp_avgs_sd
dev.off()

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
                      name = expression("log"[10]~"body mass (g)")) +
  # Customise the theme
  ylab(expression("log"[10]~"oxygen consumption (mg O "[2]~"kg"^-1~"h"^-1~")")) +
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

png(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "mass_specific_comb_spp_avgs_error_bars_temp.png"), width = 7, height = 5,
    res = 250, units = "in")
temp_avgs_error
dev.off()
