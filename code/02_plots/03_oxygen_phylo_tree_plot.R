#### Phylogeny plots ####

library(tidyverse)
library(MCMCglmm)
library(ape)
library(geiger)
library(treeplyr)
library(fishtree)
library(phytools)
library(caper)
library(ggtree)
library(RColorBrewer)
library(viridis)
library(here)
library(ggtreeExtra)
library(cowplot)

# Load new_fish_data -----------------------------------------------------------

source(here("code",
            "00_pre_processing",
            "03_oxygen_consumption",
            "02_scaling_dataset_main.R"))
# Ignore duplicated new_fish_data warning

# Standardise to common body mass and temp -------------------------------------

# Set standard body mass
standard_ln_body_mass <- log(300)
mean_ln_body_mass <- mean(new_fish_data$ln_body_mass)
sd_ln_body_mass <- sd(new_fish_data$ln_body_mass)
standard_Z_ln_body_mass <- (standard_ln_body_mass - mean_ln_body_mass) / sd_ln_body_mass

# Standard temp
boltz <- 8.62*10^-5 # Boltzmann constant
mean_temp <- mean(new_fish_data$inv_temp_mean)
sd_temp <- sd(new_fish_data$inv_temp_mean)
standard_inv_temp <- 1/((273.15+10) * boltz) # inv_temp = 1/kT (T in kelvin)
standard_Z_inv_temp <- (standard_inv_temp - mean_temp) / sd_temp # Standardize

# Load model -------------------------------------------------------------------

bmt_mod <- readRDS(here("outputs",
                        "models",
                        "main",
                        "mass_specific.rds"))

# Get parameters
bmt_mod_sum <- summary(bmt_mod$Sol)
a <- bmt_mod_sum$statistics[1]
b_ln_BM <- bmt_mod_sum$statistics[2]
b_inv_temp <- bmt_mod_sum$statistics[3]

# Get residuals ----------------------------------------------------------------

# NB: Residuals function doesn't work

# Get fitted values

new_fish_data$fitted_O2 <- a + b_ln_BM * new_fish_data$Z_ln_body_mass + b_inv_temp * new_fish_data$Z_inv_temp_mean
hist(new_fish_data$fitted_O2) # Histogram

# Get residuals
new_fish_data$resid_O2 <- new_fish_data$ln_o2_mg_kg - new_fish_data$fitted_O2
hist(new_fish_data$resid_O2) # Histogram
plot(new_fish_data$fitted_O2, new_fish_data$resid_O2) # Fitted vs. residual

# Normalize --------------------------------------------------------------------

new_fish_data$normalised_ln_O2 <- a + b_ln_BM * standard_Z_ln_body_mass + b_inv_temp * standard_Z_inv_temp + new_fish_data$resid_O2
hist(new_fish_data$normalised_ln_O2)

plot(normalised_ln_O2 ~ log(body_mass), data = new_fish_data)

# Mean oxygen consumption ------------------------------------------------------

mean_o2 <- dplyr::select(new_fish_data,
                         species,
                         normalised_ln_O2) %>% 
  dplyr::group_by(species) %>%
  summarise(mean = mean(normalised_ln_O2))

mean_o2 <- as.data.frame(mean_o2)
rownames(mean_o2) <- mean_o2$species

# Tree plot --------------------------------------------------------------------

tree_plot_base <- ggtree(new_fish_tree, layout = "fan", lwd = 1) %<+% mean_o2 + xlim(0, 200)
tree_plot_base

tree_plot_lab <- ggtree(new_fish_tree, layout = "fan") +
  geom_tiplab()
tree_plot_lab

# Add bars
MRCA(new_fish_tree, .node1 = "Beryx_splendens", .node2 = "Hoplostethus_atlanticus") # Beryciformes
MRCA(new_fish_tree, .node1 = "Triglops_nybelini", .node2 = "Champsocephalus_gunnari") # Perciformes
MRCA(new_fish_tree, .node1 = "Solea_solea", .node2 = "Hippoglossoides_platessoides") # Pleuronectiformes
MRCA(new_fish_tree, .node1 = "Ateleopus_japonicus", .node2 = "Ijimaia_dofleini") # Ateleopodiformes
MRCA(new_fish_tree, .node1 = "Seriola_dumerili", .node2 = "Coryphaena_hippurus") # Carangiformes
MRCA(new_fish_tree, .node1 = "Coryphaenoides_rupestris", .node2 = "Mora_moro") # Gadiformes
MRCA(new_fish_tree, .node1 = "Arripis_georgianus", .node2 = "Thunnus_maccoyii") # Scombrifomres


o2_tree <- ggtree(new_fish_tree, layout = "fan", lwd = 0.5) + 
  xlim(0, 280) +
  geom_fruit(data = mean_o2,
             geom = geom_bar,
             pwidth = 0.4,
             mapping = aes(y = species,
                           x = mean,
                           fill = mean),
             orientation = "y",
             stat = "identity") +
  scale_fill_viridis(option = "mako", 
                     end = 0.9,
                     name = "") +
  theme(text = element_text(size = 10),
        legend.position = "bottom") + 
  geom_cladelabel(node = 180, label = "", offset = 60, barsize = 1) + # Gadiformes
  geom_cladelabel(node = 178, label = "", offset = 60, barsize = 1) + # Beryciformes
  geom_cladelabel(node = 141, label = "", offset = 70, barsize = 1) + # Scombriformes
  geom_cladelabel(node = 219, label = "", offset = 60, barsize = 1) + # Ateleopodiformes
  geom_cladelabel(node = 135, label = "", offset = 70, barsize = 1) + # Carangiformes
  geom_cladelabel(node = 153, label = "", offset = 70, barsize = 1) + # Perciformes
  geom_cladelabel(node = 128, label = "", offset = 60, barsize = 1) # Pleuronectiformes

o2_tree

svg(here("plots",
         "oxygen_tree.svg"), width = 10, height = 10)
o2_tree
dev.off()
