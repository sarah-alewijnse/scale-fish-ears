#### Scaling dataset creator - ours only ####

# Load packages

library(tidyverse)
library(ape)
library(MCMCglmm)
library(fishtree)
library(phytools)
library(geiger)
library(treeplyr)
library(here)
library(beepr)

# Load data --------------------------------------------------------------------

oxy_data <- read.csv(here("outputs",
                          "oxygen_consumption",
                          "alewijnse_master_data_o2_consumption.csv"))

# Select only our data
oxy_data <- filter(oxy_data, group == "Ours")

# Get log body mass ------------------------------------------------------------
oxy_data$ln_body_mass <- log(oxy_data$body_mass)

# Get arrhenius temperature ----------------------------------------------------

# Boltzmann constant
boltz <- 8.62*10^-5

# Temp in kelvin
oxy_data$temp_kelvin_mean <- oxy_data$temp_mean + 273.15

# Inverse temp
oxy_data$inv_temp_mean <- 1/(oxy_data$temp_kelvin_mean * boltz)

# Get ln oxygen consumption ----------------------------------------------------

oxy_data$ln_o2_ind <- log(oxy_data$o2_mg_ind)

oxy_data$ln_o2_mg_kg <- log(oxy_data$o2_mg_kg)

# Get z scores -----------------------------------------------------------------

# Load function
load(here("functions",
          "z_score.Rdata"))

# log body mass
Z_ln_body_mass <- Z.score(oxy_data$ln_body_mass, "ln_body_mass")

# Inverse temp
Z_inv_temp <- Z.score(oxy_data$inv_temp_mean, "inv_temp_mean")

# Bind
oxy_data <- cbind(oxy_data, Z_ln_body_mass, Z_inv_temp)

# Get species ID ---------------------------------------------------------------

oxy_data$species_id <- oxy_data$species

# Get the tree -----------------------------------------------------------------

# Load tree

fish_tree <- read.nexus(here("outputs",
                             "phylogeny",
                             "fish_tree.nex"))

# Check tree (should all be true)
is.rooted(fish_tree)
is.binary(fish_tree)
is.ultrametric(fish_tree)

# Check for dropped tips
check <- name.check(phy = fish_tree,
                    data = oxy_data,
                    data.names = oxy_data$species)

# Tree not data
check$tree_not_data

# Will be several lit species

# Data not tree
missing <- unique(check$data_not_tree)
missing

# Should be (no phylogenies):
# Lepidoperca coatsii
# Coelorinchus labiatus
# Spectrunculus grandis

# Trim the tree
fish_stuff <- make.treedata(tree = fish_tree,
                            data = oxy_data,
                            name_column = "species")
# Don't worry about "duplicated data" warning 

# Save new fish tree
new_fish_tree <- fish_stuff$phy

# Remove dropped tips from data
new_fish_data <- oxy_data
new_fish_data <- dplyr::filter(new_fish_data, !species %in% missing)
glimpse(new_fish_data)
