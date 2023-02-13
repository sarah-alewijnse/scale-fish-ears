#### MCMC GLMMs - Oxygen consumption - Scaling ####

# Load packages
library(tidyverse)
library(ape)
library(MCMCglmm)
library(fishtree)
library(phytools)
library(geiger)
library(treeplyr)
library(beepr)
library(here)

# Load data --------------------------------------------------------------------

source(here("code",
            "00_pre_processing",
            "03_oxygen_consumption",
            "02_scaling_dataset_main.R"))
# Ignore duplicated data warning

# Get relatedness matrix 
inv_phylo <- inverseA(new_fish_tree, nodes = "TIPS", scale = TRUE)
str(inv_phylo)
inv_phylo$Ainv

# Set priors -------------------------------------------------------------------

priors <- list(G = list(G1 = list(V = 1, nu = 1), # Prior for phylogenetic effects
                        G2 = list(V = 1, nu = 1), # Prior for n
                        G3 = list(V = 1, nu = 1), # Prior for source
                        G4 = list(V = 1, nu = 1)), # Prior for intraspecific effects
               R = list(V = 1, nu = 1)) # Prior for residuals
str(priors) # Check

# Mass-specific ----------------------------------------------------------------

# Test model
mod_test <- MCMCglmm(ln_o2_mg_kg ~ Z_ln_body_mass + Z_inv_temp_mean, 
                     random = ~ species + # phylogenetic effects
                       species_id + # intraspecific effects
                       n + # number of samples per data point
                       source, # otolith source
                     family = "gaussian", 
                     ginverse = list(species = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000, 
                     burnin = 1000, thin = 50)

# Look at summary
summary(mod_test)

# Full model
mod_full <- MCMCglmm(ln_o2_mg_kg ~ Z_ln_body_mass + Z_inv_temp_mean, 
                     random = ~ species + # phylogenetic effects
                       species_id + # intraspecific effects
                       n + # number of samples per data point
                       source, # otolith source
                     family = "gaussian",
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000000, 
                     burnin = 10000, thin = 500);beep("mario")

# Save
saveRDS(mod_full, here("outputs",
                       "models",
                       "main",
                       "mass_specific.rds"))

# Check summary
summary(mod_full)

# Whole-organism ---------------------------------------------------------------

# Test model
mod_test <- MCMCglmm(ln_o2_ind ~ Z_ln_body_mass + Z_inv_temp_mean, 
                     random = ~ species + # phylogenetic effects
                       species_id + # intraspecific effects
                       n + # number of samples per data point
                       source, # otolith source
                     family = "gaussian", 
                     ginverse = list(species = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000, 
                     burnin = 1000, thin = 50)

# Look at summary
summary(mod_test)

# Full model
mod_full <- MCMCglmm(ln_o2_ind ~ Z_ln_body_mass + Z_inv_temp_mean, 
                     random = ~ species + # phylogenetic effects
                       species_id + # intraspecific effects
                       n + # number of samples per data point
                       source, # otolith source
                     family = "gaussian", 
                     ginverse = list(species = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000000, 
                     burnin = 10000, thin = 500);beep("mario")

# Save
saveRDS(mod_full, here("outputs",
                       "models",
                       "main",
                       "whole_organism.rds"))

# Check summary
summary(mod_full)
