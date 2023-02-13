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

new_fish_data <- read.csv(here("Outputs",
                               "07_C_resp_Models",
                               "Edited_Data",
                               "Scaling_Datasets",
                               "Scaling_Dataset.csv"))

new_fish_data$ln_O2_kg <- log(new_fish_data$O2_mg_kg)

# Load tree --------------------------------------------------------------------

new_fish_tree <- read.nexus(here("Outputs",
                                 "07_C_resp_Models",
                                 "Edited_Data",
                                 "Scaling_Datasets",
                                 "Scaling_Tree.nex"))

# Get relatedness matrix 

inv_phylo <- inverseA(new_fish_tree, nodes = "TIPS", scale = TRUE)
str(inv_phylo)
inv_phylo$Ainv

# Set priors

priors <- list(G = list(G1 = list(V = 1, nu = 0.02), # Prior for phylogenetic effects
                        G2 = list(V = 1, nu = 0.02), # Prior for n
                        G3 = list(V = 1, nu = 0.02), # Prior for source
                        G4 = list(V = 1, nu = 0.02)), # Prior for intraspecific effects
               R = list(V = 1, nu = 0.002)) # Prior for residuals
str(priors)

# Mass-specific ----------------------------------------------------------------

# Test model

mod_Test <- MCMCglmm(ln_O2_ind ~ Z_ln_BM + Z_inv_temp, 
                     random = ~ sciname + # phylogenetic effects
                       Species_ID + # intraspecific effects
                       n + # number of samples per data point
                       Source, # otolith source
                     family = "gaussian", 
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000, 
                     burnin = 1000, thin = 50)

# Look at summary

summary(mod_Test)

# Full model

mod_full <- MCMCglmm(ln_O2_ind ~ Z_ln_BM + Z_inv_temp, 
                     random = ~ sciname + # phylogenetic effects
                       Species_ID + # intraspecific effects
                       n + # number of samples per data point
                       Source, # otolith source
                     family = "gaussian", 
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000000, 
                     burnin = 10000, thin = 500);beep("mario")

saveRDS(mod_full, here("Outputs",
                       "07_C_resp_Models",
                       "Ch_6_Scaling",
                       "paper_models",
                       "whole_organism.rds"))

summary(mod_full)


# Whole-organism ---------------------------------------------------------------

# Test model

mod_Test <- MCMCglmm(ln_O2_kg ~ Z_ln_BM + Z_inv_temp, 
                     random = ~ sciname + # phylogenetic effects
                       Species_ID + # intraspecific effects
                       n + # number of samples per data point
                       Source, # otolith source
                     family = "gaussian", 
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000, 
                     burnin = 1000, thin = 50)

# Look at summary

summary(mod_Test)

# Full model

mod_mass_specific <- MCMCglmm(ln_O2_kg ~ Z_ln_BM + Z_inv_temp, 
                     random = ~ sciname + # phylogenetic effects
                       Species_ID + # intraspecific effects
                       n + # number of samples per data point
                       Source, # otolith source
                     family = "gaussian", 
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000000, 
                     burnin = 10000, thin = 500);beep("mario")

saveRDS(mod_mass_specific, here("Outputs",
                            "07_C_resp_Models",
                            "Ch_6_Scaling",
                            "paper_models",
                            "mass_specific.rds"))

summary(mod_mass_specific)

