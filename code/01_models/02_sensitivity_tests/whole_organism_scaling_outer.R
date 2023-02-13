#### MCMC GLMMs - C_resp - Scaling - ours ####

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
                               "scaling_data_outer.csv"))

# Load tree --------------------------------------------------------------------

new_fish_tree <- read.nexus(here("Outputs",
                                 "07_C_resp_Models",
                                 "Edited_Data",
                                 "Scaling_Datasets",
                                 "scaling_tree_outer.nex"))

# Get relatedness matrix 

inv_phylo <- inverseA(new_fish_tree, nodes = "TIPS", scale = TRUE)
str(inv_phylo)
inv_phylo$Ainv

# Set priors

priors <- list(G = list(G1 = list(V = 1, nu = 0.02), # Prior for phylogenetic effects
                        G2 = list(V = 1, nu = 0.02), # Prior for sciname
                        G3 = list(V = 1, nu = 0.02), # Prior for intraspecific effects
                        G4 = list(V = 1, fix = 1), # Prior for standard error of C_resp (fixed)
                        G5 = list(V = 1, fix = 1)), # Prior for standard error of temp (fixed)
               R = list(V = 1, nu = 0.002)) # Prior for residuals
str(priors) # Check

# Test model -------------------------------------------------------------------

mod_Test <- MCMCglmm(ln_O2_ind ~ Z_ln_BM + Z_inv_temp, 
                     random = ~ sciname + # phylogenetic effects
                       Species_ID + # intraspecific effects
                       Source + # otolith source
                       idh(se_C_resp):units + # C_resp uncertianty
                       idh(SE_Temp):units, # Temp uncertainty
                     family = "gaussian", 
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000, 
                     burnin = 1000, thin = 50)

summary(mod_Test)

# Full model -------------------------------------------------------------------

mod_full <- MCMCglmm(ln_O2_ind ~ Z_ln_BM + Z_inv_temp, 
                     random = ~ sciname + # phylogenetic effects
                       Species_ID + # intraspecific effects
                       Source + # otolith source
                       idh(se_C_resp):units + # C_resp uncertianty
                       idh(SE_Temp):units, # Temp uncertainty
                     family = "gaussian", 
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000000, 
                     burnin = 10000, thin = 500);beep("mario")


saveRDS(mod_full, here("Outputs",
                       "07_C_resp_Models",
                       "Ch_6_Scaling",
                       "Main_Models",
                       "whole_organism_scaling_outer.rds"))

summary(mod_full)
