#### MCMC GLMMs - C_resp - Body mass & Temperature - outer otos only ####

# Load packages

library(tidyverse)
library(ape)
library(MCMCglmm)
library(fishtree)
library(phytools)
library(geiger)
library(treeplyr)

#### Load data ####

C_resp_tidy <- read.csv("Outputs/07_C_resp_Models/Master_Data.csv")

#### Wrangle data ####

# Set source as a factor

C_resp_tidy$Source <- as.factor(C_resp_tidy$Source)
glimpse(C_resp_tidy)

# Add log BM

C_resp_tidy$log_BM <- log10(C_resp_tidy$Body_Mass)
glimpse(C_resp_tidy)

# Remove whole-otolith samples

C_resp_outer <- filter(C_resp_tidy, Sample_Type == "outer")

# Convert BM and temp to z-scores

load("~/PhD/GitHub/fish-ears/Functions/Z_score.Rdata")

Z_BM <- Z.score(C_resp_outer$Body_Mass, "Body_mass")
Z_log_BM <- Z.score(C_resp_outer$log_BM, "log_BM")
Z_temp <- Z.score(C_resp_outer$mean_Temp, "Temp")

C_resp_outer <- cbind(C_resp_outer, Z_BM, Z_log_BM, Z_temp)

glimpse(C_resp_outer)

# Histogram

ggplot(C_resp_outer, aes(x = Z_Body_mass)) +
  geom_histogram()

ggplot(C_resp_outer, aes(x = Z_log_BM)) +
  geom_histogram()

ggplot(C_resp_outer, aes(x = Z_Temp)) +
  geom_histogram()

ggplot(C_resp_outer, aes(x = mean_C_resp)) +
  geom_histogram()

#### Get the tree ####

# Extract chronogram (branch lengths = dates) of your species

fish_tree <- read.nexus("Phylogenies/My_Fish_Tree.nex")

# Check tree

is.rooted(fish_tree)
is.binary(fish_tree)
is.ultrametric(fish_tree)

# Check for dropped tips

check <- name.check(phy = fish_tree,
                    data = C_resp_outer,
                    data.names = C_resp_outer$sciname)
check$tree_not_data
check$data_not_tree

missing <- unique(check$data_not_tree)

# Trim the tree

fish_stuff <- make.treedata(tree = fish_tree,
                            data = C_resp_outer,
                            name_column = "sciname")
# Don't worry about "duplicated data" warning - we'll deal with that in a sec

# Save new fish tree

new_fish_tree <- fish_stuff$phy

# Remove dropped tips from data

new_fish_data <- C_resp_outer
new_fish_data <- dplyr::filter(new_fish_data, !sciname %in% missing)
glimpse(new_fish_data)

#### Get relatedness matrix ####

inv_phylo <- inverseA(new_fish_tree, nodes = "TIPS", scale = TRUE)
str(inv_phylo)
inv_phylo$Ainv

#### Linear model ####

# Set priors

priors <- list(G = list(G1 = list(V = 1, nu = 0.02), # Prior for phylogenetic effects
                        G2 = list(V = 1, nu = 0.02), # Prior for n
                        G3 = list(V = 1, nu = 0.02), # Prior for source
                        G4 = list(V = 1, nu = 0.02), # Prior for intraspecific effects
                        G5 = list(V = 1, fix = 1),
                        G6 = list(V = 1, fix = 1)), # Prior for standard deviation of C_resp (fixed)
               R = list(V = 1, nu = 0.002)) # Prior for residuals
str(priors)

# Test model

mod_Test <- MCMCglmm(mean_C_resp ~ Z_log_BM + Z_Temp, 
                     random = ~ sciname + # phylogenetic effects
                       Species_ID + # intraspecific effects
                       n + # number of samples per data point
                       Source + # otolith source
                       idh(se_C_resp):units + # standard deviations
                       idh(SE_Temp):units, # standard deviations
                     family = "gaussian", 
                     ginverse = list(sciname = inv_phylo$Ainv),
                     prior = priors, data = new_fish_data, nitt = 5000, 
                     burnin = 1000, thin = 50)

# Run diagnostics

autocorr(mod_Test$Sol) # Autocorrelation
plot(mod_Test$Sol) # Traceplot
# Geweke
geweke.diag(mod_Test$Sol)
geweke.plot(mod_Test$Sol)
effectiveSize(mod_Test$Sol) # Effective size
heidel.diag(mod_Test$Sol) # Heidelberg-Welch

# Look at summary

summary(mod_Test)

# Run full model

mod_BM_temp_outer <- MCMCglmm(mean_C_resp ~ Z_log_BM + Z_Temp, 
                              random = ~ sciname + # phylogenetic effects
                                Species_ID + # intraspecific effects
                                n + # number of samples per data point
                                Source + # otolith source
                                idh(se_C_resp):units + # standard deviations
                                idh(SE_Temp):units, # standard deviations
                              family = "gaussian", 
                              ginverse = list(sciname = inv_phylo$Ainv),
                              prior = priors, data = new_fish_data, nitt = 5000000, 
                              burnin = 10000, thin = 500)

autocorr(mod_BM_temp_outer$Sol) # Autocorrelation
plot(mod_BM_temp_outer$Sol) # Traceplot
# Geweke
geweke.diag(mod_BM_temp_outer$Sol)
geweke.plot(mod_BM_temp_outer$Sol)
effectiveSize(mod_BM_temp_outer$Sol) # Effective size
heidel.diag(mod_BM_temp_outer$Sol) # Heidelberg-Welch

summary(mod_BM_temp_outer$Sol)

saveRDS(mod_BM_temp_outer, file = "Outputs/07_C_resp_Models/MCMC_GLMM/01_Combined/BM_temp_outer.rds")
