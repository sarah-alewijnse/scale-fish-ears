#### Phylogeny Building ####

# Original code by Natalie Cooper

# Load libraries ---------------------------------------------------------------

library(fishtree)
library(ape)
library(tidyverse)
library(phytools)
library(here)

# Get phylogeny ----------------------------------------------------------------

# Read in species
fish <- read.csv(here("data",
                      "alewijnse_master_data.csv"))

# Add in underscore
fish <- mutate(fish, species = str_replace_all(species, " ", "_"))

# Extract chronogram (branch lengths = dates) of your species
fish_tree <- fishtree_phylogeny(species = fish$species, type = "chronogram")
str(fish_tree)

# Quick plot :)

plot(fish_tree, no.margin = TRUE, cex = 0.5, type = "fan")

# Run checks -------------------------------------------------------------------

is.rooted(fish_tree)
is.binary(fish_tree)
is.ultrametric(fish_tree)

# Force to be ultrametric
fish_tree <- force.ultrametric(fish_tree)
is.ultrametric(fish_tree)

# Get missing species ----------------------------------------------------------

missing <- setdiff(fish$species, fish_tree$tip.label)
missing

## Fix missing -----------------------------------------------------------------

# HMA - Hoplostethus mediterraneus
fish$species[fish$species == "Hoplostethus_mediterraneus"] <- "Hoplostethus_mediterraneus_mediterraneus"

# AGE - Arripis georgiana
fish$species[fish$species == "Arripis_georgiana"] <- "Arripis_georgianus"

# PAU - Pagrus auratus
fish$species[fish$species == "Pagrus_auratus"] <- "Chrysophrys_auratus"

# Substitute monophyletic species ----------------------------------------------

## IDO - Ijimaia doelfelini ====================================================
family_IDO <- "Ateleopodidae"
tree_IDO <- fishtree_phylogeny(rank = family_IDO)
plot(tree_IDO, no.margin = TRUE)

# Replace an Ijimaia with IDO
fish <- mutate(fish, species = str_replace(species, "Ijimaia_dofleini", "Ijimaia_loppei"))

## BDU - Bathypterois dubius ===================================================
family_BDU <- "Ipnopidae"
tree_BDU <- fishtree_phylogeny(rank = family_BDU)
plot(tree_BDU, no.margin = TRUE)

# Replace a Bathypterois with BDU
fish <- mutate(fish, species = str_replace(species, "Bathypterois_dubius", "Bathypterois_grallator"))

## BNI - Bathygadus nipponicus =================================================
family_BNI <- "Macrouridae"
tree_BNI <- fishtree_phylogeny(rank = family_BNI, type = "chronogram_mrca")
plot(tree_BNI, no.margin = TRUE, cex = 0.5, type = "fan")

# Replace Bathygadus with BNI
fish <- mutate(fish, species = str_replace(species, "Bathygadus_nipponicus", "Bathygadus_favosus"))

## HLE - Hymenocephalus lethonemus =============================================
family_HLE <- "Macrouridae"
tree_HLE <- fishtree_phylogeny(rank = family_HLE, type = "chronogram_mrca")
plot(tree_HLE, no.margin = TRUE, cex = 0.5, type = "fan")

# Replace a Hymenocephalus with HLE
fish <- mutate(fish, Species = str_replace(Species, "Hymenocephalus_lethonemus", "Hymenocephalus_italicus"))

# Remove those with no information ---------------------------------------------

# SGR
fish <- filter(fish, spp_code != "SGR")

# SKU
fish <- filter(fish, spp_code != "SKU")

# CLA
fish <- filter(fish, spp_code != "CLA")

# LCO
fish <- filter(fish, spp_code != "LCO")

# Add supporting additions from closely related species ------------------------

# LGR
fish$species[fish$species == "Lycodes_gracilis"] <- "Lycodes_palearis"

# CGL
fish$species[fish$species == "Cetonurus_globiceps"] <- "Ventrifossa_garmani"

# Update tree
fish_tree <- fishtree_phylogeny(species = fish$species, type = "chronogram")
str(fish_tree)
plot(fish_tree, cex = 0.5, type = "fan")

# Make ultrametric
fish_tree <- force.ultrametric(fish_tree)
is.ultrametric(fish_tree)

# View species
spp <- fish_tree$tip.label
spp

# Add species ------------------------------------------------------------------

## SLE - Saurida lessepsianus ==================================================
family_SLE <- "Synodontidae"
tree_SLE <- fishtree_phylogeny(rank = family_SLE, type = "chronogram_mrca")
plot(tree_SLE)

# Add next to SUN
fish_tree <- bind.tip(tree = fish_tree, 
                        tip.label = "Saurida_lessepsianus",
                        where = which(fish_tree$tip.label == "Saurida_undosquamis"), 
                      position = 0.5) # Use short branch length
plot(fish_tree, no.margin = T, cex = 0.5)

is.ultrametric(fish_tree)

## HMO - Helicolenus mouchezi ==================================================
family_HMO <- "Sebastidae"
tree_HMO <- fishtree_phylogeny(rank = family_HMO)
plot(tree_HMO, no.margin = TRUE)

# Add next to HDA
fish_tree <- bind.tip(fish_tree, 
                      tip.label = "Helicolenus_mouchezi",
                      where = which(fish_tree$tip.label == "Helicolenus_dactylopterus"), 
                      position = 0.5)
plot(fish_tree, no.margin = T, cex = 0.5) 

is.ultrametric(fish_tree)

## LGR - Lycodes gracilis ======================================================
family_LGR <- "Zoarcidae"
tree_LGR <- fishtree_phylogeny(rank = family_LGR, type = "chronogram_mrca")
plot(tree_LGR, no.margin = T, cex = 0.5)

# Add next to Lycodes palearis
fish_tree <- bind.tip(fish_tree, 
                      tip.label = "Lycodes_gracilis",
                      where = which(fish_tree$tip.label == "Lycodes_palearis"), 
                      position = 0.5)
plot(fish_tree, no.margin = T, cex = 0.5) 

# Remove Lycodes palearis
fish_tree <- drop.tip(fish_tree, "Lycodes_palearis")
plot(fish_tree, no.margin = T, cex = 0.5) 

is.ultrametric(fish_tree)

## CGL - Cetonurus globiceps ===================================================
family_CGL <- "Macrouridae"
tree_CGL <- fishtree_phylogeny(rank = family_CGL, type = "chronogram_mrca")
plot(tree_CGL, no.margin = T, cex = 0.5)

# Add next to Ventrifossa garmani
fish_tree <- bind.tip(fish_tree,
                      tip.label = "Cetonurus_globiceps",
                      where = which(fish_tree$tip.label == "Ventrifossa_garmani"),
                      position = 0.5)
plot(fish_tree, no.margin = T, cex = 0.5)

# Remove Ventrifossa garmani
fish_tree <- drop.tip(fish_tree, "Ventrifossa_garmani")
plot(fish_tree, no.margin = T, cex = 0.5) 

is.ultrametric(fish_tree)

# Put correct names back -------------------------------------------------------

# IDO
fish_tree$tip.label <- gsub(fish_tree$tip.label, 
                           pattern = "Ijimaia_loppei", 
                           replacement = "Ijimaia_dofleini")

# HLE
fish_tree$tip.label <- gsub(fish_tree$tip.label, 
                            pattern = "Hymenocephalus_italicus", 
                            replacement = "Hymenocephalus_lethonemus")

# BNI
fish_tree$tip.label <- gsub(fish_tree$tip.label, 
                            pattern = "Bathygadus_favosus", 
                            replacement = "Bathygadus_nipponicus")

# BDU
fish_tree$tip.label <- gsub(fish_tree$tip.label, 
                            pattern = "Bathypterois_grallator", 
                            replacement = "Bathypterois_dubius")

# View
plot(fish_tree, no.margin = T, cex = 0.7, type = "fan")

# Final checks -----------------------------------------------------------------

fish <- read.csv(here("data",
                      "alewijnse_master_data.csv"))
fish <- mutate(fish, species = str_replace_all(species, " ", "_"))

missing <- setdiff(fish$species, fish_tree$tip.label)
missing # Should only be seven: AGE, PAU and HME (spelling), and SGR, CLA, LCO and SKU (removed)

# Make ultrametric
is.ultrametric(fish_tree)
is.rooted(fish_tree)
is.binary(fish_tree)
plot(fish_tree, no.margin = T, cex = 0.7)

# Save -------------------------------------------------------------------------
write.nexus(fish_tree, file = here("outputs",
                                   "phylogeny",
                                   "scaling_fish_tree.nex"))
