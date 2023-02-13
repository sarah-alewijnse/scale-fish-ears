#### Phylogeny Building ####

# Original code by Natalie Cooper

# Load libraries

library(fishtree)
library(ape)
library(tidyverse)
library(phytools)

# Read species in list

fish <- read.csv("Data/04_Species_Metadata/Ecology_Data/Ecol_Table/Ecol_Table_Final.csv")

# Add in underscore

fish <- mutate(fish, Species = str_replace_all(Species, " ", "_"))

# Extract chronogram (branch lengths = dates) of your species

fish_tree <- fishtree_phylogeny(species = fish$Species, type = "chronogram")
str(fish_tree)

# Quick plot :)

plot(fish_tree, no.margin = TRUE, cex = 0.5)

# Run checks

is.rooted(fish_tree)
is.binary(fish_tree)
is.ultrametric(fish_tree)

fish_tree <- force.ultrametric(fish_tree)
is.ultrametric(fish_tree)

# Get missing species

missing <- setdiff(fish$Species, fish_tree$tip.label)
missing

#### Spelling Error ####

## HMA - Hoplostethus mediterraneus

fish$Species[fish$Species == "Hoplostethus_mediterraneus"] <- "Hoplostethus_mediterraneus_mediterraneus"

## AGE - Arripis georgiana

fish$Species[fish$Species == "Arripis_georgiana"] <- "Arripis_georgianus"

## PAU - Pagrus auratus

fish$Species[fish$Species == "Pagrus_auratus"] <- "Chrysophrys_auratus"

#### Monophyletic ####

## IDO - Ijimaia doelfelini

family_IDO <- fish %>% filter(Species == "Ijimaia_dofleini") %>% pull(Family)
tree_IDO <- fishtree_phylogeny(rank = family_IDO)
plot(tree_IDO, no.margin = TRUE)

# Replace an Ijimaia with IDO

fish <- mutate(fish, Species = str_replace(Species, "Ijimaia_dofleini", "Ijimaia_loppei"))

## BDU - Bathypterois dubius

family_BDU <- fish %>% filter(Species == "Bathypterois_dubius") %>% pull(Family)
tree_BDU <- fishtree_phylogeny(rank = family_BDU)
plot(tree_BDU, no.margin = TRUE)

# Replace a Bathypterois with BDU

fish <- mutate(fish, Species = str_replace(Species, "Bathypterois_dubius", "Bathypterois_grallator"))

## BNI - Bathygadus nipponicus

family_BNI <- fish %>% filter(Species == "Bathygadus_nipponicus") %>% pull(Family)
tree_BNI <- fishtree_phylogeny(rank = family_BNI, type = "chronogram_mrca")
plot(tree_BNI, no.margin = TRUE, cex = 0.5, type = "fan")

fish <- mutate(fish, Species = str_replace(Species, "Bathygadus_nipponicus", "Bathygadus_favosus"))

## HLE - Hymenocephalus lethonemus

family_HLE <- fish %>% filter(Species == "Hymenocephalus_lethonemus") %>% pull(Family)
tree_HLE <- fishtree_phylogeny(rank = family_HLE, type = "chronogram_mrca")
plot(tree_HLE, no.margin = TRUE, cex = 0.5, type = "fan")

fish <- mutate(fish, Species = str_replace(Species, "Hymenocephalus_lethonemus", "Hymenocephalus_italicus"))

#### No information ####

## Remove SGR 

fish <- filter(fish, Spp_Code != "SGR")

## Remove SKU

fish <- filter(fish, Spp_Code != "SKU")

## Remove CLA

fish <- filter(fish, Spp_Code != "CLA")

## Remove LCO

fish <- filter(fish, Spp_Code != "LCO")

#### Supporting additions ####

fish$Species[fish$Species == "Lycodes_gracilis"] <- "Lycodes_palearis"

fish$Species[fish$Species == "Cetonurus_globiceps"] <- "Ventrifossa_garmani"

#### Update tree ####

fish_tree <- fishtree_phylogeny(species = fish$Species, type = "chronogram")
str(fish_tree)
plot(fish_tree, cex = 0.5)

# Make ultrametric

fish_tree <- force.ultrametric(fish_tree)
is.ultrametric(fish_tree)

# View species

spp <- fish_tree$tip.label
view(spp)

#### Adding species ####

## SLE - Saurida lessepsianus

family_SLE <- fish %>% filter(Species == "Saurida_lessepsianus") %>% pull(Family)
tree_SLE <- fishtree_phylogeny(rank = family_SLE, type = "chronogram_mrca")
plot(tree_SLE)

# Add next to SUN

fish_tree <- bind.tip(tree = fish_tree, 
                        tip.label = "Saurida_lessepsianus",
                        where = which(fish_tree$tip.label == "Saurida_undosquamis"), 
                      position = 0.5) # Use short branch length
plot(fish_tree, no.margin = T, cex = 0.5)

is.ultrametric(fish_tree)

## HMO - Helicolenus mouchezi

family_HMO <- fish %>% filter(Species == "Helicolenus_mouchezi") %>% pull(Family)
tree_HMO <- fishtree_phylogeny(rank = family_HMO)

plot(tree_HMO, no.margin = TRUE)

# Add next to HDA

fish_tree <- bind.tip(fish_tree, 
                      tip.label = "Helicolenus_mouchezi",
                      where = which(fish_tree$tip.label == "Helicolenus_dactylopterus"), 
                      position = 0.5)
plot(fish_tree, no.margin = T, cex = 0.5) 

is.ultrametric(fish_tree)

## LGR - Lycodes gracilis

family_LGR <- fish %>% filter(Species == "Lycodes_gracilis") %>% pull(Family)
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

## CGL - Cetonurus globiceps

family_CGL <- fish %>% filter(Species == "Lycodes_gracilis") %>% pull(Family)
tree_CGL <- fishtree_phylogeny(rank = family_CGL, type = "chronogram_mrca")
plot(tree_CGL, no.margin = T, cex = 0.5)

fish_tree <- bind.tip(fish_tree,
                      tip.label = "Cetonurus_globiceps",
                      where = which(fish_tree$tip.label == "Ventrifossa_garmani"),
                      position = 0.5)
plot(fish_tree, no.margin = T, cex = 0.5)

# Remove Ventrifossa garmani

fish_tree <- drop.tip(fish_tree, "Ventrifossa_garmani")
plot(fish_tree, no.margin = T, cex = 0.5) 

is.ultrametric(fish_tree)

#### Put correct names back ####

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

#### Final checks ####

fish <- read.csv("Data/04_Species_Metadata/Ecology_Data/Ecol_Table/Ecol_Table_Final.csv")
fish <- mutate(fish, Species = str_replace_all(Species, " ", "_"))

missing <- setdiff(fish$Species, fish_tree$tip.label)
view(missing)
  # Should only be seven: AGE, PAU and HME (spelling), and SGR, CLA, LCO and SKU (removed)

# Make ultrametric

is.ultrametric(fish_tree)
is.rooted(fish_tree)
is.binary(fish_tree)
plot(fish_tree, no.margin = T, cex = 0.7)

#### Save ####

write.nexus(fish_tree, file = "Phylogenies/My_Fish_Tree.nex")
