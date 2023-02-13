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

# Load data --------------------------------------------------------------------

new_fish_data <- read.csv(here("Outputs",
                               "07_C_resp_Models",
                               "Edited_Data",
                               "Normalized_C_resp.csv"))

# Load tree --------------------------------------------------------------------

new_fish_tree <- read.nexus(here("Outputs",
                                 "07_C_resp_Models",
                                 "Edited_Data",
                                 "Scaling_Datasets",
                                 "Scaling_Tree.nex"))

# Mean oxygen consumption ------------------------------------------------------

mean_o2 <- dplyr::select(new_fish_data,
                         sciname,
                         normalised_ln_O2) %>% 
  dplyr::group_by(sciname) %>%
  summarise(mean = mean(normalised_ln_O2))

mean_o2 <- as.data.frame(mean_o2)
rownames(mean_o2) <- mean_o2$sciname

# Tree plot --------------------------------------------------------------------

tree_plot_base <- ggtree(new_fish_tree, layout = "fan", lwd = 1) %<+% mean_o2 + xlim(0, 200)
tree_plot_base

tree_plot_lab <- ggtree(new_fish_tree, layout = "fan") +
  geom_tiplab()
tree_plot_lab

# Add bars

MRCA(new_fish_tree, .node1 = "Beryx_splendens", .node2 = "Hoplostethus_atlanticus")
MRCA(new_fish_tree, .node1 = "Lethrinus_lentjan", .node2 = "Pagellus_erythrinus")
MRCA(new_fish_tree, .node1 = "Triglops_nybelini", .node2 = "Champsocephalus_gunnari")
MRCA(new_fish_tree, .node1 = "Electrona_antarctica", .node2 = "Gymnoscopelus_nicholsi")

MRCA(new_fish_tree, .node1 = "Mullus_barbatus_barbatus", .node2 = "Upeneus_moluccensis")
MRCA(new_fish_tree, .node1 = "Solea_solea", .node2 = "Hippoglossoides_platessoides")
MRCA(new_fish_tree, .node1 = "Ateleopus_japonicus", .node2 = "Ijimaia_dofleini")
MRCA(new_fish_tree, .node1 = "Makaira_nigricans", .node2 = "Kajikia_albida")
MRCA(new_fish_tree, .node1 = "Seriola_dumerili", .node2 = "Coryphaena_hippurus")
MRCA(new_fish_tree, .node1 = "Engraulis_encrasicolus", .node2 = "Sprattus_sprattus")
MRCA(new_fish_tree, .node1 = "Coryphaenoides_rupestris", .node2 = "Mora_moro")
MRCA(new_fish_tree, .node1 = "Arripis_georgianus", .node2 = "Thunnus_maccoyii")


o2_tree <- ggtree(new_fish_tree, layout = "fan", lwd = 0.5) + 
  xlim(0, 280) +
  geom_fruit(data = mean_o2,
             geom = geom_bar,
             pwidth = 0.4,
             mapping = aes(y = sciname,
                           x = mean,
                           fill = mean),
             orientation = "y",
             stat = "identity") +
  scale_fill_viridis(option = "mako", 
                     end = 0.9,
                     name = "") +
  theme(text = element_text(size = 10),
        legend.position = "bottom") + 
  geom_cladelabel(node = 179, label = "", offset = 60, barsize = 1) + # Gadiformes
  geom_cladelabel(node = 177, label = "", offset = 60, barsize = 1) + # Beryciformes
  geom_cladelabel(node = 140, label = "", offset = 70, barsize = 1) + # Scombriformes
  geom_cladelabel(node = 218, label = "", offset = 60, barsize = 1) + # Ateleopodiformes
  geom_cladelabel(node = 134, label = "", offset = 70, barsize = 1) + # Carangiformes
  geom_cladelabel(node = 152, label = "", offset = 70, barsize = 1) + # Perciformes
  geom_cladelabel(node = 127, label = "", offset = 60, barsize = 1) # Pleuronectiformes

o2_tree

svg(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "oxygen_tree_2.svg"), width = 10, height = 10)
o2_tree
dev.off()

# Mean O2 ----------------------------------------------------------------------

mean_ln_O2 <- aggregate(new_fish_data$normalised_ln_O2, by = list(sciname = new_fish_data$sciname), FUN = mean)
colnames(mean_ln_O2)[2] <- "normalised_ln_O2"
glimpse(mean_ln_O2)
glimpse(new_fish_tree)

rownames(mean_ln_O2) <- mean_ln_O2$sciname

# Tree plot --------------------------------------------------------------------

o2_tree <- ggtree(new_fish_tree, layout = "fan", lwd = 0.5) + 
  xlim(0, 300) +
  geom_fruit(data = mean_ln_O2,
             geom = geom_bar,
             pwidth = 0.5,
             mapping = aes(y = sciname,
                           x = normalised_ln_O2,
                           fill = normalised_ln_O2),
             orientation = "y",
             stat = "identity") +
  scale_fill_viridis(option = "mako", 
                     name = bquote("Normalised mass-specific \noxygen consumption"),
                     begin = 0.2,
                     end = 0.9) +
  geom_cladelabel(node = 180, label = "", offset = 40, barsize = 1) + # Gadiformes
  geom_cladelabel(node = 178, label = "", offset = 40, barsize = 1) + # Beryciformes
  geom_cladelabel(node = 142, label = "", offset = 40, barsize = 1) + # Scombriformes
  geom_cladelabel(node = 211, label = "", offset = 40, barsize = 1) + # Myctophiformes
  geom_cladelabel(node = 128, label = "", offset = 40, barsize = 1) # Pleuronectiformes

o2_tree
# Get legends 

C_resp_leg <- get_legend(
  C_resp_tree + theme(legend.title = element_blank())
)

C_resp_leg

o2_leg <- get_legend(
  o2_tree + theme(legend.title = element_blank())
)

o2_leg

comb_leg <- plot_grid(
  C_resp_leg,
  o2_leg,
  align = "v",
  ncol = 1
)

comb_leg

# Combine ----------------------------------------------------------------------

comb_tree <- plot_grid(
  C_resp_tree + theme(legend.position = "none",
                      plot.margin = unit(c(0, 0, 0, 0), "cm")),
  o2_tree + theme(legend.position = "none",
                  plot.margin = unit(c(0, 0, 0, 0), "cm")),
  labels = c("A", "B"),
  ncol = 1,
  hjust = -1,
  rel_heights = c(0.7, 0.7)
)

comb_tree 

full_plot <- plot_grid(
  comb_tree,
  comb_leg,
  align = "v",
  ncol = 2,
  hjust = -1,
  rel_widths = c(8, 1.5),
  rel_heights = c(8, 4)
)

full_plot

svg(here("Outputs",
         "08_Plots",
         "Ch_6_Scaling",
         "tree_plot.svg"), width = 5, height = 10)
full_plot
dev.off()

