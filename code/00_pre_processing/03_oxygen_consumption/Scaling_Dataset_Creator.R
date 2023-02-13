#### C_resp vs. BMT dataset ####

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

oxy_data <- read.csv(here("Outputs",
                          "07_C_resp_Models",
                          "Edited_Data",
                          "Oxygen_Consumption.csv"))

# Get lnBM

oxy_data$ln_BM <- log(oxy_data$Body_Mass)

# Get arrhenius temperature

boltz <- 8.62*10^-5
oxy_data$mean_Temp_K <- oxy_data$mean_Temp + 273.15
oxy_data$inv_Temp <- 1/(oxy_data$mean_Temp_K * boltz)
oxy_data$ln_inv_temp <- log(oxy_data$inv_Temp)

# Get ln oxygen consumption

oxy_data$ln_O2_ind <- log(oxy_data$O2_mg_ind)