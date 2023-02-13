#### Combine Posteriors ####

# Combines posteriors from all individuals into a single files

# Load packages

library(tidyverse)
library(MixSIAR) # Gibbs sampler

#### C_resp ####

# Set working directory

setwd("~/PhD/GitHub/fish-ears/Outputs/05_C_resp_Calculations/Summary_stats")

# Read all .csv files in folder

read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

tbl_with_sources <-
  list.files(pattern = "*.csv", 
             full.names = T) %>% 
  map_df(~read_plus(.))

# Remove last column

C_resp <- select(tbl_with_sources, - filename)

# Write into file

write.csv(C_resp, "All_C_resp.csv", row.names = FALSE)
