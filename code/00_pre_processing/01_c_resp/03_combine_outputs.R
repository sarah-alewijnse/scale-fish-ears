#### Combine Posteriors ####

# Combines posteriors from all individuals into a single files

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(MixSIAR) # Gibbs sampler

# Set working directory --------------------------------------------------------

setwd(here("outputs",
           "c_resp",
           "summary_stats"))

# Create table -----------------------------------------------------------------

# Read all csvs in folder
read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

# Make into table
tbl_with_sources <-
  list.files(pattern = "*.csv", 
             full.names = T) %>% 
  map_df(~read_plus(.))

# Remove last column
C_resp <- select(tbl_with_sources, - filename)

# Write into file
write.csv(C_resp, "all_c_resp.csv", row.names = FALSE)
