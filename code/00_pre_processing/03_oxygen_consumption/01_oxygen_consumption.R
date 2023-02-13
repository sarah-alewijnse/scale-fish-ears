#### Conversion to O2 consumption ####

library(tidyverse)
library(here)

# Load data --------------------------------------------------------------------

C_resp <- read.csv(here::here("data",
                              "alewijnse_master_data.csv"))

# Convert to mass-specific O2 consumption --------------------------------------

# Create function
mg_O2_kg <- function(C_resp, C, k){
  a <- C_resp/C
  mg_kg <- -(log(1-a)/k)
  print(mg_kg)
}

# Set parameters
k <- 0.004 # Decay constant - Martino et al. 2020
C <- max(C_resp$c_resp_mean) + 0.01 # Upper bound

# Convert
C_resp$o2_mg_kg <- mg_O2_kg(C_resp$c_resp_mean, C, k)

# Explore outcomes -------------------------------------------------------------

# Histogram
ggplot(data = C_resp, aes(x = o2_mg_kg)) +
  geom_histogram()

# Plot with body mass
ggplot(data = C_resp, aes(x = log10(body_mass), y = o2_mg_kg)) +
  geom_point()

# Plot with temperature
ggplot(data = C_resp, aes(x = temp_mean, y = o2_mg_kg)) +
  geom_point()

# Convert to whole organism oxygen consumption ---------------------------------

C_resp$o2_mg_ind <- C_resp$o2_mg_kg * (C_resp$body_mass/1000)

# Explore outcomes -------------------------------------------------------------

# Histogram
ggplot(data = C_resp, aes(x = log10(o2_mg_ind))) +
  geom_histogram()

# Plot with body mass
ggplot(data = C_resp, aes(x = log10(body_mass), y = log10(o2_mg_ind),
                          col = temp_mean)) +
  geom_point()

# Plot with temperature
ggplot(data = C_resp, aes(x = temp_mean, y = log10(o2_mg_ind))) +
  geom_point()

# Save ------------------------------------------------------------------------

write.csv(C_resp, here("outputs",
                       "oxygen_consumption",
                       "alewijnse_master_data_o2_consumption.csv"),
          row.names = F)
