#### Conversion to O2 consumption ####

library(tidyverse)
library(here)

# Load data --------------------------------------------------------------------

C_resp <- read.csv(here::here("Outputs",
                              "07_C_resp_Models",
                              "Edited_Data",
                              "master_data_no_outliers.csv"))

# Convert to mass-specific O2 consumption --------------------------------------

# Create function

mg_O2_kg <- function(C_resp, C, k){
  a <- C_resp/C
  mg_kg <- -(log(1-a)/k)
  print(mg_kg)
}

# Set parameters

k <- 0.004 # Decay constant - Martino et al. 2020
C <- max(C_resp$mean_C_resp) + 0.01 # Upper bound

# Convert

C_resp$O2_mg_kg <- mg_O2_kg(C_resp$mean_C_resp, C, k)

# Explore outcomes -------------------------------------------------------------

# Histogram

ggplot(data = C_resp, aes(x = O2_mg_kg)) +
  geom_histogram()

# Plot with body mass

ggplot(data = C_resp, aes(x = log10(Body_Mass), y = O2_mg_kg)) +
  geom_point()

# Plot with temperature

ggplot(data = C_resp, aes(x = mean_Temp, y = O2_mg_kg)) +
  geom_point()

# Convert to whole organism oxygen consumption ---------------------------------

C_resp$O2_mg_ind <- C_resp$O2_mg_kg * (C_resp$Body_Mass/1000)

# Explore outcomes -------------------------------------------------------------

# Histogram

ggplot(data = C_resp, aes(x = log10(O2_mg_ind))) +
  geom_histogram()

# Plot with body mass

ggplot(data = C_resp, aes(x = log10(Body_Mass), y = log10(O2_mg_ind),
                          col = mean_Temp)) +
  geom_point()

# Plot with temperature

ggplot(data = C_resp, aes(x = mean_Temp, y = log10(O2_mg_ind))) +
  geom_point()

# Save ------------------------------------------------------------------------

write.csv(C_resp, here("Outputs",
                       "07_C_resp_Models",
                       "Edited_Data",
                       "Oxygen_Consumption.csv"),
          row.names = F)
