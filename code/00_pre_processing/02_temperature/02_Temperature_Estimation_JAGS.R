#### Temperature estimation ####

# Estimates experienced temperature based on oxygen isotopes in RJAGS

# Load required packages -------------------------------------------------------

library(tidyverse)
library(rjags) # Gibbs sampler
library(coda) # Summarises MCMC outputs
library(here)

# Load data --------------------------------------------------------------------

component <- read.csv(here("Outputs",
                           "04_C_resp_Components",
                           "Compiled",
                           "All_Data.csv"))

# Calculate difference between d18O_oto and d18O_water -------------------------

component$d18O_diff <- component$oto_d18O - component$d18O_SW_mean
glimpse(component)

# Filter out those without d18O

component <- filter(component, !is.na(d18O_diff))

glimpse(component)

# Create temperature function --------------------------------------------------

# Enables you to loop the temperature function over the whole dataset

Temp_est <- function(data_set, Data_ID){
  fish_1 <- dplyr::filter(data_set, ID == Data_ID)
  # Set priors/data
  
  iso_list <- list(
    iso = fish_1$d18O_diff, # Difference between d18O_oto and d18O_water
    sigma = 1/(fish_1$oto_d18O_SD)^2,
    
    # Parameters from Hoie et al. 2004
    
    a_obs = 3.90,
    a_var = 1/(0.24^2),
    b_obs = -0.20,
    b_var = 1/(0.019^2),
    N = 1
  )
  
  inits <- list(Temp = 0.0)
  
  # Write the JAGS model
  
  cat("model
      {
      for (i in 1:N){
      mu[i] <- a_est + Temp * b_est
      iso[i] ~ dnorm(mu[i], sigma)
      }
      a_est ~ dnorm(a_obs, a_var)
      b_est ~ dnorm(b_obs, b_var)
      Temp ~ dunif(-5, 35)
      }", file = "Outputs/06_Temp_Calculations/Temp_Jags.txt")

  # Run JAGS model
  
  jags_mod <- jags.model(file = "Outputs/06_Temp_Calculations/Temp_JAGS.txt", data = iso_list, inits = inits, n.chains = 3, n.adapt = 50000)
  
  output <- coda.samples(jags_mod,
                         c("Temp"),
                         n.iter = 100000,
                         thin = 50)
  
  ## Get the outputs
  
  # View summary
  
  summary(output)
  
  # Get posterior
  
  post_1 <- as.data.frame(output[[1]])
  post_2 <- as.data.frame(output[[2]])
  post_3 <- as.data.frame(output[[3]])
  
  post_full <- rbind(post_1, post_2, post_3)
  
  # Get posterior summary
  
  mean <- mean(post_full$Temp, na.rm = TRUE)
  sd <- sd(post_full$Temp, na.rm = TRUE)
  se <- sd/sqrt(length(post_full$Temp))
  
  summary_stat <- data.frame(Data_ID = Data_ID,
                             mean_Temp = mean,
                             sd_Temp = sd,
                             SE_Temp = se)
  summary_stat
  
  ## Get diagnostics
  
  # ESS
  
  samp_size <- as.data.frame(effectiveSize(output))
  samp_size
  
  # Gelman-Rubin diagnostic
  
  r_hat <- gelman.diag(output)$psrf
  r_hat
  
  # Geweke's diagnostic
  
  geweke <- geweke.diag(output)
  
  geweke_1 <- geweke[[1]]
  geweke_2 <- geweke[[2]]
  geweke_3 <- geweke[[3]]
  
  geweke_1_z <- geweke_1$z
  geweke_2_z <- geweke_2$z
  geweke_3_z <- geweke_3$z
  
  geweke_z_scores <- c(geweke_1_z, geweke_2_z, geweke_3_z)
  geweke_z_scores
  # Want these to all be below 1.96 (or -1.96)
  
  capture.output(c("r_hat", r_hat, 
                   "geweke", geweke_z_scores, 
                   "ESS", samp_size), file = paste("Outputs/06_Temp_Calculations/Diagnostics/Diagnostics_", Data_ID, ".txt", sep = ""))
  
  # Traceplot
  
  svg(file = paste("Outputs/06_Temp_Calculations/Traceplots/Traceplot_", Data_ID, ".svg", sep = ""), width = 12)
  plot(output)
  dev.off()
  
  return(summary_stat)
  print(summary_stat)
}

# Test with a single individual ------------------------------------------------

Temp_est(component, "AAG_1")

# Loop over whole dataset

Temp_summary <- data.frame()

for(i in 1:nrow(fish)){
  Temp <- with(fish[i,],
       Temp_est(Data_ID))
  Temp_summary <- rbind(Temp_summary, Temp)
}

ggplot(Temp_summary, aes(x = mean_Temp)) +
  geom_histogram()

ggplot(Temp_summary, aes(x = SE_Temp)) +
  geom_histogram()

#### Join to data ####

fish_Temp <- left_join(fish, Temp_summary, by = "Data_ID")

# Filter out those without d18O

fish_Temp_tidy <- filter(fish_Temp, !is.na(oto_d18O))

# Write into file

write.csv(fish_Temp_tidy, "Outputs/06_Temp_Calculations/All_Temp.csv", row.names = F)
