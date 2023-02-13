#### C_resp estimations with MixSIAR ####

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(MixSIAR)
library(here)
library(tidyr)

# Load modified functions ------------------------------------------------------

load(here("Functions",
          "load_discr_data_mod.Rdata"))
load(here("Functions",
          "load_mix_data_mod.Rdata"))
load(here("Functions",
          "load_source_data_mod.Rdata"))

# Load data --------------------------------------------------------------------

component <- read.csv(here("Outputs",
                           "04_C_resp_Components",
                           "Compiled",
                           "All_Data.csv"))

# Create C_resp function -------------------------------------------------------

C.resp <- function(data_sheet, ID_number) {
  
  # Create mixture
  
  data_f <- dplyr::filter(data_sheet, ID == ID_number)
  
  mixture <- data.frame(ID = data_f$ID,
                        d13C = data_f$oto_d13C)
  
  mix <- load_mix_data_mod(mixture,
                           iso_names = "d13C",
                           factors = "ID",
                           fac_random = FALSE,
                           fac_nested = FALSE,
                           cont_effects = NULL)
  
  # Create source data
  
  sources_means <- dplyr::select(data_f, ID, DIC_mean, diet_mean)
  colnames(sources_means) <- c("ID", "DIC", "diet")
  sources_means <- tidyr::pivot_longer(sources_means, 2:3)
  
  sources_sds <- dplyr::select(data_f, ID, DIC_SD, diet_SD)
  colnames(sources_sds) <- c("ID", "DIC", "diet")
  sources_sds <- tidyr::pivot_longer(sources_sds, 2:3)
  
  sources <- data.frame(Sources = sources_means$name,
                        Meand13C = sources_means$value,
                        SDd13C = sources_sds$value,
                        ID = sources_means$ID,
                        n = 10000)
  
  source <- load_source_data_mod(sources,
                                 source_factors = NULL,
                                 conc_dep = FALSE,
                                 data_type = "means",
                                 mix)
  
  # Create discrimination data
  
  disc <- data.frame(Source = c("DIC", "diet"),
                     Meand13C = c(0, 0),
                     SDd13C = c(0, 0))
  
  discr <- load_discr_data_mod(disc, mix)
  
  # Write JAGS model
  
  model_filename <- "Outputs/05_C_resp_Calculations/C_resp_Model.txt"
  resid_err <- FALSE
  process_err <- TRUE
  MixSIAR::write_JAGS_model(model_filename, resid_err, process_err, mix, source)
  
  # Run the model
  
  test_mod <- run_model(run = "normal", mix, source, discr, model_filename,
                        alpha.prior = 1, resid_err, process_err)
  
  # Get output 
  
  R2jags::attach.jags(test_mod)
  post_full_1 <- p.fac1[,1,2]
  post_full_1 <- as.data.frame(post_full_1)
  colnames(post_full_1) <- "C_resp"
  
  # Get posterior summary
  
  mean_1 <- mean(post_full_1$C_resp, na.rm = TRUE)
  se <- function(x) sd(x)/sqrt(length(x))
  sd_1 <- sd(post_full_1$C_resp, na.rm = TRUE)
  se_1 <- se(post_full_1$C_resp)
  
  summary_stat_1 <- data.frame(My_Number = data_f$ID,
                               Species = data_f$Species,
                               mean_C_resp = mean_1,
                               sd_C_resp = sd_1,
                               se_C_resp = se_1)
  
  write.csv(summary_stat_1, paste("Outputs/05_C_resp_Calculations/Summary_stats/", ID_number, ".csv"), row.names = F)
  
  # Get diagnostics
  
  output_options <- list(summary_save = FALSE,
                         summary_name = "sum_stat",
                         sup_post = TRUE,
                         plot_post_save_pdf = FALSE,
                         plot_post_name = "post",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = TRUE,
                         plot_pairs_name = paste("Outputs/05_C_resp_Calculations/Pairs_plots/", ID_number),
                         sup_xy = FALSE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = paste("Outputs/05_C_resp_Calculations/Traceplots/", ID_number),
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = paste("Outputs/05_C_resp_Calculations/Diagnostics/", ID_number),
                         indiv_effect = FALSE,
                         plot_post_save_png = FALSE,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)
  
  output_JAGS(test_mod, mix, source, output_options)
  
  return(summary_stat_1)

}

# Test function ----------------------------------------------------------------

C.resp(data_sheet = component,
       ID_number = "AAG_2")

#### Run for all ####

C_resp_summary <- data.frame()

for(i in 1200:nrow(component)){
  C_resp <- C.resp(component, component$ID[i])
  C_resp_summary <- rbind(C_resp, C_resp_summary)
}

# Write into file

write.csv(C_resp_summary, "Outputs/05_C_resp_Calculations/C_resp_summary.csv", row.names = F)
