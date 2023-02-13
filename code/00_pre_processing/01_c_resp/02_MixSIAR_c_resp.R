#### C_resp estimations with MixSIAR ####

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(MixSIAR)
library(here)
library(tidyr)

# Load modified functions ------------------------------------------------------

source(here("code",
            "00_pre_processing",
            "01_c_resp",
            "01_modified_MixSIAR_funcs.R"))

# Load data --------------------------------------------------------------------

component <- read.csv(here("data",
                           "alewijnse_master_data.csv"))

# Create C_resp function -------------------------------------------------------

C.resp <- function(data_sheet, ID_number) {
  
  # Create mixture
  data_f <- dplyr::filter(data_sheet, specimen_id == ID_number)
  
  mixture <- data.frame(ID = data_f$specimen_id,
                        d13C = data_f$oto_d13c_mean)
  
  mix <- load_mix_data_mod(mixture,
                           iso_names = "d13C",
                           factors = "ID",
                           fac_random = FALSE,
                           fac_nested = FALSE,
                           cont_effects = NULL)
  
  # Create source data
  sources_means <- dplyr::select(data_f, specimen_id, dic_d13c_mean, diet_d13c_mean)
  colnames(sources_means) <- c("ID", "DIC", "diet")
  sources_means <- tidyr::pivot_longer(sources_means, 2:3)
  
  sources_sds <- dplyr::select(data_f, specimen_id, dic_d13c_sd, diet_d13c_sd)
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
  
  model_filename <- here("outputs",
                         "c_resp",
                         "c_resp_jags_model.txt")
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
  
  summary_stat_1 <- data.frame(My_Number = data_f$specimen_id,
                               Species = data_f$species,
                               mean_C_resp = mean_1,
                               sd_C_resp = sd_1,
                               se_C_resp = se_1)
  
  write.csv(summary_stat_1, paste0(here("outputs",
                                       "c_resp",
                                       "summary_stats",
                                       "summary_stats_"), ID_number, ".csv"), row.names = F)
  
  # Get diagnostics
  
  output_options <- list(summary_save = FALSE,
                         summary_name = "sum_stat",
                         sup_post = TRUE,
                         plot_post_save_pdf = FALSE,
                         plot_post_name = "post",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = TRUE,
                         plot_pairs_name = paste0(here("outputs",
                                                      "c_resp",
                                                      "pairs_plot_"), ID_number),
                         sup_xy = FALSE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = paste0(here("outputs",
                                                   "c_resp",
                                                   "traceplots_"), ID_number),
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = paste0(here("outputs",
                                                "c_resp",
                                                "diagnostics_"), ID_number),
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

for(i in 1:nrow(component)){
  try(C_resp <- C.resp(component, component$specimen_id[i])) # Add try so it keeps going in case of errors
}
