#### Modified MixSIAR functions ####

# These MixSIAR functions have been modified to enable them to run in a loop

# Load package -----------------------------------------------------------------

library(MixSIAR) # Bayesian mixing model for stable isotopes

# Edit MIXSIAR functions to allow for Looping ----------------------------------

# load_mix_data function =======================================================

load_mix_data_mod <- function (filename, iso_names, factors, fac_random, fac_nested, 
                               cont_effects) 
{
  X <- filename
  n.iso <- length(iso_names)
  if (n.iso == 0) {
    stop(paste("*** Error: No isotopes/tracers selected. Please select 1 or more\n        isotopes/tracers, and then load your consumer/mixture data again. ***", 
               sep = ""))
  }
  if (length(fac_random) != length(factors)) {
    stop(paste("*** Error: You have specified factors to include without saying\n        if they are random or fixed effects (length of fac_random should\n        match length of factors). Please check your load_mix_data line and try again. ***", 
               sep = ""))
  }
  if (length(factors) == 2 && length(fac_nested) != 2) {
    stop(paste("*** Error: You have specified factors to include without saying\n        if they are nested or independent (length of fac_nested should\n        match length of factors). Please check your load_mix_data line and try again. ***", 
               sep = ""))
  }
  if (length(factors) == 2) {
    if (!is.na(fac_nested[1])) {
      if (fac_nested[1] == TRUE && fac_nested[2] == TRUE) {
        stop(paste("*** Error: Both factors cannot be nested within each other. Please check\n              the fac_nested argument in your load_mix_data line and try again. ***", 
                   sep = ""))
      }
    }
  }
  n.effects <- length(factors)
  n.re <- sum(fac_random)
  n.fe <- n.effects - n.re
  fere <- ifelse(n.effects == 2 & n.re < 2, TRUE, FALSE)
  if (n.effects == 1) 
    fac_nested <- FALSE
  if (n.effects > 2) {
    stop(paste("*** Error: More than 2 random/fixed effects selected (MixSIAR can only\n        currently handle 0, 1, or 2 random/fixed effects). Please choose 0, 1,\n        or 2 random/fixed effects and then load your consumer/mixture data again. ***", 
               sep = ""))
  }
  n.ce <- length(cont_effects)
  if (n.ce > 1) {
    stop(paste("*** Error: More than 1 continuous effect selected (MixSIAR can only\n        currently handle 0 or 1 continuous effects). Please choose 0 or 1\n        continuous effects and then load your consumer/mixture data again. ***", 
               sep = ""))
  }
  if (sum(is.na(match(iso_names, colnames(X)))) > 0) {
    stop(paste("*** Error: Your 'iso_names' do not match column names in your\n        mixture data file (case sensitive). Please check your mix .csv data \n        file and load_mix_data line, then try again. ***", 
               sep = ""))
  }
  if (sum(is.na(match(factors, colnames(X)))) > 0) {
    stop(paste("*** Error: Your 'factors' do not match column names in your\n        mixture data file (case sensitive). Please check your mix .csv data\n        file and load_mix_data line, then try again. ***", 
               sep = ""))
  }
  if (sum(is.na(match(cont_effects, colnames(X)))) > 0) {
    stop(paste("*** Error: Your 'cont_effects' do not match column names in your\n        mixture data file (case sensitive). Please check your mix .csv data \n        file and load_mix_data line, then try again. ***", 
               sep = ""))
  }
  N <- dim(X)[1]
  X_iso_cols <- match(iso_names, colnames(X))
  X_iso <- as.matrix(X[, X_iso_cols[]])
  MU_names <- paste("Mean", iso_names, sep = "")
  SIG_names <- paste("SD", iso_names, sep = "")
  FAC <- replicate(n.effects, NULL)
  if (n.effects > 0) {
    for (i in 1:n.effects) {
      re <- fac_random[i]
      fac_values <- X[, factors[i]]
      fac_name <- factors[i]
      fac_levels <- length(unique(fac_values))
      if (is.numeric(fac_values)) {
        fac_labels <- paste(rep(factors[i], fac_levels), 
                            levels(factor(fac_values)), sep = " ")
      }
      else {
        fac_labels <- levels(factor(fac_values))
      }
      fac_values <- as.numeric(factor(fac_values))
      FAC[[i]] <- list(values = fac_values, levels = fac_levels, 
                       labels = fac_labels, lookup = NULL, re = re, 
                       name = fac_name)
    }
    if (n.re == 2 & !is.na(fac_nested[1])) {
      if (n.re == 2 & fac_nested[2]) {
        for (lev in 1:FAC[[2]]$levels) {
          FAC[[2]]$lookup[lev] <- FAC[[1]]$values[which(FAC[[2]]$values == 
                                                          lev)][1]
        }
      }
      if (n.re == 2 & fac_nested[1]) {
        for (lev in 1:FAC[[1]]$levels) {
          FAC[[1]]$lookup[lev] <- FAC[[2]]$values[which(FAC[[1]]$values == 
                                                          lev)][1]
        }
      }
    }
    if (n.fe == 1 & n.re == 1 & fac_random[1]) {
      tmp <- FAC[[1]]
      FAC[[1]] <- FAC[[2]]
      FAC[[2]] <- tmp
      factors <- rev(factors)
      fac_random <- rev(fac_random)
      fac_nested <- rev(fac_nested)
    }
  }
  CE_orig <- replicate(n.ce, NULL)
  CE <- replicate(n.ce, NULL)
  CE_center <- rep(NA, n.ce)
  CE_scale <- rep(NA, n.ce)
  if (n.ce > 0) {
    for (i in 1:n.ce) {
      CE_orig[[i]] <- X[, cont_effects[i]]
      CE[[i]] <- scale(X[, cont_effects[i]], center = TRUE, 
                       scale = TRUE)
      CE_center[i] <- attributes(CE[[i]])$"scaled:center"
      CE_scale[i] <- attributes(CE[[i]])$"scaled:scale"
    }
  }
  return(list(data = X, data_iso = X_iso, n.iso = n.iso, n.re = n.re, 
              n.ce = n.ce, FAC = FAC, CE = CE, CE_orig = CE_orig, CE_center = CE_center, 
              CE_scale = CE_scale, cont_effects = cont_effects, MU_names = MU_names, 
              SIG_names = SIG_names, iso_names = iso_names, N = N, 
              n.fe = n.fe, n.effects = n.effects, factors = factors, 
              fac_random = fac_random, fac_nested = fac_nested, fere = fere))
}

## Edit load_source_data function ==============================================

load_source_data_mod <-  function (filename, source_factors = NULL, conc_dep, data_type, 
                                   mix) 
{
  SOURCE <- filename
  source.fac <- length(source_factors)
  if (source.fac > 1) {
    stop(paste("More than one source factor.\n    MixSIAR can only fit source data by up to ONE factor.\n    Please specify 0 or 1 source factor and try again.", 
               sep = ""))
  }
  if (sum(is.na(match(source_factors, colnames(SOURCE)))) > 
      0) {
    stop(paste("Your 'source_factors' do not match column names in your\n        source data file (case sensitive). Please check your source .csv data\n        file and load_source_data line, then try again. ***", 
               sep = ""))
  }
  test_fac <- match(source_factors, mix$factors)
  if (source.fac == 1 && (length(test_fac) == 0 || is.na(test_fac))) {
    stop(paste("Source factor not in mix$factors.\n    You cannot model a source random effect that is not included\n    as a random/fixed effect for the mixture/consumer. Either\n     1) remove the source factor (reload source data), or\n     2) include the random/fixed effect in the mixture (reload mix data).\n    Could be a mismatch between column headings in the mix and source data files.", 
               sep = ""))
  }
  if (source.fac == 0) 
    by_factor <- NA
  else by_factor <- match(source_factors, mix$factors)
  src <- as.factor(SOURCE[, 1])
  source_names <- levels(src)
  n.sources <- length(source_names)
  levels(src) <- 1:n.sources
  SOURCE[, 1] <- as.numeric(src)
  source_factor_cols <- match(source_factors, colnames(SOURCE))
  S_factor_levels <- rep(0, length(source_factor_cols))
  if (!is.na(by_factor)) {
    SOURCE <- SOURCE[order(SOURCE[, 1], SOURCE[, source_factor_cols]), 
                     ]
  }
  else {
    SOURCE <- SOURCE[order(SOURCE[, 1]), ]
  }
  if (conc_dep) {
    CONC_names <- paste("Conc", mix$iso_names, sep = "")
    if (sum(is.na(match(CONC_names, colnames(SOURCE)))) > 
        0) {
      stop(paste("Concentration dependence column names mislabeled.\n    Should be 'Conc' + iso_names, e.g. 'Concd13C' if iso_names = 'd13C'.\n    Please ensure Conc headings in source data file match iso_names\n    in mix data file and try again. Alternatively, you may have set conc_dep=T\n    by mistake.", 
                 sep = ""))
    }
    CONC_iso_cols <- match(CONC_names, colnames(SOURCE))
    conc <- do.call(rbind, lapply(split(SOURCE[, CONC_iso_cols], 
                                        list(SOURCE[, 1])), colMeans))
  }
  else conc <- NULL
  col_sd <- function(x) {
    sds <- apply(x, 2, sd)
    return(sds)
  }
  if (data_type == "raw") {
    if (sum(is.na(match(mix$iso_names, colnames(SOURCE)))) > 
        0) {
      stop(paste("With raw source data, the iso_names in mix data\n        file must be in the column headings of source data file. Please check\n        your source and mix .csv data files and try again. ***", 
                 sep = ""))
    }
    S_iso_cols <- match(mix$iso_names, colnames(SOURCE))
    if (!is.na(by_factor)) {
      for (fac in 1:length(source_factor_cols)) {
        S_factor_levels[fac] <- length(levels(SOURCE[, 
                                                     source_factor_cols[fac]]))
      }
      list.sources.bylev <- vector("list", S_factor_levels)
      test.sources.identical <- rep(TRUE, S_factor_levels)
      S_factor_levels_names <- levels(as.factor(SOURCE[, 
                                                       source_factor_cols[fac]]))
      for (lev in 1:S_factor_levels) {
        list.sources.bylev[[lev]] <- SOURCE[SOURCE[, 
                                                   source_factors] == S_factor_levels_names[lev], 
                                            1]
      }
      if (length(unique(list.sources.bylev)) != 1) {
        stop(paste("Sources for each level of *", source_factors, 
                   "* do not match.\n        If you have different sources (or # of sources) for levels of *", 
                   source_factors, "*,\n        you must fit separate MixSIAR models for each level.", 
                   sep = ""))
      }
      if (mix$n.iso > 1) 
        mu <- do.call(rbind, lapply(split(SOURCE[, S_iso_cols], 
                                          list(SOURCE[, source_factor_cols[1]], SOURCE[, 
                                                                                       1])), colMeans))
      if (mix$n.iso == 1) 
        mu <- do.call(rbind, lapply(split(SOURCE[, S_iso_cols], 
                                          list(SOURCE[, source_factor_cols[1]], SOURCE[, 
                                                                                       1])), mean))
      S_MU <- cbind(mu, rep(1:S_factor_levels, n.sources))
      if (mix$n.iso > 1) 
        sig <- do.call(rbind, lapply(split(SOURCE[, S_iso_cols], 
                                           list(SOURCE[, source_factor_cols[1]], SOURCE[, 
                                                                                        1])), col_sd))
      if (mix$n.iso == 1) 
        sig <- do.call(rbind, lapply(split(SOURCE[, S_iso_cols], 
                                           list(SOURCE[, source_factor_cols[1]], SOURCE[, 
                                                                                        1])), sd))
      S_SIG <- cbind(sig, rep(1:S_factor_levels, n.sources))
      colnames(S_MU) <- c(mix$iso_names, source_factors)
      S_MU_factor_col <- match(source_factors, colnames(S_MU))
      S_factor1 <- S_MU[, S_MU_factor_col]
    }
    if (is.na(by_factor)) {
      if (mix$n.iso > 1) 
        S_MU <- do.call(rbind, lapply(split(SOURCE[, 
                                                   S_iso_cols], list(SOURCE[, 1])), colMeans))
      if (mix$n.iso == 1) 
        S_MU <- do.call(rbind, lapply(split(SOURCE[, 
                                                   S_iso_cols], list(SOURCE[, 1])), mean))
      if (mix$n.iso > 1) 
        S_SIG <- do.call(rbind, lapply(split(SOURCE[, 
                                                    S_iso_cols], list(SOURCE[, 1])), col_sd))
      if (mix$n.iso == 1) 
        S_SIG <- do.call(rbind, lapply(split(SOURCE[, 
                                                    S_iso_cols], list(SOURCE[, 1])), sd))
      colnames(S_MU) <- c(mix$iso_names)
      S_factor1 <- NULL
      S_factor_levels <- NULL
    }
    n.rep <- array(0, dim = c(n.sources, S_factor_levels))
    if (!is.na(by_factor)) {
      for (src in 1:n.sources) {
        for (f1 in 1:S_factor_levels) {
          n.rep[src, f1] <- table(SOURCE[which(SOURCE[, 
                                                      1] == src), source_factor_cols])[f1]
        }
      }
    }
    else {
      for (src in 1:n.sources) {
        n.rep[src] <- length(which(SOURCE[, 1] == src))
      }
    }
    max.rep <- max(n.rep)
    SOURCE_array <- array(NA, dim = c(n.sources, mix$n.iso, 
                                      S_factor_levels, max.rep))
    count <- 1
    if (!is.na(by_factor)) {
      for (src in 1:n.sources) {
        for (f1 in 1:S_factor_levels) {
          for (r in 1:n.rep[src, f1]) {
            for (iso in 1:mix$n.iso) {
              SOURCE_array[src, iso, f1, r] <- SOURCE[count, 
                                                      mix$iso_names[iso]]
            }
            count <- count + 1
          }
        }
      }
    }
    else {
      for (src in 1:n.sources) {
        for (r in 1:n.rep[src]) {
          for (iso in 1:mix$n.iso) {
            SOURCE_array[src, iso, r] <- SOURCE[count, 
                                                mix$iso_names[iso]]
          }
          count <- count + 1
        }
      }
    }
    MU_array <- NULL
    SIG2_array <- NULL
    n_array <- NULL
  }
  if (data_type == "means") {
    if (sum(is.na(match(mix$MU_names, colnames(SOURCE)))) > 
        0) {
      stop(paste("Source mean column names mislabeled.\n    Should be 'Mean' + iso_names from mix data file, e.g. 'Meand13C' if\n    mix$iso_names = 'd13C'. Please ensure headings in source data file match\n    this format and try again. Alternatively, if you have raw source data,\n    you should set data_type='raw'.", 
                 sep = ""))
    }
    if (sum(is.na(match(mix$SIG_names, colnames(SOURCE)))) > 
        0) {
      stop(paste("Source SD column names mislabeled.\n    Should be 'SD' + iso_names from mix data file, e.g. 'SDd13C' if\n    mix$iso_names = 'd13C'. Please ensure headings in source data file match\n    this format and try again. Alternatively, if you have raw source data,\n    you should set data_type='raw'.", 
                 sep = ""))
    }
    S_MU_iso_cols <- match(mix$MU_names, colnames(SOURCE))
    S_SIG_iso_cols <- match(mix$SIG_names, colnames(SOURCE))
    sample_size_col <- match("n", colnames(SOURCE))
    if (is.na(sample_size_col)) {
      stop(paste("Source sample sizes missing or entered incorrectly.\n        Check your sources .csv data file to be sure you have a column\n        titled \"n\" with the sample sizes for each source isotope estimate.", 
                 sep = ""))
    }
    S_sample_size <- SOURCE[, sample_size_col]
    if (is.na(by_factor)) {
      S_MU <- as.matrix(SOURCE[, c(S_MU_iso_cols[], source_factor_cols)])
      S_SIG <- as.matrix(SOURCE[, c(S_SIG_iso_cols[], source_factor_cols)])
      MU_array <- S_MU
      SIG2_array <- S_SIG * S_SIG
      n_array <- S_sample_size
      S_factor1 <- NULL
      S_factor_levels <- NULL
    }
    else {
      for (fac in 1:length(source_factor_cols)) {
        S_factor_levels[fac] <- length(levels(as.factor(SOURCE[, 
                                                               source_factor_cols[fac]])))
      }
      list.sources.bylev <- vector("list", S_factor_levels)
      test.sources.identical <- rep(TRUE, S_factor_levels)
      S_factor_levels_names <- levels(as.factor(SOURCE[, 
                                                       source_factor_cols[fac]]))
      for (lev in 1:S_factor_levels) {
        list.sources.bylev[[lev]] <- SOURCE[SOURCE[, 
                                                   source_factors] == S_factor_levels_names[lev], 
                                            1]
      }
      if (length(unique(list.sources.bylev)) != 1) {
        stop(paste("Sources for each level of *", source_factors, 
                   "* do not match.\n        If you have different sources (or # of sources) for levels of *", 
                   source_factors, "*,\n        you must fit separate MixSIAR models for each level.", 
                   sep = ""))
      }
      if (length(S_sample_size) != (n.sources * S_factor_levels)) {
        stop(paste("Source sample sizes missing or entered incorrectly.\n        Check your source_means.csv data file to be sure you have a column\n        titled \"n\" with the sample sizes for each source isotope estimate.", 
                   sep = ""))
      }
      S_factor1 <- factor(SOURCE[, source_factor_cols])
      SOURCE[, source_factor_cols] <- as.numeric(factor(SOURCE[, 
                                                               source_factor_cols]))
      S_MU <- as.matrix(SOURCE[, c(S_MU_iso_cols[], source_factor_cols)])
      S_SIG <- as.matrix(SOURCE[, c(S_SIG_iso_cols[], source_factor_cols)])
      MU_array <- array(NA, dim = c(n.sources, mix$n.iso, 
                                    S_factor_levels))
      SIG2_array <- array(NA, dim = c(n.sources, mix$n.iso, 
                                      S_factor_levels))
      n_array <- array(NA, dim = c(n.sources, S_factor_levels))
      count <- 1
      for (src in 1:n.sources) {
        for (f1 in 1:S_factor_levels) {
          for (iso in 1:mix$n.iso) {
            MU_array[src, iso, f1] <- S_MU[count, iso]
            SIG2_array[src, iso, f1] <- S_SIG[count, 
                                              iso] * S_SIG[count, iso]
          }
          n_array[src, f1] <- S_sample_size[count]
          count <- count + 1
        }
      }
    }
    SOURCE_array <- NULL
    n.rep <- NULL
  }
  if (length(which(S_SIG == 0)) > 0) {
    stop(paste("You have at least one source SD = 0.\n    Check your source data file to be sure each SD entry is non-zero.", 
               sep = ""))
  }
  return(list(n.sources = n.sources, source_names = source_names, 
              S_MU = S_MU, S_SIG = S_SIG, S_factor1 = S_factor1, S_factor_levels = S_factor_levels, 
              conc = conc, MU_array = MU_array, SIG2_array = SIG2_array, 
              n_array = n_array, SOURCE_array = SOURCE_array, n.rep = n.rep, 
              by_factor = by_factor, data_type = data_type, conc_dep = conc_dep))
}

## Edit load_discr_data function ===============================================

load_discr_data_mod <- function (filename, mix) 
{
  DISCR <- filename
  row.names(DISCR) <- DISCR[, 1]
  DISCR <- as.matrix(DISCR[-1])
  DISCR <- DISCR[order(rownames(DISCR)), ]
  if (sum(is.na(match(mix$MU_names, colnames(DISCR)))) > 0) {
    stop(paste("*** Error: Discrimination mean column names mislabeled.\n    Should be 'Mean' + iso_names from mix data file, e.g. 'Meand13C' if\n    mix$iso_names = 'd13C'. Please ensure headings in discr data file match\n    this format and try again.", 
               sep = ""))
  }
  if (sum(is.na(match(mix$SIG_names, colnames(DISCR)))) > 0) {
    stop(paste("*** Error: Discrimination SD column names mislabeled.\n    Should be 'SD' + iso_names from mix data file, e.g. 'SDd13C' if\n    mix$iso_names = 'd13C'. Please ensure headings in discr data file match\n    this format and try again.", 
               sep = ""))
  }
  discr_mu_cols <- match(mix$MU_names, colnames(DISCR))
  discr_sig_cols <- match(mix$SIG_names, colnames(DISCR))
  discr_mu <- as.matrix(DISCR[, discr_mu_cols])
  discr_sig2 <- as.matrix(DISCR[, discr_sig_cols] * DISCR[, 
                                                          discr_sig_cols])
  return(list(mu = discr_mu, sig2 = discr_sig2))
}
