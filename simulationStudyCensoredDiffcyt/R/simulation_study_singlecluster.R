#' Wrapper for running simulations
#'
#' @param reps number of repetitions for each condition
#' @param nr_cores number of cores to use for parallelization
#' @param n_obs_c vector of data sizes
#' @param simDataTypes vector of simulation types, options are "lm", "glm","glmer"
#' @param formulas formulas specifing data, details see  \code{\link{simulate_data}}
#' @param n_levels_fixeff number of levels for each fixed effect covariate
#'  specified in formulas, set to 2 if binary. has to be 2 if run with
#'  'run_with_covariate_dependent_censoring' = TRUE.
#' @param n_levels_raneff number of levels for each random effect covariate
#'  specified in formulas
#' @param betas cofficients in linear term. list of elements, elements are
#'  b0,b1,... each of which is a vector of values.
#' @param censoring_rates vector of approximate censoring rates used. effective
#'  censoring rate is enforced to be plus minus 5 percent of this value.
#' @param imputation_types vector of \code{\link{conditional_multiple_imputation}}
#'  types used for inference, possible options are "cc", "mrl", "km", "rs", "pmm",
#'  "ppd"
#' @param run_with_covariate_dependent_censoring logical, if an additional
#'  alternative censoring mechanismus is used where censoring depends on the
#'  level of a covariate
#' @param mi_reps number of repetitions in \code{\link{conditional_multiple_imputation}}
#' @param error_variance positive double. Variance of additional gaussian
#'  noise to add in the linear sum of the predictors. For linear regression
#'  this is the only error added.
#' @param variance_raneff positive double. The variance of the gaussian distributed
#'  random effect covariates
#' @param transform_fn vector of functions to transform censored covariate or
#'  'identity' (no transformation) or 'boxcox' (box-cox transformation).
#'   default = 'identity'
#' @param verbose verbose
#' @param seed positive integer, random seed
#' @param size_of_random_subset_for_testing number of simulations for test run
#' @importFrom magrittr %>%
#' @export
run_simulations_wrapper_singlecluster <-
  function(reps = 2,
           nr_cores = 1,
           n_obs_c = c(50),
           simDataTypes = c("lm","glm","glmer"),
           formulas = list(lm = list(formula(Y~Surv(X,I) + Z)),
                           glm = list(formula(Y~Surv(X,I) + Z)),
                           glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
           n_levels_fixeff = c(2),
           n_levels_raneff = NULL,
           # signal_to_noise_ratio = c(1),
           betas = list(b0=c(0),b1=c(0.5),b2=c(0.5)),
           censoring_rates = c(0.5,0.6),
           imputation_types = c("mrl","km","rs"),
           run_with_covariate_dependent_censoring = TRUE,
           mi_reps = 10,
           weibull_params = list(shape_x = 0.5,
                                 shape_c1 = 1,
                                 shape_c2 = 0.5,
                                 scale_x = 0.25),
           error_variance = 1,
           variance_raneff = 1,
           transform_fn = "identity",
           verbose = FALSE,
           seed = 123,
           size_of_random_subset_for_testing = NULL) {
    # directory to save tempfiles
    dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
    dir.create(dir_save, showWarnings = FALSE)
    filenames <- c()

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)

    if (run_with_covariate_dependent_censoring){
      cov_dep_cens <- c(FALSE,TRUE)
    } else{
      cov_dep_cens <- c(FALSE)
    }
    if (verbose) cat("start censoring parameter estimation ... ")
    tmp_df_new <- scales_for_censoring(
      censoring_val = censoring_rates,
      log_ratio_val = c(0),
      # seq_to_test = c(seq(0.001,0.99,length.out = 100), seq(1,5,by=0.25)),
      seq_to_test = c(seq(0.001,1,length.out = 100),seq(1,20000,by=10)),
      # seq_to_test = c(seq(0.01,0.99,length.out = 10), seq(1,5,by=0.5)),
      n = 10000,
      shapes = list(shape_x=weibull_params[["shape_x"]],
                    shape_c1=weibull_params[["shape_c1"]],
                    shape_c2=weibull_params[["shape_c2"]]),
      scale_x = weibull_params[["scale_x"]]
      # shapes = list(shape_x = 0.5,
      #               shape_c1 = 1,
      #               shape_c2 = 0.5),
      # scale_x = 0.25
      # shapes = list(shape_x = 1.81854,
      #               shape_c1 = 0.74,
      #               shape_c2 = 0.5),
      # scale_x = 2628
    ) %>%
      dplyr::rename(censoring_rate = "censoring",
                    cov_dep_cens = "log_ratio") %>%
      dplyr::mutate(cov_dep_cens = ifelse(cov_dep_cens == 0, FALSE, TRUE))
    if (verbose) cat("\t done\n")
    print(tmp_df_new)
    # create all simulation conditions to test
    params_to_test <- expand.grid(
      simType = simDataTypes,
      censoring_rate = censoring_rates,
      n_obs = n_obs_c,
      b0 = betas[["b0"]],
      b1 = betas[["b1"]],
      b2 = betas[["b2"]],
      mi_rep = mi_reps,
      cov_dep_cens = cov_dep_cens,
      error_variance = error_variance,
      variance_raneff = variance_raneff,
      transform_fn = transform_fn,
      stringsAsFactors = FALSE
    )
    params_to_test <- suppressWarnings(
      dplyr::inner_join(params_to_test,
                        tmp_df_new,
                        by = c("censoring_rate" = "censoring_rate",
                               "cov_dep_cens" = "cov_dep_cens")
      ) %>%
        dplyr::inner_join(
          tibble::tibble(
            simType=simDataTypes,
            formula = unlist(purrr::map(formulas,~as.character(.x[[1]])[3]))),
          by = c("simType" = "simType")
        )
    )
    params_to_test <- suppressWarnings(
      purrr::reduce(1:(reps + 1), ~ rbind(.x, params_to_test))[-1, ] %>%
        dplyr::mutate(simType = as.character(simType)))
    params_to_test$sim_id <- seq_len(dim(params_to_test)[1])
    # random shuffling to (hopefully) reduce overhead
    params_to_test <- params_to_test[sample(seq_len(dim(params_to_test)[1]),dim(params_to_test)[1], replace = FALSE), ]
    if (!is.null(size_of_random_subset_for_testing)){
      size_of_random_subset_for_testing <- as.integer(size_of_random_subset_for_testing)
      rsam <- sample(x = seq_len(dim(params_to_test)[1]),
                     size = size_of_random_subset_for_testing,
                     replace = FALSE)
      params_to_test <- params_to_test[rsam, ]
    }
    cat("number of simulations: ",dim(params_to_test)[1],"\n")

    ##############################################################################
    ### Start Simulations
    ##############################################################################
    if (verbose) cat("start simulations ... ")
    if (verbose) {
      # sim_counter <- 1
      cat("\n")
    }

    t0 <- Sys.time()
    outs <- bioc_par_pmap(as.list(params_to_test),
                          BPPARAM=BiocParallel::MulticoreParam(workers = nr_cores,progressbar = TRUE,RNGseed = seed),
                          function(simType, censoring_rate, n_obs,
                                   b0, b1, b2, mi_rep,
                                   cov_dep_cens, error_variance,
                                   variance_raneff, transform_fn,C1,C2,
                                   formula,sim_id){
      tmp_betas <- c(b0=b0,b1=b1,b2=b2)
      tmp_formula <- formula(paste("Y~",formula))
      variables_formula <- extract_variables_from_formula(tmp_formula)
      ############################################################################
      # simulate data
      args <- list(n = n_obs,
                   formula = tmp_formula,
                   n_levels_fixeff = n_levels_fixeff,
                   n_levels_raneff = n_levels_raneff,
                   type = simType,
                   b = tmp_betas,
                   censoring_dependent_on_covariate = cov_dep_cens,
                   error_variance = error_variance,
                   variance_raneff = variance_raneff,
                   number_of_clusters = NULL,
                   number_of_differential_clusters = 1,
                   weibull_params = list(X = list(shape = weibull_params[["shape_x"]], scale = weibull_params[["scale_x"]]),
                                         C = list(shape = weibull_params[["shape_c1"]], scale = C1)),
                   weibull_params_covariate_dependent_censoring = list(shape = weibull_params[["shape_c2"]], scale = C2),
                   transform_fn = transform_fn)

      ## potential errors if multiple differential clusters are wanted, therefore
      ## try multiple times
      max_tries <- 20
      while (max_tries>0) {
        data_sim <- do.call(simulate_singlecluster,args = args)

        cens_rate_eff <- (n_obs-sum(data_sim[[variables_formula$censoring_indicator]]))/n_obs
        # if (cens_rate_eff >= censoring_rate-0.05 & cens_rate_eff <= censoring_rate+0.05){
        if (cens_rate_eff >= censoring_rate-0.2 & cens_rate_eff <= censoring_rate+0.2){
            break
        }
        max_tries <- max_tries-1
      }
      if (max_tries<20){
        cat("max_tries:\t ",max_tries,"\n")
      }
      if (max_tries==0){
        warning("max tries of 20 reached",call. = FALSE)
      }
      if (is.null(data_sim)) stop("no valid dataset generated",call. = FALSE)
      # transform covariates to factor
      for (i in seq_along(variables_formula$covariates)) {
        data_sim[[variables_formula$covariates[i]]] <-
          factor(data_sim[[variables_formula$covariates[i]]],
                 labels = 0:(length(unique(data_sim[[variables_formula$covariates[i]]]))-1))
      }
      # transform random effects covariates to factor
      for (i in seq_along(variables_formula$random_covariates)) {
        data_sim[[variables_formula$random_covariates[i]]] <-
          factor(data_sim[[variables_formula$random_covariates[i]]],
                 labels = 0:(length(unique(data_sim[[variables_formula$random_covariates[i]]]))-1))
      }
      tmp_cov_var <- as.factor(data_sim %>%
                                 dplyr::select(!!variables_formula$covariates) %>%
                                 purrr::as_vector())
      levels(tmp_cov_var) <- seq_along(unique(tmp_cov_var))
      tmp_rancov_var <- as.factor(data_sim %>%
                                    dplyr::select(!!variables_formula$random_covariates) %>%
                                    purrr::as_vector())
      levels(tmp_rancov_var) <- seq_along(unique(tmp_rancov_var))
      data_sim <- data_sim %>%
        dplyr::mutate(!!variables_formula$covariates := tmp_cov_var)
      if (!is.null(variables_formula$random_covariates)){
        data_sim <- data_sim %>%
          dplyr::mutate(!!variables_formula$random_covariates := tmp_rancov_var)
      }
      tmp_formula_2 <- create_glmm_formula(tmp_formula)

      # estimate signal-to-noise ratio
      tmp <- stringr::str_split(string = as.character(tmp_formula_2),pattern = "\\+")
      tmp2 <- stringr::str_replace_all(tmp[[3]]," ","")
      formula_full <- stringr::str_c(tmp[[2]],tmp[[1]],
                                     stringr::str_replace(as.character(tmp_formula_2)[[3]],tmp2[1],"TrVal"))
      formula_reduced <-  stringr::str_c(tmp[[2]],tmp[[1]],stringr::str_c(tmp2[-1],collapse = "+"))
      snr <- tryCatch(general_signal_to_noise(formula_full, formula_reduced, data_sim,
                                              type = simType, weights = data_sim$n_cells),
                      error = function(e) NA)
      effect_of_cov_dep_cens <-
        log_ratio_censoring_per_covariate_level(data_sim,
                                                variables_formula$covariates[1],
                                                variables_formula$censoring_indicator)
      ############################################################################
      ### Imputation
      ############################################################################

      ############################################################################

      data_sim <- data_sim %>% dplyr::mutate_at(dplyr::vars(variables_formula$covariates), ~as.numeric(.)-1)
      ### Run all methods
      out_ls <- map(imputation_types,function(imputation_type){
        out_tib <-
        tryCatch({suppressMessages(suppressWarnings({
            ### Run single method
            tmp_out <- conditional_multiple_imputation(
              data = data_sim,
              id = "ID",
              formula = tmp_formula,
              regression_type = simType,
              mi_reps = mi_rep,
              imputation_method = imputation_type ,
              weights = "n_cells",
              family = "binomial"
            )}))
            censored_variable_est_name <- ifelse(imputation_type %in% c("cc","pmm"),
                                                 variables_formula$censored_variable,
                                                 paste0(variables_formula$censored_variable,"_est"))
            pval_ci <- pval_and_ci_from_fits(tmp_out$fits)[censored_variable_est_name, ]
            # log ratio of number of observed events per level of the covariate
            effect_of_cov_dep_cens <-
              log_ratio_censoring_per_covariate_level(data_sim,
                                                      variables_formula$covariates[1],
                                                      variables_formula$censoring_indicator)
            mean_mse_fits <- mean_mse_of_fits(tmp_out$fits,variables_formula$response)
            out_tib <-
              dplyr::tibble(
                simType = simType,
                formula = formula,
                impType = imputation_type,
                n_dat = n_obs,
                n_cens = length(data_sim[[variables_formula$censoring_indicator]]) -
                  sum(data_sim[[variables_formula$censoring_indicator]]),
                censoring_rate = censoring_rate,
                b0_True = b0,
                b1_True = b1,
                b0 = tmp_out[["betasMean"]][["b0"]],
                b1 = tmp_out[["betasMean"]][["b1"]],
                Varb0 = tmp_out[["betasVar"]][[1]],
                Varb1 = tmp_out[["betasVar"]][[2]],
                pval = pval_ci[["pval"]],
                ci_lower = pval_ci[["ci_lower"]],
                ci_upper = pval_ci[["ci_upper"]],
                total_var = pval_ci[["total_var"]],
                within_var = pval_ci[["within_var"]],
                between_var = pval_ci[["between_var"]],
                df_adj = pval_ci[["df_adj"]],
                riv = pval_ci[["riv"]],
                lambda = pval_ci[["lambda"]],
                fmi = pval_ci[["fmi"]],
                std_err = pval_ci[["std_err"]],
                mean_mse_fits = mean_mse_fits,
                mi_rep = mi_rep,
                cov_dep_cens = cov_dep_cens,
                error_variance = error_variance,
                variance_raneff = variance_raneff,
                transform_fn = transform_fn,
                snr = snr,
                sim_id = sim_id,
                effect_of_cov_dep_cens = effect_of_cov_dep_cens
              )
            for (i in (seq_along(variables_formula$covariates) + 1)) {
              out_tib <- out_tib %>%
                dplyr::mutate(!!paste0("b", i) := tmp_out[["betasMean"]][[paste0("b", i)]],
                              !!paste0("Varb", i) := tmp_out[["betasVar"]][[i +1]],
                              !!paste0("b", i,"_True") := tmp_betas[i+1])
            }
            return(out_tib)
          },
          error = function(e) {})
        return(out_tib)
      })

      out <- tryCatch(dplyr::bind_rows(out_ls), error = function(e){})
      if(verbose) {
        if(sim_id%%(round(dim(params_to_test)[1]/100))==0){
          print(cat("proportion done:",round(sim_id/dim(params_to_test)[1],3),"\t--\tTime:",format(round(Sys.time()-t0,2))))
        }
      }

      return(out)
    })
    if (verbose) cat("\t done\n")
    print(Sys.time()-t0)
    ##############################################################################
    ### End Simulations
    ##############################################################################

    ## combine results
    params <- tryCatch(
      dplyr::bind_rows(outs)
      ,error  = function(e) {
        saveRDS(outs, paste0(dir_save,"/test_results_params.rds"))
        message(paste0("Raw data saved as: ", dir_save,"/test_results_params.rds"))
        outs
      }
    )

    ##############################################################################
    ### data processing
    ##############################################################################
    params_stats <-
      tryCatch(
        simulations_data_processor(params),
        error = function(cnd) {
          NULL
        }
      )
    return(list(params_stats, params,arguments = as.list(match.call())))
  }



simulations_data_processor <- function(params){
  grouping <- c("simType", "formula", "impType", "n_dat", "censoring_rate",
                "b0_True", "b1_True","b2_True","mi_rep", "cov_dep_cens",
                "error_variance", "variance_raneff","transform_fn")
  grouping_quo <- dplyr::quos(simType, formula, impType, n_dat, censoring_rate,
                              b0_True, b1_True,b2_True,mi_rep, cov_dep_cens,
                              error_variance, variance_raneff,transform_fn)
  unselect <- dplyr::quos(pval,snr,ci_lower,ci_upper,total_var,within_var,between_var,
                          df_adj,riv,lambda,fmi,std_err,mean_mse_fits,sim_id,effect_of_cov_dep_cens)
  # params <<-params
  # count occurences
  params_count <- params %>%
    # dplyr::select(-c(!!!unselect)) %>%
    dplyr::group_by(!!!grouping_quo) %>%
    dplyr::tally() %>%
    dplyr::rename(n_sim = n)
  # calculate means for each simulation condition
  params_mean <- params %>%
    dplyr::select(-c(n_cens),-c(!!!unselect)) %>%
    dplyr::group_by(!!!grouping_quo) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::select(-dplyr::starts_with("Varb"))

  params_mean_m <- reshape2::melt(
    params_mean,
    id.vars = grouping,
    variable.name = "parameter",
    value.name = "mean"
  )

  # calculate variance for each simulation condition
  params_var <- params %>%
    dplyr::select(-c(n_cens),-c(!!!unselect)) %>%
    dplyr::group_by(!!!grouping_quo) %>%
    dplyr::summarise_all(var) %>%
    dplyr::select(-dplyr::starts_with("Varb"))
  params_var_m <- reshape2::melt(
    params_var,
    id.vars = grouping,
    variable.name = "parameter",
    value.name = "variance"
  )
  # calculate mean proportion censored for each simulation condition
  params_ncens <- params %>%
    dplyr::mutate(proportion_censored = n_cens/n_dat) %>%
    dplyr::select(-c(n_cens), -dplyr::starts_with("Varb"), -dplyr::starts_with("b"), dplyr::ends_with("True")) %>%
    dplyr::group_by(!!!grouping_quo) %>%
    dplyr::summarise_all(mean)
  # number of distinct simulations
  n_disti_types <- params_count %>% dplyr::n_distinct()
  # combine tables
  params_stats <- dplyr::inner_join(
    params_mean_m,
    params_var_m,
    by = c(grouping,"parameter")
  ) %>%
    dplyr::inner_join(
      params_count,
      by = grouping
    )

  params_stats <- suppressMessages(dplyr::inner_join(params_stats, params_ncens))
  # calculate bias, mse and se
  bias_tib <- params_stats %>%
    dplyr::mutate(parameter = paste0(parameter, "_True")) %>%
    dplyr::select(dplyr::starts_with("b"), parameter, mean)
  bias_vec <- bias_tib %>%
    purrr::pmap_dbl(
      ~ if (..4 == "b0_True") {
        ..5 - ..1
      } else if (..4 == "b1_True") {
        ..5 - ..2
      } else if (..4 == "b2_True") {
        ..5 - ..3
      }
    )
  params_stats <- params_stats %>%
    dplyr::mutate(
      bias = bias_vec,
      mse = variance + bias^2,
      se = sqrt(variance)/sqrt(n_sim),
      ci_low = mean + qt(0.025, n_sim - 1) * se,
      ci_high = mean + qt(0.975, n_sim - 1) * se
    )
  # calculate covarage probability
  tmpmelt <- reshape2::melt(
    params %>%
      dplyr::select(-n_cens, -dplyr::starts_with("Varb"), -c(!!!unselect)),
    id.vars = grouping
  )
  tmpmelt <- dplyr::inner_join( dplyr::select(params_stats, !!!grouping_quo, simType,
                                              parameter, ci_low, ci_high),
                                tmpmelt,
                                by = c(grouping, "parameter" = "variable")
  )
  tmpmelt <- tmpmelt %>%
    dplyr::mutate(
      in_ci = (value >= ci_low & value <= ci_high)
    ) %>%
    dplyr::group_by(!!!grouping_quo,parameter, in_ci) %>%
    dplyr::tally() %>%
    dplyr::ungroup() %>%
    dplyr::distinct(!!!grouping_quo,parameter, .keep_all = TRUE
    ) %>%
    dplyr::rename(n_out_ci = n)
  params_stats <- params_stats %>%
    dplyr::inner_join(
      tmpmelt,
      by = c(grouping, "parameter")
    ) %>%
    dplyr::mutate(
      covProb = ifelse(in_ci, n_out_ci, n_sim - n_out_ci)/n_sim
    ) %>%
    dplyr::select(-c(n_out_ci, in_ci))
  return(params_stats)
}

