#' Wrapper for running simulations
#'
#' @param reps number of repetitions for each condition
#' @param nr_cores number of cores to use for parallelization
#' @param n_obs_c vector of data sizes
#' @param simDataTypes vector of simulation types, options are "lm", "glm","glmer"
#' @param reference_data reference data set to use, from which alpha and size arguments are inferred
#'  see \code{\link[diffcyt]{simulate_multicluster}}.
#' @param nr_diff Positive Integer. Total number of clusters
#'  with a true signal. see \code{\link[diffcyt]{simulate_multicluster}}.
#' @param alpha parameter of dirichlet multinomial distribution see \code{\link[diffcyt]{simulate_multicluster}}.
#' @param sizes parameter of dirichlet multinomial distribution see \code{\link[diffcyt]{simulate_multicluster}}.
#' @param slope list of slope arguments to test, see \code{\link[diffcyt]{simulate_multicluster}}.
#' @param group list of group arguments to test, see \code{\link[diffcyt]{simulate_multicluster}}.
#' @param group_slope list of group_slope arguments to test, see \code{\link[diffcyt]{simulate_multicluster}}.
#' @param diff_cluster see \code{\link[diffcyt]{simulate_multicluster}}
#' @param enforce_sum_alpha see \code{\link[diffcyt]{simulate_multicluster}}
#' @param simDataTypes default = "glmer". The linear model to test.
#' @param formulas formulas specifing data, details see  \code{\link[diffcyt]{simulate_singlecluster}}
#' @param censoring_rates vector of approximate censoring rates used. effective
#'  censoring rate is enforced to be plus minus 5 percent of this value.
#' @param weibull_params list with elements shape_x,shape_c1,shape_2,scale_x for the weibull parameters
#'  to use for sampling the censored covariate. scale_c1 is inferred based on the given censoring rate.
#' @param imputation_types vector of \code{\link[diffcyt]{conditional_multiple_imputation}}
#'  types used for inference, possible options are "cc", "mrl", "km", "rs", "pmm","km_exp","km_os","km_wei"
#' @param mi_reps number of repetitions in \code{\link[diffcyt]{conditional_multiple_imputation}}
#' @param transform_fn vector of functions to transform censored covariate or
#'  'identity' (no transformation) or 'boxcox' (box-cox transformation).
#'   default = 'identity'
#' @param contrast default = c(0,1,0). What to test. first element is the intercept, second the censored covariate,
#'  third the second covariate (if specified).
#' @param verbose verbose
#' @param seed positive integer, random seed
#' @param size_of_random_subset_for_testing number of simulations for test run
#' @importFrom magrittr %>%
#' @export
run_simulations_wrapper_multicluster <-
  function(reps = 2,
           nr_cores = 1,
           n_obs_c = c(50),
           reference_data = NULL,
           nr_diff = 6,
           alpha = NULL,
           sizes = NULL,
           slope = FALSE,
           group = FALSE,
           group_slope = FALSE,
           diff_cluster = list(c(7,5),c(20,4),c(11,13)),
           enforce_sum_alpha = FALSE,
           simDataTypes = c("glmer"),
           formulas = list(glmer = list(formula(Y~Surv(X,I) +(1|R)))),
           censoring_rates = c(0.3),
           weibull_params = list(shape_x = 0.5,
                                 shape_c1 = 1,
                                 shape_c2 = 0.5,
                                 scale_x = 0.25),
           imputation_types = c("mrl","km","rs"),
           mi_reps = 10,
           transform_fn = "identity",
           contrast = NULL,
           verbose = FALSE,
           seed = 123,
           size_of_random_subset_for_testing = NULL) {
    # directory to save tempfiles
    dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
    dir.create(dir_save, showWarnings = FALSE)
    filenames <- c()

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)

    # infer scale parameter of censoring distribution dependent on censoring rate
    if (verbose) cat("start censoring parameter estimation ... ")
    tmp_df_new <- scales_for_censoring(
      censoring_val = censoring_rates,
      log_ratio_val = c(0),
      seq_to_test = c(seq(0.001,1,length.out = 100),seq(1,20000,by=10)),
      n = 10000,
      shapes = list(shape_x=weibull_params[["shape_x"]],
                    shape_c1=weibull_params[["shape_c1"]],
                    shape_c2=weibull_params[["shape_c2"]]),
      scale_x = weibull_params[["scale_x"]]
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
      nr_diff_1 = nr_diff,
      slope_1=slope,
      group_1=group,
      group_slope_1=group_slope,
      mi_rep = mi_reps,
      transform_fn = transform_fn,
      stringsAsFactors = FALSE
    )
    params_to_test <- suppressWarnings(
      dplyr::inner_join(params_to_test,
                        tmp_df_new,
                        by = c("censoring_rate" = "censoring_rate")
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

    # random subset for testing
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
                          function(simType, censoring_rate, n_obs, nr_diff_1,
                                   slope_1, group_1,group_slope_1,mi_rep,
                                   transform_fn, cov_dep_cens,C1,C2,formula,sim_id){
      # some changes for sampling
      if (is.na(slope_1)) slope_1 <- NULL
      if (is.na(group_1)) group_1 <- NULL
      if (is.na(group_slope_1)) group_slope_1 <- NULL
      if (length(slope_1)>1) slope_1 <- as.list(slope_1)
      if (length(group_slope_1)>1) group_slope_1 <- as.list(group_slope_1)

      # uncensored formula
      tmp_formula <- formula(paste("Y~",formula))
      # variable names used
      variables_formula <- extract_variables_from_formula(tmp_formula)
      ############################################################################
      # simulate data
      # simulate covariate
      ## potential errors if multiple differential clusters are wanted, therefore
      ## try multiple times
      if (transform_fn=="div_100"){
        transform_fn_insim <- function(x){x/100}
      } else {
        transform_fn_insim <- transform_fn
      }
      max_tries <- 100
      while (max_tries>0) {
        one_cluster_data <- simulate_singlecluster(n = n_obs,
                                          formula = tmp_formula,
                                          type = simType,
                                          weibull_params = list(X = list(shape = weibull_params[["shape_x"]], scale = weibull_params[["scale_x"]]),
                                                                C = list(shape = weibull_params[["shape_c1"]], scale = C1)),
                                          censoring_dependent_on_covariate = cov_dep_cens,
                                          weibull_params_covariate_dependent_censoring = list(shape = weibull_params[["shape_c2"]], scale = C2),
                                          transform_fn = transform_fn_insim)

        cens_rate_eff <- (n_obs-sum(one_cluster_data[[variables_formula$censoring_indicator]]))/n_obs
        if (cens_rate_eff >= censoring_rate-0.05 & cens_rate_eff <= censoring_rate+0.05){
          break
        }
        max_tries <- max_tries-1
      }
      if (max_tries<20){
        cat("max_tries:\t ",max_tries,"\n")
      }
      if (max_tries==0){
        warning("max tries of 100 reached",call. = FALSE)
      }

      # resample sizes if not correct size
      if (!is.null(sizes) & length(sizes) != n_obs){
        sizes <- sample(x = sizes,size = n_obs,replace = TRUE)
      }
      # arguments list for simulating data
      args_ls <- list(counts = reference_data,
                   nr_diff = nr_diff_1,
                   nr_samples = n_obs,
                   alphas = alpha,
                   sizes = sizes,
                   covariate = one_cluster_data[["TrVal"]],
                   slope = slope_1,
                   group = group_1,
                   group_slope = group_slope_1,
                   diff_cluster = diff_cluster,
                   enforce_sum_alpha = enforce_sum_alpha)

      # simulate counts
      data_sim_comp <- do.call(simulate_multicluster,args=args_ls)

      # some data manipulations for input into testing
      colnames(data_sim_comp$row_data)[1] <- "cluster_id"
      d_counts <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts=data_sim_comp$counts),
        rowData = data_sim_comp$row_data,
        colData = data_sim_comp$col_data)
      data_sim <- as.data.frame(SummarizedExperiment::colData(d_counts))
      colnames(data_sim)[1] <- variables_formula$random_covariates
      colnames(data_sim)[2] <- paste0(variables_formula$censored_variable,"_True")
      if (dim(data_sim)[2] > 2){
        colnames(data_sim)[3:dim(data_sim)[2]] <- variables_formula$covariates
      }
      data_sim[[variables_formula$censoring_indicator]] <- one_cluster_data[[variables_formula$censoring_indicator]]
      data_sim[[variables_formula$censored_variable]] <- one_cluster_data[[variables_formula$censored_variable]]

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

      # log ratio of censoring
      if (!is.null(variables_formula$covariates)){
        effect_of_cov_dep_cens <-
          log_ratio_censoring_per_covariate_level(data_sim,
                                                  variables_formula$covariates[1],
                                                  variables_formula$censoring_indicator)
        }else {
        effect_of_cov_dep_cens <- NA
        }

      da_formula <- diffcyt::createFormula(
        data_sim,
        cols_fixed = c(variables_formula$censored_variable, variables_formula$covariates),
        cols_random = c(variables_formula$random_covariates),event_indicator = variables_formula$censoring_indicator
      )
      if (is.null(contrast)){
        contrast <- diffcyt::createContrast(c(0, 1,rep(0,length(variables_formula$covariates))))
      }

      ############################################################################
      ### Imputation
      ############################################################################

      ############################################################################
      # conditions of simulation
      conditions_tib <- dplyr::tibble(
        simType = simType,
        formula = formula,
        n_dat = n_obs,
        nr_diff = nr_diff_1,
        n_cens = length(data_sim[[variables_formula$censoring_indicator]]) -
          sum(data_sim[[variables_formula$censoring_indicator]]),
        censoring_rate = censoring_rate,
        mi_rep = mi_rep,
        cov_dep_cens = cov_dep_cens,
        transform_fn = transform_fn,
        effect_of_cov_dep_cens = effect_of_cov_dep_cens,
        sim_id = sim_id)
      if (!is.null(slope_1)){
        conditions_tib <- dplyr::mutate(conditions_tib,slope=paste(unlist(slope_1),collapse = ","))
      }
      if (!is.null(group_1)){
        conditions_tib <- dplyr::mutate(conditions_tib,group=paste(unlist(group_1),collapse = ","))
      }
      if (!is.null(group_slope_1)){
        conditions_tib <- dplyr::mutate(conditions_tib,group_slope=paste(unlist(group_slope_1),collapse = ","))
      }


      ### Run all methods
      out_ls <- purrr::map(imputation_types,function(imputation_type){
        out <- tryCatch({
          ## Run single method
          out <-
            testDA_censoredGLMM(
              d_counts = d_counts,
              formula = da_formula,
              contrast = contrast,
              imputation_method = imputation_type,
              verbose = FALSE,
              mi_reps = mi_rep
            )
          # add all metadata
          out_rowdata <- dplyr::as_tibble(SummarizedExperiment::rowData(out))
          out_conditions <- conditions_tib %>%
            dplyr::mutate(impType = imputation_type)
          cbind(out_rowdata, out_conditions)
        },
        error = function(e){e})
        return(out)
      })
      names(out_ls) <- imputation_types
      ## test coxph if no second variable
      if (is.null(group_1)){
        out_ls[["coxph"]] <- tryCatch({
          counts <- assay(d_counts)
          proportions <- apply(counts,2,function(x) x/sum(x))
          covariate <- data_sim[[variables_formula$censored_variable]]
          obs_ind <- data_sim[[variables_formula$censoring_indicator]]
          p_val <- sapply(seq_len(dim(proportions)[1]), function(i){
            coef(summary(survival::coxph(formula = formula(survival::Surv(covariate,obs_ind)~proportions[i,]))))[5]
          })
          p_adj <- p.adjust(p_val)
          out_conditions <- cbind(conditions_tib,dplyr::tibble(p_val = p_val,p_adj = p_adj,impType = "coxph"))
          cbind(out_conditions,as.data.frame(SummarizedExperiment::rowData(d_counts)))
        },error = function(e) {e})
      }
      # run uncensored method
      da_formula_uncens <- diffcyt::createFormula(
        data_sim,
        cols_fixed = c(paste0(variables_formula$censored_variable,"_True"), variables_formula$covariates),
        cols_random = c(variables_formula$random_covariates))

      out_ls[["GLMM"]] <- tryCatch({
        out_GLMM <-testDA_GLMM(
          d_counts = d_counts,
          formula = da_formula_uncens,
          contrast = contrast,
          min_samples = 5)
        cbind(cbind(dplyr::as_tibble(SummarizedExperiment::rowData(out_GLMM)),data_sim_comp$row_data),
                                  dplyr::mutate(conditions_tib,impType = "GLMM"))
      }, error = function(e) {e})

      da_design <- createDesignMatrix(
        data_sim,
        c(paste0(variables_formula$censored_variable,"_True"), variables_formula$covariates)
      )
      out_ls[["edgeR"]] <- tryCatch({
        out_edgeR <- testDA_edgeR(
          d_counts = d_counts,
          design = da_design,
          contrast = contrast,
          min_samples = 5)
        cbind(cbind(dplyr::as_tibble(SummarizedExperiment::rowData(out_edgeR)),data_sim_comp$row_data),
                                  dplyr::mutate(conditions_tib,impType = "edgeR"))
      }, error = function(e) {e})

      out_ls[["voom"]] <-  tryCatch({
        out_voom <- testDA_voom(
          d_counts = d_counts,
          design = da_design,
          contrast = contrast,
          min_samples = 5)
        cbind(cbind(dplyr::as_tibble(SummarizedExperiment::rowData(out_voom)),data_sim_comp$row_data),
                                   dplyr::mutate(conditions_tib,impType = "voom"))
      }, error = function(e) {e})

      # if testing for group assoc
      if (c(contrast)[2] == 0){
        # run uncensored method
        da_formula_uncens_og <- diffcyt::createFormula(
          data_sim,
          cols_fixed = c(variables_formula$covariates),
          cols_random = c(variables_formula$random_covariates))

        contrast_og  <- diffcyt::createContrast(c(0, 1))

        out_ls[["GLMM_og"]] <- tryCatch({
          out_GLMM <-testDA_GLMM(
            d_counts = d_counts,
            formula = da_formula_uncens_og,
            contrast = contrast_og,
            min_samples = 5)
          cbind(cbind(dplyr::as_tibble(SummarizedExperiment::rowData(out_GLMM)),data_sim_comp$row_data),
                dplyr::mutate(conditions_tib,impType = "GLMM_og"))
        }, error = function(e) {e})
      }


      return(out_ls)
    })
    if (verbose) cat("\t done\n")
    print(Sys.time()-t0)
    ##############################################################################
    ### End Simulations
    ##############################################################################

    ## combine results
    params <- tryCatch(
      if (is.null(number_of_clusters)){
        dplyr::bind_rows(outs)
      } else {
        tryCatch(dplyr::bind_rows(outs),error = function(e) outs)
      },
      error  = function(e) {
        saveRDS(outs, paste0(dir_save,"/test_results_params.rds"))
        message(paste0("Raw data saved as: ", dir_save,"/test_results_params.rds"))
        outs
      }
    )

    return(list(NULL, params,arguments = as.list(match.call())))
  }


# pmap parallel with BiocParallel
bioc_par_pmap <- function(.l, .f, ..., BPPARAM = BiocParallel::MulticoreParam(workers = 1)) {
  .f <- purrr::as_mapper(.f, ...)
  do.call(BiocParallel::bpmapply, c(.l, list(FUN = .f, MoreArgs = list(...), SIMPLIFY = FALSE, BPPARAM = BPPARAM)))
}
