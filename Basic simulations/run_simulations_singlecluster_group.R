# script to run single cluster simulations
dir_save <-"./"

suppressPackageStartupMessages({
  library(diffcyt)
  library(magrittr)
  library(purrr)
  library(SummarizedExperiment)
})
devtools::load_all("simulationStudyCensoredDiffcyt")
RhpcBLASctl::blas_set_num_threads(1)
covariates_type <- "params_FlowCAP"

if (covariates_type == "params_FlowCAP"){
  # parameters of weibull fit
  weifit_ls <- readRDS(paste0(dir_save,"weibull_fits_FlowCAP.rds"))
  weibull_params <- list(shape_x = weifit_ls$weifit$estimate[1],
                         shape_c1 = weifit_ls$weifit_cens$estimate[1],
                         shape_c2 = weifit_ls$weifit_cens$estimate[1],
                         scale_x = weifit_ls$weifit$estimate[2])
} else {
  weibull_params <-  list(shape_x = 0.5,
                        shape_c1 = 1,
                        shape_c2 = 0.5,
                        scale_x = 0.25)
}
cat("start simulations ")

standard_params = list(sample_size = 100,
                       censoring_rate = 0.5,
                       beta1 = -1e-4,
                       mi_repetitions = 50,
                       random_effect_variance = 1)

expanded_params = list(sample_size = c(50, 100, 200),
                       censoring_rate = c(0.3, 0.5, 0.7),
                       beta1 = c(-5e-5,-1e-4,-5e-4),
                       mi_repetitions = c(10,50,100),
                       random_effect_variance = c(0.5, 1, 2))

random_seeds <- list(sample_size = 123,
                     censoring_rate = 124,
                     beta1 = 125,
                     mi_repetitions = 126,
                     random_effect_variance = 127)

for(current_sim_param in names(standard_params)){
  standard_params_copy <- standard_params
  standard_params_copy[[current_sim_param]] <- expanded_params[[current_sim_param]]
  out <- run_simulations_wrapper_singlecluster(reps = 100,
                                               nr_cores = 20,
                                               n_obs_c = standard_params_copy$sample_size,
                                               simDataTypes = "glmer",
                                               formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
                                               n_levels_fixeff = c(2),
                                               n_levels_raneff = NULL,
                                               betas = list(b0=c(-2),b1=standard_params_copy$beta1,b2=c(1)),
                                               censoring_rates = standard_params_copy$censoring_rate,
                                               imputation_types = c("mrl","km","rs","cc","km_exp","km_wei","pmm"),
                                               run_with_covariate_dependent_censoring = FALSE,
                                               mi_reps = standard_params_copy$mi_repetitions,
                                               weibull_params = weibull_params,
                                               error_variance = 0,
                                               variance_raneff = standard_params_copy$random_effect_variance,
                                               transform_fn = c("identity"),
                                               verbose = FALSE,
                                               seed = random_seeds[[current_sim_param]],
                                               size_of_random_subset_for_testing = NULL)
  saveRDS(out, paste0(dir_save,"/res_sclu_output_",covariates_type,"_complete_group_",current_sim_param,"_100reps.rds"))
}

