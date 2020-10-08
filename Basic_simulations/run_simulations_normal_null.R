# script to run null simulations of single cluster simulations
# save here
dir_save <-"./"

suppressPackageStartupMessages({
  library(diffcyt)
  library(magrittr)
  library(purrr)
  library(SummarizedExperiment)
})
devtools::load_all("../simulationStudyCensoredDiffcyt")
RhpcBLASctl::blas_set_num_threads(1)
covariates_type <- "params_FlowCAP"

if (covariates_type == "params_FlowCAP"){
  # parameters of weibull fit
  weifit_ls <- readRDS(paste0(dir_save,"../FlowCAP_IV/results/weibull_fits_FlowCAP.rds"))
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

for (i in 1:3){
  out <- run_simulations_wrapper_singlecluster(reps = 1000,
                                               nr_cores = 10,
                                               n_obs_c = c(100),
                                               simDataTypes = "glmer",
                                               formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
                                               n_levels_fixeff = c(2),
                                               n_levels_raneff = NULL,
                                               betas = list(b0=c(-2),b1=c(0),b2=c(0)),
                                               censoring_rates = c(0.3,0.5,0.7),
                                               imputation_types = c("mrl","km","rs","cc","km_exp","pmm"),
                                               run_with_covariate_dependent_censoring = FALSE,
                                               mi_reps = c(50),
                                               weibull_params = weibull_params,
                                               error_variance = 0,
                                               variance_raneff = c(1),
                                               transform_fn = c("identity"),
                                               verbose = FALSE,
                                               seed = 123+(i-1),
                                               size_of_random_subset_for_testing = NULL)
  saveRDS(out, paste0(dir_save,"/res_sclu_output_",covariates_type,"_complete_nullsim_",i,".rds"))
}

