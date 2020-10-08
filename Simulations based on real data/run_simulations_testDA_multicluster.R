# multiple cluster simulations, one (censored) covariate

suppressPackageStartupMessages({
  library(diffcyt)
  library(magrittr)
  library(purrr)
  library(SummarizedExperiment)
})

devtools::load_all("simulationStudyCensoredDiffcyt")
RhpcBLASctl::blas_set_num_threads(1)
dir_save <-"./"
transform_fn <- "identity"
clustering <- ""
covariates_type <- "params_FlowCAP"

cat("clustering to use: \t",clustering,"\n")
if (clustering == "som100"){
  da_res <- readRDS("../FlowCAP_4/results/da_res1_som100_cc_complete.rds")
} else {
  da_res <- readRDS("../FlowCAP_4/results/da_res1_cc_complete.rds")
}
DM_fit_file_name <- ifelse(clustering=="som100","DM_fit_FlowCAP_som100.rds","DM_fit_FlowCAP.rds")
if (DM_fit_file_name %in% list.files(dir_save)){
  dirfit <- readRDS(paste0(dir_save,DM_fit_file_name))
} else {
  cat("fit dirichlet multinomial on FlowCAP_4 data\n")
  dirfit <- dirmult::dirmult(t(SummarizedExperiment::assay(da_res$d_counts)))
  saveRDS(dirfit, paste0(dir_save,DM_fit_file_name))
}

if (covariates_type == "params_FlowCAP"){
  # parameters of weibull fit
  weifit_ls <- readRDS(paste0(dir_save,"weibull_fits_FlowCAP.rds"))
  weibull_params <- list(shape_x = weifit_ls$weifit$estimate[1],
                         shape_c1 = weifit_ls$weifit_cens$estimate[1],
                         shape_c2 = weifit_ls$weifit_cens$estimate[1],
                         scale_x = weifit_ls$weifit$estimate[2])
} else {
  weibull_params <- NULL
}

clustering <- ifelse(covariates_type=="",clustering,paste0(clustering,"_",covariates_type))
clu_name <- ifelse(clustering %in% c("","_"),"meta20",stringr::str_replace(clustering,"^_+",""))
clu_name <- ifelse(clustering == "_params_FlowCAP","meta20_params_FlowCAP",clu_name)

for(alpha_factor in c(1,2,5,10)){
  plot_dir_clustering <- paste0(dir_save,"/",clu_name,"/")
  if (!dir.exists(plot_dir_clustering)){
    dir.create(plot_dir_clustering)
  }
  plot_dir <- paste0(dir_save,"/",clu_name,"/af_",alpha_factor,"/")
  if (!dir.exists(plot_dir)){
    dir.create(plot_dir)
  }
  alphas <- dirfit$gamma * alpha_factor
  sizes <- apply(t(SummarizedExperiment::assay(da_res$d_counts)), 1, sum)
  cov_ls <- lapply(seq_along(sizes), function(k) cov_matrix_dirichlet_multinomial(alphas,sizes[k]))
  names(cov_ls) <- paste0("s_",round(sizes))
  saveRDS(list(alphas=alphas,sizes=sizes,cov=cov_ls), paste0(dir_save,"/",clu_name,"/af_",alpha_factor,"/res_mclu_params_af_",alpha_factor,"_",clustering,"_params_FlowCAP.rds"))

  cat("start simulations for alpha_factor:\t",alpha_factor,"\n")
  out <- run_simulations_wrapper_multicluster(reps=10,
                                              nr_cores = 10,
                                              n_obs_c = c(50, 100, 200, 400),
                                              nr_diff = 6,
                                              mi_reps = 50,
                                              alpha = alphas,
                                              sizes = sizes,
                                              slope=c(list(rep(0.9,3))),
                                              group=NA,
                                              group_slope=NA,
                                              diff_cluster = list(c(7,5),c(20,4),c(11,13)),
                                              enforce_sum_alpha = FALSE,
                                              censoring_rates = c(0.3,0.5,0.7),
                                              weibull_params = weibull_params,
                                              imputation_types = c("km","rs","mrl","cc","pmm","km_exp"),
                                              formulas = list(glmer = list(formula(Y~Surv(X,I) + (1|R)))),
                                              verbose = TRUE, seed = 123,
                                              transform_fn=c(transform_fn))
  saveRDS(out, paste0(dir_save,"/",clu_name,"/af_",alpha_factor,"/res_mclu_output_af_",alpha_factor,"_",clu_name,"_complete_",transform_fn,".rds"))
}




