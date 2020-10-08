################################################################################
## load packages and data
library(tidyverse)
library(yardstick)
library(iCOBRA)
library(ggpubr)
library(scales)

alpha_factor <- 5
clustering <- ""
# clustering <- "som100"
# covariates_type <- ""
covariates_type <- "params_FlowCAP"
# sim_type <- "complete"
# sim_type <- "complete_group"
sim_type <- "complete_group_test_group"
# transform_fn <- "log_positive"
# transform_fn <- "div_100"
transform_fn <- ""
is_group_bool <- ifelse(sim_type %in% c("complete_group","complete_group_test_group"),TRUE,FALSE)

clustering <- ifelse(covariates_type=="",clustering,paste0(clustering,"_",covariates_type))
clu_name <- ifelse(clustering %in% c("","_"),"meta20",stringr::str_replace(clustering,"^_+",""))
clu_name <- ifelse(clustering == "_params_FlowCAP","meta20_params_FlowCAP",clu_name)
plot_dir_clustering <- paste0("~/Documents/sherborne/Shared_portmacquarie/retger/simulation_study/results/",clu_name,"/")
if (!dir.exists(plot_dir_clustering)){
  dir.create(plot_dir_clustering)
}
plot_dir <- paste0("~/Documents/sherborne/Shared_portmacquarie/retger/simulation_study/results/",clu_name,"/af_",alpha_factor,"/")
if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
}
# dir_save <- "~/Documents/sherborne/Shared_portmacquarie/retger/simulation_study/results/"
dir_save <- plot_dir

# test_results <- readRDS(paste0(dir_save,"res_mclu_output_af_",alpha_factor,ifelse(clustering=="","","_"),clu_name,"_",sim_type,".rds"))
test_results <- readRDS( paste0(dir_save,"res_mclu_output_af_",alpha_factor,"_",clu_name,"_",sim_type,ifelse(transform_fn=="","",paste0("_",transform_fn)),".rds"))
################################################################################
## processing and plotting
all_impTypes <-names(test_results[[2]][[1]])

test_results[[2]] <- lapply(test_results[[2]],function(e) {
  names(e)[names(e)=="km_exp"] <- "kme"
  names(e)[names(e)=="GLMM_og"] <- "ncGLMM"
  return(e)
})

# select methods
if (is_group_bool){
  if (sim_type=="complete_group"){
    all_impTypes <- c("GLMM","cc","km","mrl","rs","kme","pmm")
    # all_impTypes <- c("GLMM","cc","km","mrl","rs","km_exp","km_os","km_wei")
  } else {
    all_impTypes <- c("GLMM","cc","km","mrl","rs","kme","pmm","ncGLMM")
    # all_impTypes <- c("GLMM","cc","km","mrl","rs","km_exp","km_os","km_wei","GLMM_og")
  }
  tib_selection <- quos(n_dat,censoring_rate,cov_dep_cens,transform_fn,nr_diff,slope,group_slope)
} else {
  all_impTypes <- c("GLMM","cc","km","mrl","rs","kme","pmm","coxph")
  # all_impTypes <- c("GLMM","cc","km","mrl","rs","km_exp","km_os","km_wei","coxph")
  tib_selection <- quos(n_dat,censoring_rate,cov_dep_cens,transform_fn,nr_diff,slope)
}
cat("Methods used for plotting: \n  ",paste0(all_impTypes,collapse = " - "),"\n")

# get all simulation conditions
params <- purrr::map(test_results[[2]], ~ .x[c(all_impTypes)])
conditions <- purrr::map(params, function(dfs){
  ind <- which(unlist(purrr::map(dfs,~is.data.frame(.x))))[1]
  if (is.na(ind) & is_group_bool) {
    return(tibble::tibble(n_dat=NA,censoring_rate=NA,cov_dep_cens=NA,transform_fn=NA,nr_diff=NA,slope=NA,group_slope=NA))
  } else if(is.na(ind) & !is_group_bool){
    return(tibble::tibble(n_dat=NA,censoring_rate=NA,cov_dep_cens=NA,transform_fn=NA,nr_diff=NA,slope=NA))
  }
  dfs[[ind]] %>% dplyr::select(!!!tib_selection) %>% unique()
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(ID=seq_along(n_dat)) %>%
  dplyr::arrange(n_dat,censoring_rate,transform_fn,cov_dep_cens)

# unique simulation conditions
all_conds <- conditions %>%
  select(-ID) %>%
  group_by(!!!tib_selection) %>%
  tally()
n_dat <- unique(all_conds$n_dat)
nr_diff <- unique(all_conds$nr_diff)
censrate <- unique(all_conds$censoring_rate)
transform_fn <- unique(all_conds$transform_fn)
cov_dep_cens <- unique(all_conds$cov_dep_cens)
slope <- unique(all_conds$slope)
if (is_group_bool){
  group_slope <- unique(all_conds$group_slope)
}

# unique conditions as a list
arg_list <- as.list(expand.grid(n_sam=n_dat,
                                censrate=censrate,
                                cens_mech=cov_dep_cens,
                                trans_fn=transform_fn,
                                nr_diff=nr_diff,
                                slope=slope))
if (is_group_bool){
  arg_list <- as.list(expand.grid(n_sam=n_dat,
                                  censrate=censrate,
                                  cens_mech=cov_dep_cens,
                                  trans_fn=transform_fn,
                                  nr_diff=nr_diff,
                                  slope=slope,
                                  group_slope=unique(all_conds$group_slope)))
}
attributes(arg_list) <- NULL

# unique conditions as a data frame, used for plotting
arg_df <- data.frame(id = seq_along(arg_list[[1]]),
                     n_sam=arg_list[[1]],
                     censrate=arg_list[[2]],
                     cens_mech=arg_list[[3]],
                     trans_fn=arg_list[[4]],
                     nr_diff_1=arg_list[[5]],
                     slope_1=arg_list[[6]])
if (is_group_bool){
  arg_df <- cbind(arg_df, group_slope_1=arg_list[[7]])
}

# id's of simulations per simulation condition, for aggregating the plots
if (is_group_bool){
  cond_filt <- purrr::pmap(arg_list,function(n_sam,censrate,cens_mech,trans_fn,nr_diff_1,slope_1,group_slope_1){
    conditions %>%
      dplyr::filter(
        censoring_rate==censrate,
        cov_dep_cens==cens_mech,
        transform_fn==trans_fn,
        n_dat==n_sam,
        nr_diff==nr_diff_1,
        slope==slope_1,
        group_slope==group_slope_1
      ) %>%
      dplyr::select(ID) %>%
      purrr::as_vector()
  })
} else {
  cond_filt <- purrr::pmap(arg_list,function(n_sam,censrate,cens_mech,trans_fn,nr_diff_1,slope_1){
    conditions %>%
      dplyr::filter(
        censoring_rate==censrate,
        cov_dep_cens==cens_mech,
        transform_fn==trans_fn,
        n_dat==n_sam,
        nr_diff==nr_diff_1,
        slope==slope_1
      ) %>%
      dplyr::select(ID) %>%
      purrr::as_vector()
  })
}

# number of clusters per simulation condition, to get same lenghts if NA's present
length_theo <- length(cond_filt[[1]]) * dim(test_results[[2]][[1]][[1]])[1]

# generate all plots
list_data <-   purrr::map(seq_along(cond_filt),function(i){
  # combined datasets from multiple repetitions
  params_one_comb <- suppressMessages(purrr::map(cond_filt[[i]],function(x){
    params_one <- params[[x]]
    params_one <- purrr::map(seq_along(params_one), function(x){
      if(is(params_one[[x]],"data.frame")){
        return(params_one[[x]])
      }else{
        NA_tib <- params_one$cc %>%
          mutate(p_val=NA,p_adj=NA,impType=names(params_one)[x])
        return(NA_tib)
      }
    })
    purrr::map(params_one, ~bind_rows(.x)) %>% bind_rows()
  }) %>%
    dplyr::bind_rows())

  # as tibble
  params_one_comb$impType[params_one_comb$impType=="km_exp"] <- "kme"
  params_one_comb$impType[params_one_comb$impType=="GLMM_og"] <- "ncGLMM"
  tib2plt <- params_one_comb %>%
    dplyr::mutate(proportion_censored = n_cens/n_dat,
                  truestru = !is.na(paired),
                  class = as.factor(truestru),
                  impType = factor(impType, levels = all_impTypes),
                  simType = factor(simType, levels = c("glmer")),
                  n_dat= as.factor(n_dat),
                  cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR")))

  # to cobradata for use in cobra functions
  cobratib <- tib2plt %>%
    mutate(id = seq_len(dim(tib2plt)[1]))
  cobratib$impType <- droplevels(cobratib$impType)
  tmp_impTypes <- c()
  pval <- data.frame(km = cobratib$p_val[cobratib$impType=="km"])
  padj <- data.frame(km = cobratib$p_adj[cobratib$impType=="km"])
  for (impType in unique(cobratib$impType)) {
    na_ind <- is.na(cobratib$p_adj[cobratib$impType==impType])
    print(paste(impType,":",sum(na_ind),"NA's"))
    if (sum(na_ind)<length_theo){
      pval[[impType]] <- cobratib$p_val[cobratib$impType==impType]
      padj[[impType]] <- cobratib$p_adj[cobratib$impType==impType]
      tmp_impTypes <- c(tmp_impTypes,impType)
    }
  }
  tmp_impTypes <- sort(tmp_impTypes)
  rownames(pval) <- paste0("I_",cobratib$id[cobratib$impType=="km"])
  rownames(padj) <- paste0("I_",cobratib$id[cobratib$impType=="km"])

  truth <- data.frame(status = as.integer(cobratib$class[cobratib$impType=="km"])-1,
                      cov_dep_cens = cobratib$cov_dep_cens[cobratib$impType=="km"],
                      sim_id = cobratib$sim_id[cobratib$impType=="km"],
                      row.names = paste0("I_",cobratib$id[cobratib$impType=="km"]),
                      n_sam=arg_df$n_sam[i],
                      censrate=arg_df$censrate[i]
  )
  return(list(pval=pval,padj=padj,truth=truth))
})
#

plts <-   purrr::map(seq_along(list_data),function(i){
  cobradata <- COBRAData(pval = list_data[[i]][["pval"]], truth = list_data[[i]][["truth"]],padj = list_data[[i]][["padj"]])
  # calculate performance
  cobraperf <- calculate_performance(cobradata,
                                     binary_truth = "status",
                                     aspects = c("fdrtpr", "fdrtprcurve",
                                                 "tpr", "roc", "fpr","overlap"),
                                     thrs = c(0.01, 0.05, 0.1))
  print(names(list_data[[i]]$pval))
  color_ind <- unlist(purrr::map(sort(names(list_data[[i]]$pval)), ~which(.x==all_impTypes)))
  colorscheme <- c(RColorBrewer::brewer.pal(12,"Paired")[-c(7,8,9)][color_ind])
  if ("coxph" %in% names(list_data[[i]]$pval)){
    colorscheme[which(color_ind == 8)] <- "black"
  }
  cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme,
                                     facetted = TRUE)
})
saveRDS(list(cobradata = plts, arg_df = arg_df),paste0(plot_dir,"multicluster_cobra_data_af_",alpha_factor,"_",clustering,"_",sim_type,ifelse(transform_fn=="","",paste0("_",transform_fn)),".rds"))


# pval <- purrr::map(list_data, ~.x[["pval"]]) %>% reduce(rbind)
# padj <- purrr::map(list_data, ~.x[["padj"]]) %>% reduce(rbind)
# truth <- purrr::map(list_data, ~.x[["truth"]]) %>% reduce(rbind)
# cobradata <- COBRAData(pval = pval, truth = truth,padj = padj)
# calculate performance
# cobraperf <- calculate_performance(cobradata,
#                                    binary_truth = "status",
#                                    aspects = c("fdrtpr", "fdrtprcurve",
#                                                "tpr", "roc", "fpr","overlap"),
#                                    thrs = c(0.01, 0.05, 0.1),
#                                    splv = c("n_sam"))
# print(names(pval))
# color_ind <- unlist(purrr::map(names(pval), ~which(.x==all_impTypes)))
# colorscheme <- c(RColorBrewer::brewer.pal(12,"Paired")[-c(7,8,9)][color_ind])
# cobraplot <- prepare_data_for_plot(cobraperf, colorscheme = colorscheme,
#                                    facetted = TRUE)
# plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)



