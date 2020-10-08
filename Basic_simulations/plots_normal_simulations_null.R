# plots for null simulations for single cluster


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)
covariates_type <- "params_FlowCAP"

dir_save <- "./"
plot_dir <- dir_save
plot_width <- 20
plot_height <- 15

test_results_proc <- lapply(1:3,function(i){
  tmpres <- readRDS(paste0(dir_save,"/res_sclu_output_",covariates_type,"_complete_nullsim_",i,".rds"))
  tmpres[[2]] %>%
    filter(!impType %in% c("km_os","km_wei")) %>%
    mutate(simulation_run=i)
}) %>% reduce(rbind)

test_results_proc$impType[test_results_proc$impType=="km_exp"] <- "kme"
tmp_impTypes <- test_results_proc$impType %>% unique() %>% sort()
all_impTypes <- c("GLMM","cc","km","mrl","rs","kme","pmm")
color_ind <- unlist(purrr::map(tmp_impTypes, ~which(.x==all_impTypes)))
colorscheme <- c(RColorBrewer::brewer.pal(12,"Paired")[-c(7,8,9)][color_ind])


################################################################################
## p-value plot
pval_plot <- function(tib2plt_n_dat,top_split_var=mi_rep,legend_label){
  top_split_var <- enquo(top_split_var)
  l_f_1 <- as_labeller(function(string,name="Censoring rate: ") paste0(name,as.double(string)*100,"%"))
  global_labeller <- labeller(
    censoring_rate = l_f_1,
    .default = label_value
  )
  p <- ggplot(data=tib2plt_n_dat) +
    guides(fill=FALSE)+
    geom_density(aes(pval_r1,fill=!!top_split_var),alpha=0.5) +
    geom_density(aes(pval_r2,fill=!!top_split_var),alpha=0.5) +
    geom_density(aes(pval_r3,fill=!!top_split_var),alpha=0.5) +
    labs(color = "legend_label", y = "", x = "") +
    ggh4x::facet_nested(rows=vars(!!top_split_var),cols=vars(impType),labeller = global_labeller) +
    theme(
          panel.grid.major.y = element_line(size=0.1),
          panel.spacing = unit(0.1, "lines"),
          plot.margin = margin(t=15,r=5,b=-10,l=5),
          plot.background = element_rect(color = "grey"),
          axis.line = element_blank(),
          panel.grid = element_blank())+
    scale_x_continuous(breaks = c(0,1)) +
    scale_fill_manual(values=c("lightgrey","grey","darkgrey"))

  g <- ggplot_gtable(ggplot_build(p))
  strip_t <- stringr::str_which(g$layout$name,'strip-t-[:digit:]+')
  fills <- colorscheme
  k <- 1
  for (i in strip_t) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  g
}

################################################################################
#### regression coeff example
cens_rate <- 0.7
n_samples <- 100
tib2plt_reg_coef <- test_results_proc %>%
  filter(error_variance==0,
         mi_rep==50,
         n_dat==n_samples,
         transform_fn =="identity"
  ) %>% mutate(censoring_rate=factor(censoring_rate))

tib2plt_reg_coef_wide <- tib2plt_reg_coef %>% select(impType,censoring_rate,pval,simulation_run,sim_id) %>%
  pivot_wider(id_cols=c("impType","censoring_rate","simulation_run","sim_id"),
              values_from = "pval",
              names_from="simulation_run",
              names_prefix="pval_r")



plt_cens_pval_all <- pval_plot(tib2plt_reg_coef_wide,censoring_rate,"")
ggsave(paste0(plot_dir,"simulation_normal_effect_cens_null.png"), plt_cens_pval_all, width = plot_width, height = plot_height, units = "cm")
