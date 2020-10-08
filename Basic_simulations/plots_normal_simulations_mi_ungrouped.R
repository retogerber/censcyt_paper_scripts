# plot for single cluster simulations


library(tidyverse)
library(ggpubr)
library(ggh4x)

covariates_type <- "params_FlowCAP"

dir_save <- "./"
plot_dir <- dir_save

plot_width <- 20
plot_height <- 15


################################################################################
## b1 plot
est_plot <- function(tib2plt_n_dat, facetting_var,legend_label, breaks=NULL, param_to_plot="RB", ylab="RB",remove_top=FALSE, include_zero=TRUE){
  if (is.null(breaks)){
    breaks <- c(0,signif(seq(max(tib2plt_n_dat[[param_to_plot]]),min(tib2plt_n_dat[[param_to_plot]]),length.out = 6),2))
    print(param_to_plot)
    tmpwhich <- (which(abs(breaks[-1])< (max(breaks[-1])-min(breaks[-1]))/20))
    print(breaks[-1][tmpwhich])
    if(length(tmpwhich) != 0){
      breaks <- c(0,breaks[-1][-tmpwhich])
    }
  }
  l_f_1 <- as_labeller(function(string,name="Sample size:") paste(name,string))
  l_f_2 <- as_labeller(function(string,name="Censoring rate: ") paste0(name,as.double(string)*100,"%"))
  l_f_3 <- as_labeller(function(string,name=expression(beta["1"])) paste0(name,": ",string), default = label_parsed)
  l_f_4 <- as_labeller(function(string,name="Mulitple imputations: ") paste(name,string))
  l_f_5 <- as_labeller(function(string,name="Variance random effect: ") paste(name,string))

  global_labeller <- labeller(
    n_dat = l_f_1,
    censoring_rate = l_f_2,
    b1_True = l_f_3,
    mi_rep = l_f_4,
    variance_raneff = l_f_5,
    .default = label_value
  )
  facetting_var <- enquo(facetting_var)
  p <- ggplot(data=tib2plt_n_dat) +
    geom_rect(aes(fill=factor(!!facetting_var),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)) +
    geom_violin(aes(1, !!sym(param_to_plot),fill=impType), size=0.2,draw_quantiles = c(0.5))+
    scale_fill_manual(values=c("lightgrey","grey","darkgrey",colorscheme)) +
    labs(color = legend_label, y = ylab, x = "") +
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size=0.3),
          panel.spacing = unit(0.1, "lines"),
          legend.position="none",
          axis.ticks.x = element_blank()) +
    scale_y_continuous(breaks = breaks,labels = function(x) format(x, scientific = TRUE))
  if (include_zero){
    p <- p + geom_hline(yintercept = 0,size=0.2)
  }
  if(remove_top){
    p <- p + facet_nested(cols=vars(!!facetting_var,impType),labeller = global_labeller) +
      theme(strip.background = element_blank(),
            strip.text = element_text(colour = "white",size = 0.0001),
            plot.margin = margin(r=5,b=-10,l=5,t = -20))
    g <- ggplot_gtable(ggplot_build(p))
  } else {
    p <- p + facet_nested(cols=vars(!!facetting_var,impType),labeller = global_labeller, nest_line = TRUE) +
      theme(plot.margin = margin(r=7,b=-10,l=5))
    g <- ggplot_gtable(ggplot_build(p))

    strip_t <- stringr::str_which(g$layout$name,'strip-t-[:digit:]+-1')
    fills <- c("lightgrey","grey","darkgrey")
    k <- 1
    for (i in strip_t) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    strip_t <- stringr::str_which(g$layout$name,'strip-t-[:digit:]+-2')
    fills <- rep(colorscheme,3)
    k <- 1
    for (i in strip_t) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
  g
}


mse_plot <- function(tib2plt_n_dat, facetting_var,legend_label, breaks=NULL, param_to_plot="value", ylab="RB", remove_top=TRUE){
  if (is.null(breaks)){
    breaks <- c(0,signif(max(tib2plt_n_dat[[param_to_plot]]),2))
  }
  limits <- c(min(breaks),1.1*max(breaks))
  l_f_1 <- as_labeller(function(string,name="Sample size:") paste(name,string))
  l_f_2 <- as_labeller(function(string,name="Censoring rate: ") paste0(name,as.double(string)*100,"%"))
  l_f_3 <- as_labeller(function(string,name=expression(beta["1"])) paste0(name,": ",string), default = label_parsed)
  l_f_4 <- as_labeller(function(string,name="Mulitple imputations: ") paste(name,string))
  l_f_5 <- as_labeller(function(string,name="Variance random effect: ") paste(name,string))

  global_labeller <- labeller(
    n_dat = l_f_1,
    censoring_rate = l_f_2,
    b1_True = l_f_3,
    mi_rep = l_f_4,
    variance_raneff = l_f_5,
    .default = label_value
  )
  facetting_var <- enquo(facetting_var)
  p <- ggplot(data=tib2plt_n_dat) +
    geom_rect(aes(fill=factor(!!facetting_var),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)) +
    guides(fill=FALSE)+
    geom_col(aes(x=1,y=!!sym(param_to_plot)))+
    geom_label(aes(x=1,y=!!sym(param_to_plot),label=..y..),size=2,
               label.padding = unit(0.1, "lines"),
               label.r = unit(0.05, "lines")) +
    scale_fill_manual(values=c(rep(c("lightgrey","grey","darkgrey"),1))) +
    labs(color = legend_label, y = ylab, x = "") +
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size=0.3),
          panel.spacing = unit(0.1, "lines"),
          axis.ticks.x = element_blank(),
          legend.position="none",
          plot.margin = margin(r=5,b=-10,l=5,t = -20)
    ) +
    scale_y_continuous(breaks = breaks,limits=limits, labels = function(x) format(x, scientific = TRUE))

  if(remove_top){
    p + facet_nested(cols=vars(!!facetting_var,impType),labeller = global_labeller) +
      theme(strip.background = element_blank(),
            strip.text = element_text(colour = "white",size = 0.0001),
            plot.margin = margin(r=5,b=-10,l=5,t = -20))
  } else {
    p + facet_nested(cols=vars(!!facetting_var,impType),labeller = global_labeller, nest_line = TRUE) +
      theme(plot.margin = margin(r=7,b=-10,l=5))
  }
}

alpha <- 0.05

params_lookup = list(sample_size = "n_dat",
                       censoring_rate = "censoring_rate",
                       beta1 = "b1_True",
                       mi_repetitions = "mi_rep",
                       random_effect_variance = "variance_raneff")

sim_conditions <- c("sample_size", "censoring_rate", "beta1", "mi_repetitions", "random_effect_variance")
for(sim_cond in sim_conditions){
  tib2plt_cens <- readRDS(paste0(dir_save,"/res_sclu_output_",covariates_type,"_complete_group_",sim_cond,"_100reps.rds"))[[2]]

  tib2plt_cens <- tib2plt_cens %>% mutate(ci_lower = ifelse(impType == "cc",b1 - qt(p = 1-alpha/2,df = (n_dat-n_cens))*std_err,ci_lower),
                                          ci_upper = ifelse(impType == "cc",b1 + qt(p = 1-alpha/2,df = (n_dat-n_cens))*std_err,ci_upper))


  tib2plt_cens_ungrouped <- tib2plt_cens %>%
    mutate(RB=b1-b1_True,
           AW=ci_upper-ci_lower) %>%
    filter(impType != "km_wei") %>%
    mutate(impType = ifelse(impType == "km_exp","kme",impType)) %>%
    select(censoring_rate,impType,n_dat,transform_fn,mi_rep,variance_raneff, b1_True, RB, AW)

  tibplt <- tib2plt_cens %>%
    group_by(censoring_rate,impType,n_dat,transform_fn,mi_rep,b1_True,variance_raneff) %>%
    summarise(RB=mean(b1)-mean(b1_True),
              PB=abs((mean(b1)- mean(b1_True))/mean(b1_True)),
              CR=mean(ci_lower<mean(b1_True) & mean(b1_True) < ci_upper),
              AW=mean(ci_upper-ci_lower),
              RMSE=sqrt(mean((b1-mean(b1_True))^2)),
              b1=mean(b1))


  tibplt <- tibplt %>%
    filter(impType != "km_wei") %>%
    mutate(impType = ifelse(impType == "km_exp","kme",impType)) %>%
    pivot_longer(cols=c("RB","PB","CR","AW","RMSE")) %>%
    mutate(value=signif(value,3),
           censoring_rate=factor(censoring_rate)) %>%
    filter(name != "PB")



  tmp_impTypes <- tib2plt_cens$impType %>% unique() %>% sort()
  all_impTypes <- c("GLMM","cc","km","mrl","rs","km_exp","pmm")
  color_ind <- unlist(purrr::map(tmp_impTypes, ~which(.x==all_impTypes)))
  colorscheme <- c(RColorBrewer::brewer.pal(12,"Paired")[-c(7,8,9)][color_ind])

  param_to_plot <- sym(params_lookup[[sim_cond]])

  tib2plt_n_dat <- tib2plt_cens_ungrouped #%>% filter(!!!filt_quos) #%>%
  p1 <- est_plot(tib2plt_n_dat, !!param_to_plot,"", breaks=NULL, param_to_plot="RB", ylab=expression(hat(beta)["1"]-beta["1"]))
  p3 <- est_plot(tib2plt_n_dat, !!param_to_plot,"", breaks=NULL, param_to_plot="AW", ylab="CI width", remove_top=TRUE, include_zero = FALSE)
  tib2plt_n_dat <- tibplt %>% filter(name=="CR") #n_dat %in% n_dat_filter, censoring_rate %in% censoring_rate_filter, transform_fn=="identity")
  p2 <- mse_plot(tib2plt_n_dat, !!param_to_plot,"", breaks=NULL, param_to_plot="value", ylab="CR")
  tib2plt_n_dat <- tibplt %>% filter(name=="RMSE")#n_dat %in% n_dat_filter, censoring_rate %in% censoring_rate_filter, transform_fn=="identity")
  p4 <- mse_plot(tib2plt_n_dat, !!param_to_plot,"", breaks=NULL, param_to_plot="value", ylab="RMSE")
  plt_cens <- ggarrange(p1,p2,p3,p4,
                        common.legend = FALSE,ncol = 1,align = "v",heights = c(5,1,4,1))
  ggsave(paste0(plot_dir,"simulation_normal_mi_eval_",sim_cond,"_group.png"), plt_cens, width = plot_width, height = plot_height, units = "cm")
}

