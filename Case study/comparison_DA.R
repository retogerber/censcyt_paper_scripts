library(SummarizedExperiment)
library(tidyverse)
library(diffcyt)

# setwd("/home/retger/censored_diffcyt/FlowCAP_4")
setwd("/home/reto/Documents/sherborne/Shared_portmacquarie/retger/FlowCAP_4/")

main_dir <- "."
run_type <- "FloReMi_complete"
# run_type <- "complete"

use_full_data <- TRUE
subset_sce <- FALSE
transform_fn <- TRUE
# type <- "minimal_test"
# data_dir_meta <- paste0(main_dir,"/data/")
plot_dir <- paste0(main_dir,"/results/plots/",run_type,"/FloReMi_Preprocessed/")
result_dir <- paste0(main_dir,"/results/data/",run_type,"/FloReMi_Preprocessed/")
# data_out_dir <- paste0(main_dir,"/data/Spike_in/")

matching_dir <- paste0(main_dir,"/results/data/",run_type,"/Preprocessed/")

# cluster_matching <- readRDS(paste0(matching_dir, "cluster_matching.rds")) %>%
#   mutate(meta20_normal = as.character(meta20_normal),
#          meta20_signal = as.character(meta20_signal))
res_ls <- readRDS(paste0(result_dir,"assay_sce_",run_type,
                         ifelse(subset_sce,"_subset",""),
                         ifelse(use_full_data,"_full",""),".rds"))
clustering <- "meta100"
clustering_vec <- c("meta20","meta50","meta100","som400")
comtests_ls <- lapply(clustering_vec, function(clustering){

  n_cluster <- stringr::str_extract(clustering,"[:number:]+") %>% as.integer()
  comtests <- lapply(c("cc","rs","km","mrl"),function(method_used){
    da_res1_normal <- readRDS(paste0(result_dir,"da_res1_", method_used, "_",run_type,
                                     ifelse(subset_sce,"_subset",""),
                                     ifelse(transform_fn,"_log",""),
                                     ifelse(use_full_data,"_full",""),"_",clustering,".rds"))
    as.data.frame(rowData(da_res1_normal$res)) %>%
      arrange(p_val) %>%
      mutate(!!paste0(method_used,"_rank") := row_number()) %>%
      dplyr::rename(!!method_used := p_val,
             !!paste0(method_used,"_adj") := p_adj)
  }) %>%
    purrr::reduce(inner_join) %>%
    # select(-starts_with("cc")) %>%
    mutate(rank_mean=rowMeans(select(., ends_with("_rank")))) %>%
    arrange(rank_mean)


  # comtests %>% mutate(rank_median=rowMedians(select(., ends_with("_rank"))))
  # for (method_used in c("cc","rs","km","mrl")){
  #   da_res1_normal <- readRDS(paste0(result_dir,"Preprocessed/da_res1_", method_used, "_",run_type,"_som400.rds"))
  #   print(method_used)
  #   print(as.data.frame(rowData(da_res1_normal$res)) %>% arrange(p_val) %>% slice_head(n=5))
  # }



  d_counts <- res_ls[[1]]
  code_id <- rowData(d_counts)$cluster_id
  cluster_id <- metadata(d_counts)$cluster_codes[, clustering][code_id]
  rowData(d_counts)$code_id <- code_id
  rowData(d_counts)$cluster_id <- cluster_id

  form <- res_ls[[2]]
  str(form$data)
  tmptib1 <- as_tibble(assay(d_counts)) %>%
    mutate(cluster_id=cluster_id) %>%
    group_by(cluster_id) %>%
    summarise(across(starts_with("Sa"),sum))

  tmptib2 <- tmptib1 %>% ungroup() %>%
    pivot_longer(starts_with("Sa"),names_to="Sample",values_to="Abundance") %>%
    arrange(Sample) %>%
    mutate(survival_time=rep(form$data$survival_time,each=as.integer(str_extract(clustering,"[:digit:]+"))),
           Status =rep(form$data$Status,each=as.integer(str_extract(clustering,"[:digit:]+"))),
           condition =rep(form$data$condition,each=as.integer(str_extract(clustering,"[:digit:]+"))),
           patient_id = rep(form$data$patient_id,each=as.integer(str_extract(clustering,"[:digit:]+")))) %>%
    group_by(Sample) %>%
    summarise(Proportion=Abundance/sum(Abundance),cluster_id=cluster_id,survival_time=survival_time,
              Status=Status, patient_id=patient_id, condition=condition, Abundance=Abundance)

  comtests <- comtests %>%
    inner_join(tmptib2 %>% group_by(cluster_id) %>% summarise(Average_prop=mean(Proportion)))
  return(list(comtests,tmptib2))
})
names(comtests_ls) <- clustering_vec


comtests <- comtests_ls[[3]][[1]]
# tmptib <- as_tibble(t(assay(d_counts))) %>%
#   rename_with(function(x){paste0("C",x)}) %>%
#   # rowwise() %>%
#   mutate(sum=rowSums(select(., starts_with("C")))) %>%
#   mutate(survival_time=form$data$survival_time, Status =form$data$Status)
head(comtests)
ranktib <- tibble(
  mrl=head(comtests %>% arrange(mrl),10) %>% pull(cluster_id),
  rs=head(comtests %>% arrange(rs),10)%>% pull(cluster_id),
  km=head(comtests %>% arrange(km),10)%>% pull(cluster_id),
  cc=head(comtests %>% arrange(cc),10)%>% pull(cluster_id)
)
topranks <- ranktib %>% as_vector() %>% unique()
ranktib_comb <- tibble(
  cluster_id=topranks,
  mrl=topranks%in%ranktib$mrl,
  rs=topranks%in%ranktib$rs,
  km=topranks%in%ranktib$km,
  cc=topranks%in%ranktib$cc
) %>% group_by(cluster_id) %>%
  summarise(mrl=mrl,rs=rs,km=km,cc=cc,
            rank_sum=sum(mrl,rs,km,cc)) %>%
  arrange(desc(rank_sum))
ranktib_comb


methods <- c("mrl_adj","rs_adj","km_adj","cc_adj")
# methods <- c("mrl","rs","km","cc")

lapply()
padj_limit_tib <- lapply(c(0.01,0.05,0.1),function(p_level){
  tibble(
    !!methods[1] := comtests %>% arrange(!!sym(methods[1])) %>%
      mutate(!!methods[1] := !!sym(methods[1])<p_level, !!methods[1] := ifelse(is.na(!!sym(methods[1])),FALSE,!!sym(methods[1]))) %>%
      summarise(!!methods[1] := sum(!!sym(methods[1]))/n_cluster) %>% pull(!!sym(methods[1])),
    !!methods[2] := comtests %>% arrange(!!sym(methods[2])) %>%
      mutate(!!methods[2] := !!sym(methods[2])<p_level, !!methods[2] := ifelse(is.na(!!sym(methods[2])),FALSE,!!sym(methods[2]))) %>%
      summarise(!!methods[2] := sum(!!sym(methods[2]))/n_cluster) %>% pull(!!sym(methods[2])),
    !!methods[3] := comtests %>% arrange(!!sym(methods[3])) %>%
      mutate(!!methods[3] := !!sym(methods[3])<p_level, !!methods[3] := ifelse(is.na(!!sym(methods[3])),FALSE,!!sym(methods[3]))) %>%
      summarise(!!methods[3] := sum(!!sym(methods[3]))/n_cluster) %>% pull(!!sym(methods[3])),
    !!methods[4] := comtests %>% arrange(!!sym(methods[4])) %>%
      mutate(!!methods[4] := !!sym(methods[4])<p_level, !!methods[4] := ifelse(is.na(!!sym(methods[4])),FALSE,!!sym(methods[4]))) %>%
      summarise(!!methods[4] := sum(!!sym(methods[4]))/n_cluster) %>% pull(!!sym(methods[4])),
    p_level=p_level
  )
}) %>% reduce(rbind)

padj_limit_tib



# comtests %>% filter(cluster_id == 9)
# comtests %>% filter(cluster_id == 98)
#
#
# used_clusters <- switch(clustering,
#                         som400=c(138,153,158,223,133),
#                         meta100=c(9,98,38,53,89,94),
#                         meta50=c(12,43,49,7))
# comtests %>% filter(cluster_id %in% used_clusters)
#
#
#
#
#
#
#
# padjtib <- tibble(
#   mrl=head(comtests %>% arrange(mrl),10) %>% pull(mrl_adj),
#   rs=head(comtests %>% arrange(rs),10)%>% pull(rs_adj),
#   km=head(comtests %>% arrange(km),10)%>% pull(km_adj),
#   cc=head(comtests %>% arrange(cc),10)%>% pull(cc_adj)
# )
#
#
#
#
# tmpwide <- pivot_wider(comtests_ls[[3]][[2]] %>% filter(cluster_id==9) %>% ungroup() %>% select(-Sample),
#                        names_from = "condition", values_from = "Proportion") %>%
#   mutate(diff=Stim-Unstim)
#
# ggplot(tmpwide ) +
#   # geom_point(aes(survival_time,Proportion,color=condition)) +
#   geom_point(aes(survival_time,diff,color=Status))
#
#
# ggplot(comtests_ls[[3]][[2]] %>% filter(cluster_id==9) ) +
#   # geom_point(aes(survival_time,Proportion,color=condition)) +
#   geom_point(aes(survival_time,Proportion,color=Status)) +
#   geom_line(aes(x=survival_time,y=zoo::rollmedian(Proportion, 21, na.pad=TRUE)))
#   # facet_wrap(~condition)

l_f_1 <- as_labeller(function(string,name="Cluster number:") paste(name,string))

global_labeller <- labeller(
  cluster_id = l_f_1,
  .default = label_value
)
# zoo::rollmedian(comtests_ls[[3]][[2]] %>% filter(cluster_id==9) %>% pull(Proportion), 21, na.pad=TRUE)

# tmppltdat <- comtests_ls[[3]][[2]] %>%
#   arrange(survival_time) %>%
#   filter(cluster_id %in% c(9)) %>%
#   ungroup() %>%
#   mutate(Abundance =  comtests_ls[[3]][[2]] %>% group_by(Sample) %>% summarise(Abundance=sum(Abundance)) %>% pull(Abundance)) %>%
#   filter(Status==1)
# # ggplot(tmppltdat ) +
# #   geom_point(aes(survival_time,Proportion,color=Status,shape=condition)) +
# #   ylim(c(0,0.1)) +
# #   geom_smooth(aes(x=survival_time,y=Proportion),formula = y~x+condition,span=0.8, method="glm", method.args = list(family="binomial"))+
# #   # geom_line(aes(x=exp(survival_time),y=zoo::rollmedian(Proportion, 51, na.pad=TRUE))) +
# #   facet_wrap(~cluster_id,scales = "free_y",labeller = global_labeller)
#
#
#
#
# glmerfit <- lme4::glmer(Proportion ~ survival_time + condition + (1|Sample),
#                         data = tmppltdat,
#                         family = "binomial",
#                         weights = tmppltdat$Abundance)
#
# glmerfitsum <- summary(glmerfit)
# coef(glmerfitsum)
#
# glmerfun <- function(x,z){
#   dlogis(coef(glmerfitsum)[1,1]+coef(glmerfitsum)[2,1]*x+coef(glmerfitsum)[3,1]*z)
# }
# ggplot(tmppltdat ) +
#   geom_point(aes(survival_time,Proportion,color=Status,shape=condition)) +
#   ylim(c(0,0.1)) +
#   geom_function(fun = glmerfun, colour = "red", args=list(z=condition))+
#   # geom_smooth(aes(x=survival_time,y=Proportion),formula = y~x+condition,span=0.8, method="glm", method.args = list(family="binomial"))+
#   # geom_line(aes(x=exp(survival_time),y=zoo::rollmedian(Proportion, 51, na.pad=TRUE))) +
#   facet_wrap(~cluster_id,scales = "free_y",labeller = global_labeller)




comred <- comtests_ls[[3]][[2]] %>%
  arrange(survival_time) %>%
  group_by(cluster_id) %>%
  mutate(mean_Prop=mean(Proportion)) %>%
  filter(cluster_id %in% c(38,71)) %>%
  # filter(Proportion/mean_Prop < 10) %>%
  mutate(Status=factor(Status,levels=c(0,1),labels=c("censored","observed")),
         is_DA = if_else(cluster_id %in% c(9,38),"DA","non DA"),
         log_prop=log(Proportion+0.00001),
         exp_surv=exp(survival_time)-11)

out <- conditional_multiple_imputation(comred %>% filter(cluster_id==38) %>%
                                  mutate(Status=as.integer(Status)-1,
                                         patient_id=as.character(patient_id),
                                         condition=as.integer(condition)-1),
                                formula(Proportion~Surv(survival_time,Status)+condition+(1|Sample)+(1|patient_id)),
                                regression_type = "glmer",mi_reps = 20,imputation_method = "rs")
medparams <- colMedians(out[[4]]$betas)

comred <- comred %>% mutate(glmerfit = dlogis(medparams[1]+medparams[2]*survival_time))

p <- ggplot(comred %>% arrange(survival_time) %>% filter(cluster_id==38)) +
  # geom_point(aes(survival_time,Proportion,color=Status,shape=condition)) +
  geom_line(aes(survival_time,glmerfit))+
  # ylim(c(0,0.02)) +
  # coord_trans(y = "log") +
  # geom_smooth(aes(x=survival_time,y=Proportion), method="lm")+
  # geom_line(aes(x=exp(survival_time),y=zoo::rollmedian(Proportion, 51, na.pad=TRUE))) +
  # facet_grid(rows=vars(is_DA,cluster_id))
  scale_color_brewer(palette="Set1",direction=-1) +
  facet_wrap(is_DA~cluster_id,labeller = global_labeller,ncol = 1) +
  labs(shape="",color="",x="log(Survival time + 11)") +
  theme_bw()
p
ggsave(filename = paste0(plot_dir,"example_proportionvssurvival_plot.png"),plot = p, width = 7,height = 9)



comredsum <- comtests_ls[[3]][[2]] %>%
  inner_join(comtests_ls[[3]][[1]],by="cluster_id") %>%
  arrange(survival_time) %>%
  group_by(cluster_id) %>%
  mutate(mean_Prop=mean(Proportion)) #%>%
  # filter(cluster_id %in% c(9,38,68,71)) %>%
  # filter(Proportion/mean_Prop < 10)


comredsum_wide <- comredsum %>% pivot_longer(cols=c("rs","rs_adj","cc","cc_adj","mrl","mrl_adj","km","km_adj"))

ggplot(comredsum_wide %>% select(name,value,mean_Prop) %>% unique() #%>% filter(stringr::str_detect(name,"_adj"))
       ) +
  geom_density(aes(x=value)) +
  facet_wrap(~factor(name),scales = "free_y")


ranpval <- tibble(id=1:100,x=runif(100),x_adj=p.adjust(x)) %>% pivot_longer(cols = c("x","x_adj"))
combpvaldat <- rbind(comredsum_wide %>% select(name,value) %>% unique(),ranpval)
ggplot(combpvaldat) +
  geom_density(aes(x=value))+
  facet_wrap(~name, scales = "free_y")



# geom_density(aes(x=Proportion,color=Status))
# geom_density(aes(x=survival_time))
# geom_line(aes(x=survival_time,y=zoo::rollmedian(Proportion, 21, na.pad=TRUE)))

# zoo::rollmean(tmptib2 %>% filter(cluster_id==9) %>% pull(Proportion), 100, na.pad=TRUE)
#
#
# ggplot(tmptib2 %>% filter(cluster_id==9) ) +
#   geom_point(aes(survival_time,Proportion,color=Status)) +
#   geom_density(aes(x=survival_time))

# tmpfilename <- "/home/reto/Documents/sherborne/Shared_portmacquarie/retger/FlowCAP_4/results/data/FloReMi_complete/FloReMi_Preprocessed/sce_FloReMi_complete_subset_full.rds"
# sce <- readRDS(tmpfilename)
#
# topclus <- comtests_ls$meta100 %>% arrange(mrl) %>% head(12) %>% pull(cluster_id) %>% as.integer()
#
# cluco <- cluster_codes(sce)
# m100 <- lapply(topclus, function(x) which(cluco$meta100==x))
#
# m50_of_m100 <- sapply(m100, function(x) unique(cluco$meta50[x]))
#
# unlist(lapply(m100,function(x) length(x))) / sapply(m50_of_m100, function(x) sum(cluco$meta50==x))
#
#
# sum(cluco$meta100==100)
# m50_of_m100 <- cluco$meta50[topclus]
#
# sapply(unique(m50_of_m100), function(x) sum(cluco$meta50==x))
