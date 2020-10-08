# script to obtain weibull distribution from FlowCAP IV dataset


library(tidyverse)
library(survival)

# where the metadata file is located
data_dir_meta <- "./FlowCAP_4/data/"
# where to save the data and the plots
plot_dir <- "./FlowCAP_4/results/"

MetaDataFile <- paste0(data_dir_meta, "MetaDataTrain.csv")
MetaData <- as_tibble(read.csv(MetaDataFile, sep = ",", header = TRUE)) %>%
  dplyr::filter(!is.na(Status)) %>%
  dplyr::select(-Prediction) %>%
  mutate(patient_id = as.factor(seq(1:length(Status)))) %>%
  dplyr::rename(survival_time = "Survival.Time")

cat("Preprocessing:\n")
MetaData_melt <- reshape2::melt(MetaData, id.vars = c("Status","survival_time","patient_id"), value.name = "sample_name", variable.name = "condition") %>%
  arrange(survival_time, patient_id,condition)

# remove all censored values higher than the highest observed one for fitting
maxobs <- max(which(MetaData_melt$Status==1))
MetaData_melt_fit <- MetaData_melt[seq_len(maxobs),]

# data format for fitting in fitdistrplus
in_fit_cens <- data.frame(left=MetaData_melt_fit$survival_time,right=ifelse(MetaData_melt_fit$Status==1,MetaData_melt_fit$survival_time,NA))
weifit <- fitdistrplus::fitdistcens(in_fit_cens[in_fit_cens$left>0,],"weibull")


## fit separately
# censored:
MetaData_melt_fit_cens <- MetaData_melt[MetaData_melt$Status==0,]
weifit_cens <- fitdistrplus::fitdist(MetaData_melt_fit_cens[MetaData_melt_fit_cens$survival_time>0,"survival_time"],
                                     "weibull")
weifit_cens
ggplot(MetaData_melt_fit_cens) +
  geom_histogram(aes(survival_time,after_stat(density)),fill="darkgreen") +
  stat_function(fun=dweibull,args=list(scale=weifit_cens$estimate[2],shape=weifit_cens$estimate[1]))
# observed:
MetaData_melt_fit_obs <- MetaData_melt[MetaData_melt$Status==1,]
weifit_obs <- fitdistrplus::fitdist(MetaData_melt_fit_cens[MetaData_melt_fit_obs$survival_time>0,"survival_time"],"weibull")
weifit_obs
ggplot(MetaData_melt_fit_obs) +
  geom_histogram(aes(survival_time,after_stat(density)),fill="darkgreen") +
  stat_function(fun=dweibull,args=list(scale=weifit_obs$estimate[2],shape=weifit_obs$estimate[1]))


## compare two censoring mechanisms
# uniform censoring
n <- dim(MetaData_melt_fit)[1]
cens_ind <- which(MetaData_melt_fit$Status==0)
T <- rweibull(n,weifit$estimate[1],weifit$estimate[2])
T_min <- 0
T_actual <- T
T_actual[cens_ind]  <- sapply(cens_ind,function(i) runif(1,T_min,T[i]))

# kolmogorov smirnov test
ks_test_fit_comp <- ks.test(T_actual,MetaData_melt_fit$survival_time)
ks_test_fit_cens <- ks.test(T_actual[cens_ind],MetaData_melt_fit$survival_time[cens_ind])
ks_test_fit_obs <- ks.test(T_actual[-cens_ind],MetaData_melt_fit$survival_time[-cens_ind])

tmpdf <- cbind(MetaData_melt_fit,T_actual) %>%
  dplyr::mutate( ks_pval=ifelse(Status==1,ks_test_fit_obs$p.value,ks_test_fit_cens$p.value),
                 Status=factor(Status,levels = c(0,1),labels = c("censored","observed")))
tmpdf <- rbind(tmpdf,tmpdf %>% mutate(Status="combined",ks_pval=ks_test_fit_comp$p.value))
plt_unif <- ggplot(tmpdf) +
  geom_histogram(aes(T_actual,after_stat(density),fill="simulated"),alpha=0.5) +
  geom_histogram(aes(survival_time,after_stat(density),fill="real"),alpha=0.5) +
  facet_grid(col=vars(Status)) +
  labs(fill="",x="") +
  xlim(c(0,6000)) +
  ylim(c(0,0.00125)) +
  geom_text(aes(x=2000,y=0.0012,label=paste0("ks pvalue: ",signif(ks_pval,3))))

# parametric censoring
T <- rweibull(n,weifit$estimate[1],weifit$estimate[2])
C <- rweibull(n,weifit_cens$estimate[1],weifit_cens$estimate[2])
sum(T>C)/n
T_actual <- ifelse(T>C,C,T)

# kolmogorov smirnov test
ks_test_fit_comp <- ks.test(T_actual,MetaData_melt_fit$survival_time)
ks_test_fit_cens <- ks.test(T_actual[cens_ind],MetaData_melt_fit$survival_time[cens_ind])
ks_test_fit_obs <- ks.test(T_actual[-cens_ind],MetaData_melt_fit$survival_time[-cens_ind])

tmpdf <- cbind(MetaData_melt_fit,T_actual) %>%
  dplyr::mutate( ks_pval=ifelse(Status==1,ks_test_fit_obs$p.value,ks_test_fit_cens$p.value),
                 Status=factor(Status,levels = c(0,1),labels = c("censored","observed")))
tmpdf <- rbind(tmpdf,tmpdf %>% mutate(Status="combined",ks_pval=ks_test_fit_comp$p.value))
plt_param <- ggplot(tmpdf) +
  geom_histogram(aes(T_actual,after_stat(density),fill="simulated"),alpha=0.5) +
  geom_histogram(aes(survival_time,after_stat(density),fill="real"),alpha=0.5) +
  facet_grid(col=vars(Status))+
  labs(fill="",x="survival time") +
  xlim(c(0,6000)) +
  ylim(c(0,0.00125)) +
  geom_text(aes(x=2000,y=0.0012,label=paste0("ks pvalue: ",signif(ks_pval,3))))

# all parameters of the fitting
weifit_ls <- list(weifit=weifit,weifit_cens=weifit_cens,weifit_obs=weifit_obs)
saveRDS(weifit_ls,paste0(plot_dir,"weibull_fits_FlowCAP.rds"))

# comparison plot of different censoring mechanisms
plt_comb <- ggpubr::ggarrange(plt_unif,plt_param,ncol = 1,common.legend = TRUE,labels = c("uniform","parametric"),label.y = 1.05,label.x = -0.01)
ggsave(paste0(plot_dir,"censoring_mechansim_comparison.png"),plt_comb,width = 12,height = 8)













