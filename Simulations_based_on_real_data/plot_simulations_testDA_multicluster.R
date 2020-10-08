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
slope_selected <- "0.9,0.9,0.9"
group_slope_selected <- "0.2,0.2,0.2"
n_sam_selected <- 50
censrate_selected <- 0.5
# transform_fn <- "log_positive"
# transform_fn <- "div_100"
transform_fn <- ""

# transformation_scale_x <- "log10"
# transformation_scale_x <- "identity"
transformation_scale_x <- "sqrt"


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


## load data
savelist <- readRDS(paste0(plot_dir,"multicluster_cobra_data_af_",alpha_factor,"_",clustering,"_",sim_type,ifelse(transform_fn=="","",paste0("_",transform_fn)),".rds"))
plts <- savelist$cobradata
arg_df <- savelist$arg_df


# ################################################################################
# ## censoring
# if (is_group_bool){
#   tmp_plt_ind <- arg_df %>%
#     filter(n_sam==n_sam_selected,slope_1==slope_selected, group_slope_1==group_slope_selected)
# }else{
#   tmp_plt_ind <- arg_df %>%
#     filter(n_sam==n_sam_selected,slope_1==slope_selected)
# }
#
# tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
#   p <- suppressMessages(plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1,stripsize = 0)  +
#     scale_x_continuous(breaks = c(0.01,0.1,1),limits = c(0,1),trans = transformation_scale_x) +
#     labs(title=paste0("censoring rate: ",100*tmp_plt_ind$censrate[i],"%"))  +
#       theme(aspect.ratio=1,
#             strip.background = element_blank(),
#             strip.text = element_blank())+
#     guides(colour = guide_legend(nrow = 1)))
#   if (i!=1){
#     p + theme(axis.text.y = element_blank(),
#               axis.title.y = element_blank())
#   } else {
#     p + theme(plot.margin = margin(t = 5,r = 0,b = 5,l = 0))
#   }
#
# })
# plt_comb <- ggarrange(plotlist = tmp_plts,nrow=1,ncol=3, common.legend = TRUE,widths = c(5,4,4))#,ncol=2)
# bottom_annot <- paste("Size:",unique(tmp_plt_ind$n_sam),"\t",
#                       "Slope:",stringr::str_split(unique(tmp_plt_ind$slope_1),",")[[1]][1])
# if (is_group_bool){
#   bottom_annot <- paste(bottom_annot, "\t","Group slope:",stringr::str_split(unique(tmp_plt_ind$group_slope_1),",")[[1]][1])
# }
# plt_comb <- annotate_figure(plt_comb,
#                             fig.lab = "Censoring rate",fig.lab.size = 14,fig.lab.face = "bold",
#                             bottom = bottom_annot)
# plt_comb
# ggsave(paste0(plot_dir,"simulation_multicluster_effect_cens_af_",alpha_factor,"_",clustering,"_",sim_type,".png"), plt_comb, width = 30, height = 13, units = "cm")
#
# ################################################################################
# ## size
# if (is_group_bool){
#   tmp_plt_ind <- arg_df %>%
#     filter(censrate==censrate_selected,slope_1==slope_selected, group_slope_1==group_slope_selected)%>%
#     arrange(desc(n_sam))
# }else{
#   tmp_plt_ind <- arg_df %>%
#     filter(censrate==censrate_selected,slope_1==slope_selected)%>%
#     arrange(desc(n_sam))
# }
# tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
#   p <- suppressMessages(plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  +
#     scale_x_continuous(breaks = c(0.01,0.1,1),limits = c(0,1),trans = transformation_scale_x) +
#     labs(title=paste0("Size: ",tmp_plt_ind$n_sam[i])) +
#     theme(aspect.ratio=1,
#           strip.background = element_blank(),
#           strip.text = element_blank()) +
#     guides(colour = guide_legend(nrow = 1)))
#   # if (i!=1){
#   #   p + theme(axis.text.y = element_blank(),
#   #             axis.title.y = element_blank())
#   # } else {
#   #   p + theme(plot.margin = margin(t = 5,r = 0,b = 5,l = 0))
#   # }
#   p
# })
# plt_comb <- ggarrange(plotlist = tmp_plts,nrow=2,ncol=2, common.legend = TRUE)#,ncol=2)
# bottom_annot <- paste("Censoring rate:",unique(tmp_plt_ind$censrate),"\t",
#                       "Slope:",stringr::str_split(unique(tmp_plt_ind$slope_1),",")[[1]][1])
# if (is_group_bool){
#   bottom_annot <- paste(bottom_annot, "\t","Group slope:",stringr::str_split(unique(tmp_plt_ind$group_slope_1),",")[[1]][1])
# }
# plt_comb <- annotate_figure(plt_comb,
#                             fig.lab = "Sample size",fig.lab.size = 14,fig.lab.face = "bold",
#                             bottom = bottom_annot)
# ggsave(paste0(plot_dir,"simulation_multicluster_effect_size_af_",alpha_factor,"_",clustering,"_",sim_type,".png"), plt_comb, width = 30, height = 33, units = "cm")
#
# ################################################################################
# ## slope
# if (!is_group_bool){
#   tmp_plt_ind <- arg_df %>%
#     filter(censrate==censrate_selected,n_sam==n_sam_selected) %>%
#     arrange(desc(n_sam))
#
#   tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
#     suppressMessages(plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  +
#       scale_x_continuous(breaks = c(0.005,0.01,0.05,0.1,1),limits = c(0.005,1),trans = transformation_scale_x) +
#       labs(title=paste0("slope: ",stringr::str_split(tmp_plt_ind$slope[i],",")[[1]][1])) +
#       theme(aspect.ratio=1) +
#       guides(colour = guide_legend(nrow = 1)))
#   })
#   plt_comb <- ggarrange(plotlist = tmp_plts,nrow=1,ncol=3, common.legend = TRUE)#,ncol=2)
#   bottom_annot <- paste("Size:",unique(tmp_plt_ind$n_sam),"\t",
#                         "Censoring rate:",unique(tmp_plt_ind$censrate),"\t")
#   plt_comb <- annotate_figure(plt_comb,
#                               fig.lab = "Slope",fig.lab.size = 14,fig.lab.face = "bold",
#                               bottom = bottom_annot)
#   ggsave(paste0(plot_dir,"simulation_multicluster_effect_slope_af_",alpha_factor,"_",clustering,"_",sim_type,".png"), plt_comb, width = 30, height = 13, units = "cm")
# }
# ################################################################################
# ## group slope
# if (is_group_bool){
#   tmp_plt_ind <- arg_df %>%
#     filter(censrate==censrate_selected,n_sam==n_sam_selected) %>%
#     arrange(desc(n_sam))
#   tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
#     suppressMessages(plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  +
#                        scale_x_continuous(breaks = c(0.005,0.01,0.05,0.1,1),limits = c(0.005,1),trans = transformation_scale_x) +
#                        labs(title=paste0("group slope: ",stringr::str_split(tmp_plt_ind$group_slope_1[i],",")[[1]][1])) +
#                        theme(aspect.ratio=1) +
#                        guides(colour = guide_legend(nrow = 1)))
#   })
#   plt_comb <- ggarrange(plotlist = tmp_plts,nrow=1,ncol=3, common.legend = TRUE)#,ncol=2)
#   bottom_annot <- paste("Size:",unique(tmp_plt_ind$n_sam),"\t",
#                         "Censoring rate:",unique(tmp_plt_ind$censrate),"\t",
#                         "Slope:",stringr::str_split(unique(tmp_plt_ind$slope_1),",")[[1]][1])
#   plt_comb <- annotate_figure(plt_comb,
#                               fig.lab = "Group slope",fig.lab.size = 14,fig.lab.face = "bold",
#                               bottom = bottom_annot)
#   ggsave(paste0(plot_dir,"simulation_multicluster_effect_group_slope_af_",alpha_factor,"_",clustering,"_",sim_type,".png"), plt_comb, width = 30, height = 13, units = "cm")
# }

################################################################################
## censoring rate and size
if (is_group_bool){
  tmp_plt_ind <- arg_df %>%
    filter(slope_1==slope_selected, group_slope_1==group_slope_selected)%>%
    arrange(desc(censrate))
}else{
  tmp_plt_ind <- arg_df %>%
    filter(slope_1==slope_selected)%>%
    arrange(desc(censrate))
}
tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
  p <- suppressMessages(plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  +
                          coord_cartesian(ylim=c(0,1),xlim=c(0,0.6))+
                          # scale_y_continuous(limits=c(0.4,1)) +
                          scale_x_continuous(breaks = c(0.01,0.1,0.6),trans = transformation_scale_x) +
                          labs(title=paste0("Size: ",tmp_plt_ind$n_sam[i],"\tcensoring rate: ",100*tmp_plt_ind$censrate[i],"%")) +
                     theme(#aspect.ratio=1,
                           strip.background = element_blank(),legend.position = "bottom",
                           strip.text = element_blank(),legend.text = element_text(size=20),legend.key = element_rect(size=50),
                           axis.title.x = element_blank(),
                           axis.title.y = element_blank()
                           )+
                     guides(colour = guide_legend(nrow = 1,override.aes = list(size=3))))
  if (i%%4!=1){
    p <- p + theme(axis.text.y = element_blank())
  }
  if (i<9){
    p <- p + theme(axis.text.x = element_blank())
  }
  # else {
  #   p <- p + theme(plot.margin = margin(t = 5,r = 0,b = 5,l = 0))
  # }
  p

})
plt_comb <- ggarrange(plotlist = tmp_plts,nrow=3,ncol=4, common.legend = TRUE,widths = c(10,9,9,9),heights = c(9,9,10),legend = "top")
# plt_comb <- ggarrange(plotlist = tmp_plts,nrow=3,ncol=4, common.legend = TRUE)#,ncol=2)
# bottom_annot <- paste("Slope:",stringr::str_split(unique(tmp_plt_ind$slope_1),",")[[1]][1])
# if (is_group_bool){
#   bottom_annot <- paste(bottom_annot, "\t","Group slope:",stringr::str_split(unique(tmp_plt_ind$group_slope_1),",")[[1]][1])
# }
plt_comb <- annotate_figure(plt_comb,left = text_grob("TPR", face = "bold", size = 20,rot = 90),
                            # fig.lab = "Size and Censoring rate",fig.lab.size = 14,fig.lab.face = "bold",
                            bottom = text_grob("FDR", face = "bold", size = 20))
# plt_comb
ggsave(paste0(plot_dir,"simulation_multicluster_effect_size_cens_af_",alpha_factor,"_",clustering,"_",sim_type,"_",transform_fn,".png"), plt_comb, width = 50*0.7, height = 40*0.7, units = "cm")

# ################################################################################
# ## all combined
# tmp_plt_ind <- arg_df %>%
#   arrange(desc(n_sam))
# if (is_group_bool){
#   tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
#     suppressMessages(plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  +
#                        scale_x_continuous(breaks = c(0.005,0.01,0.05,0.1,1),limits = c(0.005,1),trans = transformation_scale_x) +
#                        labs(title=paste0("censoring rate: ",100*tmp_plt_ind$censrate[i],"%\n",
#                                          "Group slope: ",stringr::str_split(tmp_plt_ind$group_slope[i],",")[[1]][1],"\n",
#                                          "Size: ",tmp_plt_ind$n_sam[i])) +
#                        theme(aspect.ratio=1) +
#                        guides(colour = guide_legend(nrow = 1)))
#   })
#   bottom_annot <- paste("Slope:",stringr::str_split(unique(tmp_plt_ind$slope_1),",")[[1]][1])
#
# }else{
#   tmp_plts <- purrr::map(seq_along(tmp_plt_ind$id), function(i){
#     suppressMessages(plot_fdrtprcurve(plts[tmp_plt_ind$id[i]][[1]],pointsize = 2,linewidth = 1)  +
#                        scale_x_continuous(breaks = c(0.005,0.01,0.05,0.1,1),limits = c(0.005,1),trans = transformation_scale_x) +
#                        labs(title=paste0("censoring rate: ",100*tmp_plt_ind$censrate[i],"%\n",
#                                          "slope: ",stringr::str_split(tmp_plt_ind$slope[i],",")[[1]][1],"\n",
#                                          "Size: ",tmp_plt_ind$n_sam[i])) +
#                        theme(aspect.ratio=1) +
#                        guides(colour = guide_legend(nrow = 1)))
#   })
#   bottom_annot <- ""
#
# }
# plt_comb <- ggarrange(plotlist = tmp_plts,nrow=3,ncol=12, common.legend = TRUE)#,ncol=2)
# plt_comb <- annotate_figure(plt_comb,fig.lab.size = 14,fig.lab.face = "bold",
#                             bottom = bottom_annot)
# ggsave(paste0(plot_dir,"simulation_multicluster_all_conditions_af_",alpha_factor,clustering,"_",sim_type,".png"), plt_comb, width = 100, height = 30, units = "cm")
#
#
#
