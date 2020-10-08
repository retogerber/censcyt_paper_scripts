
saveRDS(comtests,paste0(result_dir,"da_res1_comtests_",run_type,
                        ifelse(subset_sce,"_subset",""),
                        ifelse(transform_fn,"_log",""),
                        ifelse(use_full_data,"_full",""),"_",clustering,".rds"))
#### on server:
tmpfilename <- "/home/reto/Documents/sherborne/Shared_portmacquarie/retger/FlowCAP_4/results/data/FloReMi_complete/FloReMi_Preprocessed/sce_FloReMi_complete_subset_full.rds"
# sce <- readRDS(paste0(result_dir,"/sce_", preprocessing_type, type,ifelse(subset_sce,"_subset",""), ifelse(is.null(transform_fn),"","_log"),
#                    ifelse(use_full_data,"_full",""),".rds"))
# tmpfilename <- "/home/reto/Documents/sherborne/Shared_portmacquarie/retger/FlowCAP_4/results/data/FloReMi_complete/FloReMi_Preprocessed/da_res1_mrl_FloReMi_complete_log_full_meta10.rds"
sce <- readRDS(tmpfilename)

set.seed(123)
prop_to_keep <- 0.05
len <- dim(counts(sce))[2]
cols_to_remove <- sample(seq_len(len),size = len-(round(len*prop_to_keep)) ,replace = FALSE)
sce_sub <- sce[,-cols_to_remove]


tmpfilename <- "/home/retger/censored_diffcyt/FlowCAP_4/results/data/FloReMi_complete/FloReMi_Preprocessed/da_res1_comtests_FloReMi_complete_log_full_meta100.rds"
comtests <- readRDS(tmpfilename)

# u <- filterSCE(sce, k = "meta100",
#                cluster_id %in% c(78,9,6,43,77,28,17,99,71,37))
# devtools::load_all("/home/reto/ownCloud/software/CATALYST_fork/CATALYST/")
library(CATALYST)
clustering <- "meta100"
testplt <-  plotClusterExprs(sce_sub, k = clustering, features = "type")
clusters <- levels(testplt$data$cluster_id)
clusters <- clusters[clusters!="average"]
clusters_order <- sprintf("%03i",as.integer(str_extract(clusters,"[:digit:]+"))) %>% order()
for (method_used in c("cc","rs","km","mrl")){
  clusters_to_plot <- clusters[clusters_order][comtests %>% arrange(!!sym(method_used)) %>% head() %>% pull(cluster_id) %>% as.integer()]
  data_to_plot <- testplt$data[testplt$data$cluster_id %in% c(clusters_to_plot,"average"),]
  # levels(data_to_plot$cluster_id) <- droplevels(data_to_plot$cluster_id)
  clus <- as.integer(as.integer(str_extract(levels(data_to_plot$cluster_id),"[:digit:]+")))
  clumatch <- match(clus,comtests$cluster_id)
  levels(data_to_plot$cluster_id) <- paste(levels(data_to_plot$cluster_id),
                                           "- pval",signif(comtests[clumatch,method_used],2),
                                           "- padj",signif(comtests[clumatch,paste0(method_used,"_adj")],2))
  levels(data_to_plot$cluster_id)[length(levels(data_to_plot$cluster_id))] <- "average"
  testplt_sub <- testplt %+% data_to_plot
  rm(data_to_plot)
  gc()
  ggsave(paste0(plot_dir,"plotClusterExprs_", method_used, "_",run_type,
                ifelse(subset_sce,"_subset",""),
                ifelse(use_full_data,"_full",""),"_",clustering,".png"), testplt_sub,width = 20,height = 15)
}

used_clusters <- switch(clustering,
       som400=c(138,153,158,223,133),
       meta100=c(9,98,38,53,89,94))
clusters_to_plot <- clusters[clusters_order][used_clusters]
data_to_plot <- testplt$data[testplt$data$cluster_id %in% c(clusters_to_plot,"average"),]
# levels(data_to_plot$cluster_id) <- droplevels(data_to_plot$cluster_id)
clus <- as.integer(as.integer(str_extract(levels(data_to_plot$cluster_id),"[:digit:]+")))
clumatch <- match(clus,comtests$cluster_id)
# levels(data_to_plot$cluster_id) <- paste(levels(data_to_plot$cluster_id),
#                                          "- pval",signif(comtests[clumatch,method_used],2),
#                                          "- padj",signif(comtests[clumatch,paste0(method_used,"_adj")],2))
levels(data_to_plot$cluster_id)[length(levels(data_to_plot$cluster_id))] <- "average"
testplt_sub <- testplt %+% data_to_plot
rm(data_to_plot)
gc()
ggsave(paste0(plot_dir,"plotClusterExprs_topclusters_",run_type,
              ifelse(subset_sce,"_subset",""),
              ifelse(use_full_data,"_full",""),"_",clustering,".png"), testplt_sub,width = 20,height = 15)


clco <- cluster_codes(sce_sub)
clco[clco$meta50==49,c("som400","meta100")]
clco[133,"meta100"]

# # A tibble: 10 x 4
# mrl   rs    km    cc
# <chr> <chr> <chr> <chr>
# 158   93    153   233
# 153   153   93    360
# 379   133   158   66
# 189   138   189   206
# 201   223   137   130
# 138   158   66    204
# 223   388   223   222
# 398   147   138   350
# 241   244   379   62
# 390   87    133   208


sce_filt<- filterSCE(sce, k = "meta100",
                          cluster_id %in% c(9, 98, 38, 53))
exprplot <- plotExprHeatmap(sce_filt, features = "type",
                by = "cluster_id", k = "meta100", #m = "som400",
                scale = "never", q = 0.01, perc = TRUE, col_dend = FALSE)

png(paste0(plot_dir,"plotExprHeatmap_topclusters_",run_type,
           ifelse(subset_sce,"_subset",""),
           ifelse(use_full_data,"_full",""),"_",clustering,".png"))
exprplot
dev.off()
















