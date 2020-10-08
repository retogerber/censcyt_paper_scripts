# run complete case study

# command line input
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Exactly one argument must be supplied, '1' for normal data, '2' for spike-in data", call.=FALSE)
}
# iii = 1 if preprocessed data, iii = 2 if spike_in data
iii <- as.integer(args[1])

# set working directory to directory of this script
# setwd(directory_name)

# either "complete" for using all samples or "test_run" for a subset
run_type <- "complete"

# either "" for the standard preprocessing or "FloReMi_" for the processing from
# the FloReMi workflow (main difference is removal of potentially dead cells with Gating
# in FloReMi)
preprocessing_type <- "FloReMi_"
# preprocessing_type <- ""

# random subset of max 10000 cells per sample
subset_sce <- FALSE

# transformation of survival time, NULL if none
transform_fn <- function(x) {log(x+11)}
# transform_fn <- NULL

# remove samples with fewer than 10000 cells (only if subset_sce==TRUE)
filter_small_samples <- FALSE

# only Stim
second_cov <- TRUE

# full_data, (instead of only training data from flowCAP4)
use_full_data <- TRUE

set.seed(123)
################################################################################
## load packages

suppressPackageStartupMessages({
  library(tidyverse)
  library(diffcyt)
  library(CATALYST)})
RhpcBLASctl::blas_set_num_threads(1)
t_start <- Sys.time()
# functions for printing system time difference
print_total_time <- function(){
  cat("Total time:",round(as.double(Sys.time()-t_start),2), units(Sys.time()-t_start),"\n")
}
print_partial_time <- function(){
  cat("Completed in",round(as.double(Sys.time()-t_tmp),2), units(Sys.time()-t_tmp),"\n")
}

main_dir <- "."
data_type <- paste0(preprocessing_type, c("Preprocessed","Spike_in"))
if (preprocessing_type=="FloReMi_"){
  filename_pattern <- c("_preprocessed","_frm_signal")
}else{
  filename_pattern <- c("_pp","_signal")
}
if (iii == 1){
  type <- run_type
}else{
  type <- paste0(run_type,"_spike_in")
}
data_dir <- paste0(main_dir,"/data/", data_type[iii], "/")
data_dir_meta <- paste0(main_dir,"/data/")
if (!dir.exists(paste0(main_dir,"/results/plots/",preprocessing_type,run_type))){
  dir.create(paste0(main_dir,"/results/plots/",preprocessing_type,run_type))
}
plot_dir <- paste0(main_dir,"/results/plots/",preprocessing_type,run_type,"/",data_type[iii],"/")
if (!dir.exists(plot_dir)){
  dir.create(plot_dir)
}
if (!dir.exists(paste0(main_dir,"/results/data/",preprocessing_type,run_type))){
  dir.create(paste0(main_dir,"/results/data/",preprocessing_type,run_type))
}
result_dir <- paste0(main_dir,"/results/data/",preprocessing_type,run_type,"/",data_type[iii],"/")
if (!dir.exists(result_dir)){
  dir.create(result_dir)
}
if(use_full_data){
  MetaDataFile <- paste0(data_dir_meta, "MetaDataFull.csv")
} else {
  MetaDataFile <- paste0(data_dir_meta, "MetaDataTrain.csv")
}
ChannelsFile <- paste0(data_dir_meta, "FlowCAPchannels.csv")

MetaData <- as_tibble(read.csv(MetaDataFile, sep = ",", header = TRUE)) %>%
  dplyr::filter(!is.na(Status)) %>%
  dplyr::select(Status,Survival.Time,Stim, Unstim) %>%
  mutate(patient_id = as.factor(seq(1:length(Status)))) %>%
  dplyr::rename(survival_time = "Survival.Time")

if(!is.null(transform_fn)){
  MetaData <- MetaData %>% dplyr::mutate(survival_time = transform_fn(survival_time))
}


cat("Meta data preprocessing:\n")
MetaData_melt <- reshape2::melt(MetaData, id.vars = c("Status","survival_time","patient_id"), value.name = "sample_name", variable.name = "condition")
# sort according to time
MetaData_melt <- MetaData_melt %>%
  arrange(survival_time, patient_id,condition)
experiment_info <- MetaData_melt %>%
  mutate(sample_id = as.factor(paste0("Sa_", gsub(".fcs","",sample_name)))) %>%
  dplyr::select(-sample_name) #%>%

if (!second_cov){
  experiment_info <- experiment_info[experiment_info$condition=="Stim",]
}

if (run_type == "test_run"){
  experiment_info <- experiment_info[seq(1,380,10),]
}

Channels <- as_tibble(read.csv(ChannelsFile, sep = ",", header = TRUE, stringsAsFactors = FALSE))
marker_info <- Channels %>%
  mutate(marker_class = if_else((X == "Not used") | (Channel.Name %in% c("Time","FSC-A","FSC-H","SSC-A")), "none","type")) %>%
  transmute(channel_name = Channel.Name, marker_name = if_else(Reagent != "", Reagent, channel_name), marker_class) #%>%


# subset of processed files
files <- experiment_info %>%
  mutate(rowid = row_number()) %>%
  dplyr::select(sample_id) %>%
  transmute(filename = paste0(data_dir,as.character(gsub("Sa_","",sample_id)), filename_pattern[iii], ".fcs")) %>%
  as_vector()

stopifnot(all(basename(files) %in% list.files(data_dir)))

cols_to_include <- marker_info %>% transmute(marker_class_b = marker_class == "type") %>% as_vector()

print_total_time()
cat("start reading flowSet\n")
t_tmp <- Sys.time()

if (subset_sce){
  d_flowframes <- lapply(files, function(file){
    flowCore::read.FCS(file,
                       transformation = FALSE,
                       truncate_max_range = FALSE)
  })
  d_flowSet <- flowCore::flowSet(d_flowframes)
  rm(d_flowframes)
  if (!filter_small_samples){
  d_flowSet_filt <- lapply(seq_along(d_flowSet), function(i) {
    if(dim(d_flowSet[[i]])[1]<10000){
      d_flowSet[[i]]
    } else {
      d_flowSet[[i]][sample(dim(d_flowSet[[i]])[1],replace=FALSE,size=10000),]
    }
    })
  } else if (filter_small_samples){
    d_flowSet_filt <- lapply(seq_along(d_flowSet), function(i) {
      if(dim(d_flowSet[[i]])[1]<10000){
        NULL
      } else {
        d_flowSet[[i]][sample(dim(d_flowSet[[i]])[1],replace=FALSE,size=10000),]
      }
    })
  }
  d_flowSet_is_not_null <- which(!sapply(d_flowSet_filt, is.null))
  d_flowSet <- flowCore::flowSet(d_flowSet_filt[!sapply(d_flowSet_filt,is.null)])

} else {
  d_flowSet <- flowCore::read.flowSet(
    files, transformation = FALSE, truncate_max_range = FALSE, which.lines = unlist(ifelse(subset_sce,10000,list(NULL)))
  )
}

# subset of experiment_info
experiment_info <- experiment_info %>%
  mutate(rowid = row_number()) %>%
  dplyr::select(-rowid) %>%
  mutate(patient_id = droplevels(as.factor(patient_id)),
         sample_id = droplevels(as.factor(sample_id)),
         condition = droplevels(as.factor(condition)))

if(exists("d_flowSet_is_not_null")){
  experiment_info <- experiment_info[d_flowSet_is_not_null,]
}



if (data_type[iii] == "Spike_in"){
  metadata_frames <- flowCore::fsApply(d_flowSet,function(x){
    exprs(x)[ ,"cell_id", drop = FALSE]
  })
}


marker_info <- marker_info %>% dplyr::rename(channel = channel_name, antigen = marker_name, class = marker_class)
md <- data.frame(experiment_info)
md$file_name <- gsub("Sa_","",paste0(md$sample_id, filename_pattern[iii],".fcs"))

cat("prepData\n")
t_tmp <- Sys.time()
colnames(d_flowSet) <- gsub("\\.","-",colnames(d_flowSet))
sce <- prepData(x = d_flowSet,
                 panel = data.frame(marker_info),
                 md = md,
                 panel_cols = list(channel = "channel",
                                   antigen = "antigen",
                                   class = "class"),
                 md_cols = list(file = "file_name",
                                id = "sample_id",
                                factors = c("condition","survival_time","Status","patient_id")),
                 transform = FALSE)
assay(sce,"exprs") <- assay(sce,"counts")
print_partial_time()
rm(d_flowSet)
if (run_type == "test_run" | subset_sce){

  cat("plotting 5 plots\n")
  t_tmp <- Sys.time()
  tryCatch({
    p1 <- plotExprs(sce)
    ggsave(paste0(plot_dir,"plotExprs_",preprocessing_type, type,
                  ifelse(subset_sce,"_subset",""),
                  ifelse(is.null(transform_fn),"","_log"),
                  ifelse(use_full_data,"_full",""),".png"), p1, width = 20,height = 20)
  },
  error = function(x) NA)

  tryCatch({
  p2 <- plotCounts(sce, group_by = "sample_id", color_by = "condition")
  ggsave(paste0(plot_dir,"plotCounts_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
                ifelse(use_full_data,"_full",""),".png"), p2)
  },
  error = function(x) NA)


  tryCatch({
  p3 <- pbMDS(sce, color_by = "condition", label_by = "sample_id")
  ggsave(paste0(plot_dir,"pbMDS_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
                ifelse(use_full_data,"_full",""),".png"), p3, width = 30,height = 30)
  },
  error = function(x) NA)


  p4 <- tryCatch({
  plotExprHeatmap(sce, scale = "last",
                  hm_pal = rev(hcl.colors(10, "YlGnBu")))
  },
  error = function(x) NA)
  png(paste0(plot_dir,"plotExprHeatmap_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
             ifelse(use_full_data,"_full",""),".png"),width = 1000)
  p4
  dev.off()



  tryCatch({
  p5 <- plotNRS(sce, features = "type", color_by = "condition")
  ggsave(paste0(plot_dir,"plotNRS_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
                ifelse(use_full_data,"_full",""),".pdf"), p5)
  },
  error = function(x) NA)
  print_partial_time()
}

print_total_time()
cat("run clustering\n")
t_tmp <- Sys.time()


sce <- CATALYST::cluster(sce, features = "type",
               xdim = 20, ydim = 20, maxK = 100, seed = 1234)
print_partial_time()
if (run_type == "test_run" | subset_sce){

  cat("plotting 3 plots\n")
  t_tmp <- Sys.time()

  p6 <- tryCatch({
  plotExprHeatmap(sce, features = "type",
                  by = "cluster_id", k = "meta50",
                  bars = TRUE, perc = TRUE)
  },
  error = function(x) NA)
  png(paste0(plot_dir,"plotExprHeatmap_2_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
             ifelse(use_full_data,"_full",""),".png"),width = 1000)
  p6
  dev.off()

  p8 <- tryCatch({
  plotMultiHeatmap(sce,
                   hm1 = "type", hm2 = "CD27", k = "meta50",
                   row_anno = FALSE, bars = TRUE, perc = TRUE)
  },
  error = function(x) NA)
  png(paste0(plot_dir,"plotMultiHeatmap_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
             ifelse(use_full_data,"_full",""),".png"),width = 1500,height = 1000)
  p8
  dev.off()

  print_partial_time()
}
cat("run umap\n")
t_tmp <- Sys.time()
sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")
print_partial_time()

if (data_type[[iii]] == "Spike_in"){
  colData(sce) <- cbind(colData(sce),metadata_frames)
}

cat("save sce as rds\n")
saveRDS(sce,paste0(result_dir,"/sce_", preprocessing_type, type,ifelse(subset_sce,"_subset",""), ifelse(is.null(transform_fn),"","_log"),
                   ifelse(use_full_data,"_full",""),".rds"))
print_total_time()
if (run_type == "test_run" | subset_sce){

  cat("plotting 4 plots\n")
  t_tmp <- Sys.time()

  tryCatch({
  p9 <- plotDR(sce, "UMAP", color_by = "meta50",facet_by = "condition")
  ggsave(paste0(plot_dir,"plotDR_UMAP_",preprocessing_type,ifelse(subset_sce,"_subset",""), type,ifelse(is.null(transform_fn),"","_log"),
                ifelse(use_full_data,"_full",""),".png"), p9)
  },
  error = function(x) NA)


  tryCatch({
  p10 <- plotCodes(sce, k = "meta50")
  ggsave(paste0(plot_dir,"plotCodes_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
                ifelse(use_full_data,"_full",""),".png"), p10)
  },
  error = function(x) NA)

  p11 <- tryCatch({
  plotMultiHeatmap(sce,
                   hm1 = "type", hm2 = "CD27", k = "meta50", m = "meta20",
                   row_anno = FALSE, col_anno = FALSE, bars = TRUE, perc = TRUE)
  },
  error = function(x) NA)
  png(paste0(plot_dir,"plotMultiHeatmap_2_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
             ifelse(use_full_data,"_full",""),".png"),width = 1000)
  p11
  dev.off()




  tryCatch({
  p12 <- plotAbundances(sce, k = "meta20", by = "sample_id")
  ggsave(paste0(plot_dir,"plotAbundances_",preprocessing_type, type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
                ifelse(use_full_data,"_full",""),".png"), p12)
  },
  error = function(x) NA)
  print_partial_time()
}

ei <- metadata(sce)$experiment_info
ei$survival_time <- as.numeric(as.character(ei$survival_time))

if (!second_cov){
  contrast <- createContrast(c(0, 1))
  da_formula1 <- createFormula(ei,
                               cols_fixed = c("survival_time"),
                               cols_random = c("patient_id"),
                               event_indicator = "Status")
} else{
  contrast <- createContrast(c(0, 1, 0))
  da_formula1 <- createFormula(ei,
                               cols_fixed = c("survival_time","condition"),
                               cols_random = c("sample_id","patient_id"),
                               event_indicator = "Status")
}

cs_by_s <- split(seq_len(ncol(sce)), colData(sce)$sample_id)
experiment_info <- metadata(sce)$experiment_info
cs <- unlist(cs_by_s[as.character(experiment_info$sample_id)])
es <- t(assays(sce)[["exprs"]])[cs, , drop = FALSE]

# create SummarizedExperiment (in transposed format compared to SingleCellExperiment)
d_se <- SummarizedExperiment(
  assays = list(exprs = es),
  rowData = colData(sce)[cs, ],
  colData = rowData(sce),
  metadata =  metadata(sce)
)
d_counts <- calcCounts(d_se)
metadata(d_counts) <- metadata(sce)
saveRDS(list(d_counts,da_formula1),paste0(result_dir,"/assay_sce_", preprocessing_type, type,ifelse(subset_sce,"_subset",""),
                                          ifelse(use_full_data,"_full",""), ".rds"))
print_total_time()

for (clu in c("som400","meta100","meta50","meta20","meta10")){
  for (imp_meth in c("cc","rs","km","mrl")) {
    cat("start DA for", imp_meth,"\n")
    t_tmp <- Sys.time()
    da_res1 <- tryCatch(diffcyt(sce,
                                formula = da_formula1, contrast = contrast,
                                analysis_type = "DA", method_DA = "diffcyt-DA-censored-GLMM",
                                clustering_to_use = clu, verbose = TRUE, mi_reps = 200,
                                imputation_method = imp_meth, BPPARAM = BiocParallel::MulticoreParam(workers=5)),
                        error = function(e) e)
    print_partial_time()
    saveRDS(da_res1,paste0(result_dir,"/da_res1_",imp_meth,"_", preprocessing_type,
                           type,ifelse(subset_sce,"_subset",""),ifelse(is.null(transform_fn),"","_log"),
                           ifelse(use_full_data,"_full",""),"_",clu, ".rds"))
  }
}



cat("script finished in \n")
print_total_time()
