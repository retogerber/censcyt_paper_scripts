# Scripts for censored diffcyt (censcyt) manuscript

This repository contains all the scripts to recreate the analysis and plots in the **_censcyt_: censored covariates in differential abundance analysis in cytometry** manuscript.


It contains five subdirectories:
1. [Basic simulations](./Basic_simulations) 

   simulations to recreate Figures 1 and 2 plus Supplementary Figures S1-S4.

2. [Simulations based on real data](./Simulations_based_on_real_data) 

   simulations to recreate Figure 4 and 5 and Supplementary Figure S5.

3. [Case study](./Case_study) 

   scripts to recreate Figure 6 and Supplementary Figure S6.

4. [FlowCAP\_IV](./FlowCAP_IV) 

   contains some metadata of the FlowCAP 4 dataset, the raw data can be downloaded from the [flowrepository FR-FCM-ZZ99](http://flowrepository.org/id/FR-FCM-ZZ99).

5. [simulationStudyCensoredDiffcyt](./simulationStudyCensoredDiffcyt) 

   additional scripts for running the simulations.




## Rerun analysis

All simulations scripts were run using [Singularity](https://sylabs.io) containers. More specifically: the container can be pulled using 

```
singularity pull singularity_diffcyt.sif shub://retogerber/singularity_diffcyt    
```

the Rscripts can then be run with:

```
singularity exec --bind .:/home/retger/censored_diffcyt singularity_diffcyt.sif Rscript script.R

```

The bind points (after `--bind`) maybe need to be adjusted and swap `script.R` with the script you would like to run.



## Run order
All scripts rely on weibull parameters obtained from the FlowCAP 4 dataset ([weibull\_fits\_FlowCAP.rds](./Simulations_based_on_real_data)) which can be generated using the script [survival\_time\_simulation.R](survival_time_simulation.R)

The scripts in [Simulations based on real data](./Simulations_based_on_real_data) additionally depend on [da\_res1\_cc\_complete.rds](./Simulations_based_on_real_data) that can be obtained by running the Case study script [censoredGLMM\_complete.R](./Case_study/censoredGLMM_complete.R).


## Session Info for running the simulations
```
> sessionInfo()
R version 4.0.0 (2020-04-24)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-openmp/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] scales_1.1.1                flowCore_2.1.0             
 [3] SummarizedExperiment_1.19.5 DelayedArray_0.15.3        
 [5] matrixStats_0.56.0          Matrix_1.2-18              
 [7] Biobase_2.49.0              GenomicRanges_1.41.5       
 [9] GenomeInfoDb_1.25.1         IRanges_2.23.9             
[11] S4Vectors_0.27.12           BiocGenerics_0.35.4        
[13] magrittr_1.5                diffcyt_1.9.6              
[15] forcats_0.5.0               stringr_1.4.0              
[17] dplyr_1.0.1                 purrr_0.3.4                
[19] readr_1.3.1                 tidyr_1.1.1                
[21] tibble_3.0.3                ggplot2_3.3.2              
[23] tidyverse_1.3.0            

loaded via a namespace (and not attached):
  [1] TH.data_1.0-10              minqa_1.2.4                
  [3] colorspace_1.4-1            rjson_0.2.20               
  [5] ellipsis_0.3.1              circlize_0.4.10            
  [7] cytolib_2.1.6               XVector_0.29.2             
  [9] GlobalOptions_0.1.2         base64enc_0.1-3            
 [11] fs_1.4.1                    clue_0.3-57                
 [13] rstudioapi_0.11             hexbin_1.28.1              
 [15] CytoML_2.1.0                mvtnorm_1.1-1              
 [17] fansi_0.4.1                 lubridate_1.7.9            
 [19] xml2_1.3.2                  codetools_0.2-16           
 [21] splines_4.0.0               jsonlite_1.7.0             
 [23] nloptr_1.2.2.2              broom_0.5.6                
 [25] cluster_2.1.0               dbplyr_1.4.4               
 [27] png_0.1-7                   graph_1.67.1               
 [29] compiler_4.0.0              httr_1.4.1                 
 [31] backports_1.1.8             assertthat_0.2.1           
 [33] limma_3.45.6                cli_2.0.2                  
 [35] tools_4.0.0                 ncdfFlow_2.35.1            
 [37] igraph_1.2.5                gtable_0.3.0               
 [39] glue_1.4.1                  GenomeInfoDbData_1.2.3     
 [41] flowWorkspace_4.1.3         reshape2_1.4.4             
 [43] ggcyto_1.17.0               Rcpp_1.0.5                 
 [45] cellranger_1.1.0            vctrs_0.3.2                
 [47] nlme_3.1-148                lme4_1.1-23                
 [49] rvest_0.3.5                 lifecycle_0.2.0            
 [51] statmod_1.4.34              XML_3.99-0.5               
 [53] edgeR_3.31.4                zlibbioc_1.35.0            
 [55] MASS_7.3-51.6               zoo_1.8-8                  
 [57] RProtoBufLib_2.1.0          hms_0.5.3                  
 [59] RBGL_1.65.0                 sandwich_2.5-1             
 [61] RColorBrewer_1.1-2          ComplexHeatmap_2.5.3       
 [63] yaml_2.2.1                  gridExtra_2.3              
 [65] latticeExtra_0.6-29         stringi_1.4.6              
 [67] boot_1.3-25                 shape_1.4.4                
 [69] rlang_0.4.7                 pkgconfig_2.0.3            
 [71] bitops_1.0-6                lattice_0.20-41            
 [73] tidyselect_1.1.0            plyr_1.8.6                 
 [75] R6_2.4.1                    generics_0.0.2             
 [77] multcomp_1.4-13             DBI_1.1.0                  
 [79] pillar_1.4.6                haven_2.3.1                
 [81] withr_2.2.0                 survival_3.1-12            
 [83] RCurl_1.98-1.2              FlowSOM_1.21.0             
 [85] tsne_0.1-3                  modelr_0.1.8               
 [87] crayon_1.3.4                jpeg_0.1-8.1               
 [89] GetoptLong_1.0.2            locfit_1.5-9.4             
 [91] grid_4.0.0                  readxl_1.3.1               
 [93] data.table_1.13.0           blob_1.2.1                 
 [95] Rgraphviz_2.33.0            ConsensusClusterPlus_1.53.0
 [97] reprex_0.3.0                digest_0.6.25              
 [99] RcppParallel_5.0.2          munsell_0.5.0   
```

