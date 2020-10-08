# Scripts for censored diffcyt (censcyt) manuscript

This repository contains all the scripts to recreate the analysis and plots in the **_censcyt_: censored covariates in differential abundance analysis in cytometry** manuscript.


It contains five subdirectories:
1. [Basic simulations](./'Basic simulations') 

   simulations to recreate Figures 1 and 2 plus Supplementary Figures S1-S4.

2. [Simulations based on real data](./Simulations\ based\ on\ real\ data) 

   simulations to recreate Figure 4 and 5 and Supplementary Figure S5.

3. [Case study](./Case\ study) 

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
All scripts rely on weibull parameters obtained from the FlowCAP 4 dataset ([weibull\_fits\_FlowCAP.rds](./Simulations\ based\ on\ real\ data)) which can be generated using the script [survival\_time\_simulation.R](survival_time_simulation.R)

The scripts in [Simulations based on real data](./Simulations\ based\ on\ real\ data) additionally depend on [da\_res1\_cc\_complete.rds](./Simulations\ based\ on\ real\ data) that can be obtained by running the Case study script [censoredGLMM\_complete.R](./Case\ study/censoredGLMM_complete.R).





