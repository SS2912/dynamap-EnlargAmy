# dynamap-EnlargAmy
Functions and pipeline for the Amygdala enlargment project (SS, JT, 2024)

**Main pipeline script**
 - Amygdala_Enlarg_pipelineFC.m  --> calls all other functions below

**Main functions**
  - FC_amy_enlarg.m     --> Calculates node strength for each channel and subj in dir_data, outputs a table for analyses
  - graph_analysis.m    --> Calculates different garph measures for each channel/subjects and gives a result table for each channel

**Support functions**
  - reduce_h2matrix.m   --> for reduce to only some ROIs
  - FCmatrix_no_rep.m   --> same as in toolbox (1 FC value mean/median per each anat lable VEP)
  - ins_findmaxh2_sara.m
  - ins_countlinks_sara.m
  - apply_topology_measures_on_matrix.m (SMV)   --> used in graph_analysis
  - fdr.m
