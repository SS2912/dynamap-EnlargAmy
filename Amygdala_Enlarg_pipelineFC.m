%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%     JOY's AMY pipeline     %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sara simula, February 2024
% pipeline of analysis for AMYGDALA ENLARGEMENT project with Joy T. (FC only)

%%%%%%%%%%%%%  1. FC analysis  %%%%%%%%%%%%%%%%%%%
clear
close all

% Varargin:
%   thr_h2, thr_lag - optional thershold for h2 (thr_h2=0 for strength, >0 for degrees) and lag. Default: thr_h2 = 0; thr_lag = 0
%   norm            - 1 if you want to nromalise each subj by number of channels (Default= 1)
%   graph           - 1 (default) to show connectivity matrices of zvalues compared to baseline chosen in base, 0 to not show any graph
%   roi             - "all" (default, h2 between all channels) or "EZ", "EZPZ", "NI" to calculate node strength only in subset of EZ channels (or non-inv channels)
%   meth            - how to obtain signle value of h2 for each region ("median", "mean". default = "median")

% dir_info = "\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\big_table_draft_V2.xlsx";
dir_info = "\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\big_table_draft_V5_final.xlsx";
date = datestr(clock,'YYYY-mm-dd_HH.MM');

% category = "ipsi/contra"; % used for PTSD but not used here, can be adapted later

%% Broad %%
close all
method = 'mean';
dir_data = "\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\1.rawdata\broadband";
[FCtable, ~] = FC_amy_enlarg(dir_data, dir_info, 'meth', method, 'graph', 0);

name = strcat('\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\2.analysis\FCtables\', 'FCtable_broad', date, '-', method, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name)  

%% Other bands %%
close all
method = 'mean';
band = "theta";
dir_data = strcat("\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\1.rawdata\", band);
[FCtable, ~] = FC_(dir_data, dir_info, 'meth', method, 'graph', 1);

name = strcat('\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\2.analysis\FCtables\', 'FCtable_', band, date, '-', method, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(FCtable, name) 

%% to do for other bands

%%%%%%%%%%%%%  2. graph measures DRAFT  %%%%%%%%%%%%%%%%%
clear
cfg = struct(); % or cfg.deg_thr = 0.2; % or other number
symmetric = 1;  % symmetrize (1) or not (0) the connectivity matrix to make it Undirected (WU,BU)
graphs = 0;

dir_info = "\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\big_table_draft_V3_mergedROI.xlsx";
dir_data = "\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\1.rawdata\broadband";

[topo_table_broad, topo_broad] = graph_analysis(dir_data, 'dir_info', dir_info);

name = strcat('\\dynaserv\Shared\Maria\Joy_project\Sara_analyses\2.analysis\FCtables\', 'graph_measures_broad', date, '.xlsx'); % when saving to csv, it doesnt save variable names
writetable(topo_table_broad, name)  