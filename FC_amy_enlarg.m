function [FCtable, h2_amyhpc] = FC_Joy(dir_data, dir_info, varargin)
%
% Calculates node strength for each channel and subj in dir_data
% Syntax:  
%    FCtable = FC_PTSD(dir_data, dir_info, varargin)
%
% Inputs:
%   dir_data    - input adress of folder containing h2 results from bids pipeline with sub-code at beginning 1 file per subject with subj code, all subjects in the same folder 
%   dir_info    - dir of info table with info on each channel
% 
% Varargin:
%   thr_h2, thr_lag - optional thershold for h2 (thr_h2=0 for strength, >0 for degrees) and lag. Default: thr_h2 = 0; thr_lag = 0
%   norm            - 1 if you want to normalise each subj by number of channels (Default= 1)
%   graph           - 1 (default) to show connectivity matrices of zvalues compared to baseline chosen in base, 0 to not show any graph
%   roi             - "all" (default, h2 between all channels) or "EZ", "EZPZ", "NI" to calculate node strength only in subset of EZ channels (or non-inv channels)
%   meth            - how to obtain signle value of h2 for each region ("median", "mean". default = "median")
%
% Output:
%   FCtable     - 1 row per channel (and per subject), mean node stength OUT and TOT
%
% Required functions: 
%   - ins_countlinks.m,
%   - FCmatrix_no_rep.m (TOOLBOX_Sara, Onedrive Phd projects)
%   - optional: reduce_h2matrix.m (if intraconnectivity only in "roi" regions)
%
% Authors: Sara Simula (original: Feb 2024. Last version: March 24)

% 1. Optional variables: default values
thr_h2 = 0;  % thr_h2=0 for strength, >0 for degrees. Default: strength
thr_lag = 0; % value in ms, optional input. Default: 0 % TO CHECK???
norm =1; %% normalise strength by number of channels
graph = 0;
roi = "all";
category = "ipsi/contra";
meth = "median";

for ii = 1:2:nargin-2
        if strcmp('thr_h2', varargin{ii})
            thr_h2 = varargin{ii+1}; 
        elseif strcmp('thr_lag', varargin{ii})
            thr_lag = varargin{ii+1};
        elseif strcmp('norm', varargin{ii})
            norm = varargin{ii+1};
        elseif strcmp('graph', varargin{ii})
            graph = varargin{ii+1};
        elseif strcmp('roi', varargin{ii})
            roi = varargin{ii+1};
        elseif strcmp('category', varargin{ii})
            category = varargin{ii+1};
        elseif strcmp('meth', varargin{ii})
            meth = varargin{ii+1};
        elseif strcmp('band', varargin{ii})
            band = varargin{ii+1};
        end
end

%% 2. Set directory for h2 files to analyse
cd(dir_data)
myfiles = dir('*.mat');

% create new folder for h2 graphs matrices
newfold = strcat("matrices_figs_", meth);
mkdir(newfold);
%% 3. Calculate the strength/degrees for each channel and each patient and each channel
subj_info= [];
widedata=[];
FCtable = [];
h2_amyhpc = struct();

chan_infoAll = readtable(dir_info);
varnames = chan_infoAll.Properties.VariableNames;

for index=1:length(myfiles)
%     index=8; %debug
    rest = load(myfiles(index).name);
    subj = string(extractBetween(myfiles(index).name, 'sub-', '_ses'));
    channels = upper(string(rest.electrode_names));
   
    % read and add info on subject and channels
    chan_infoSub = chan_infoAll(chan_infoAll.subject == subj, :);

    b = table2cell(chan_infoSub);
    subj_info = string(b);

    amy_enlarg = subj_info(1, strcmp(varnames,"AE_1_5"));

    %% 3.1 Only take 1 value per each structure (varragin for max, median or mean)
    anat_column = contains(varnames,"brain_area");
    channel_col = contains(varnames,"channel");

    clearvars new_size new_row new_col A nbrwin

    subj_anat = subj_info(:, anat_column);
    subj_chan = subj_info(:, channel_col);
   
    [h2_mean, h2_median] = FCmatrix_no_rep(rest.aw_h2, rest.aw_lag, subj_anat, subj_chan, channels);
   
    % substitute new matrices and chan names instead of old rest.aw_h2 etc
    rest.aw_h2 = [];
    rest.aw_lag = [];
    rest.electrode_names = {};

    switch meth
        case "mean"
            rest.aw_h2  = h2_mean.h2;
            rest.aw_lag = h2_mean.lag;
        case "median"
            rest.aw_h2  = h2_median.h2;
            rest.aw_lag = h2_median.lag;
    end

    rest.electrode_names = h2_mean.chan; %mean or median channels are the same
    channels = upper(string(rest.electrode_names));

    %% 3.2 Reduce matrix to only roi regions if wanted
    % to change with correct inputs 
%     if ~strcmp("all", roi)
%         [h2_roi, subj_info_new] = reduce_h2matrix(h2_raw, roi, subj_info, roi_col);
%         rest = [];
%         rest = h2_roi;
%         clearvars subj_info
%         subj_info = subj_info_new;
%     end

    
    %% 3.3 Calculate the in, out, tot strengths/degrees using the func countlinks in graphcompare
    if thr_h2==0 && norm==1
        
        [linksrest,~]=ins_countlinks(rest,thr_h2);     % linkspre contains 1 row per channel, 1 column per window of h2
        OUTrest = mean(linksrest.outstrength_norm,2);  % median across windows of calculation h2
        TOTrest = mean(linksrest.totstrength_norm,2);  % median across windows of calculation h2
  
    else %degrees YET TO BE NORMALIZED!!!!! NOT TO USE FOR NOW
        [linksrest,~]=ins_countlinks(rest,thr_h2);
        OUTrest = mean(linksrest.outdegree,2);
        TOTrest = mean(linksrest.totdegree,2);

    end

    tmp = [channels', OUTrest, TOTrest];
    
    result_col = length(subj_info)+1;
    tmp(:,1) = upper(tmp(:,1));

    for chan = 1:length(tmp(:,1))
        match = sum(contains(upper(subj_info(:,channel_col)), tmp(chan,1)));
        if match 
            subj_info(contains(upper(subj_info(:,channel_col)), tmp(chan,1)), result_col:result_col+1) = tmp(chan,2:3);
        else 
            subj_info(contains(upper(subj_info(:,channel_col)), tmp(chan,1)), result_col:result_col+1) = NaN;
        end
    end
   
    clearvars tmp match
    subj_info = subj_info(~ismissing(subj_info(:,result_col)),:);
        
    FCtable = [FCtable; subj_info];

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% 4. Graphs
   if graph

       for j=1:length(channels)
           match = strcmp(upper(subj_info(:,channel_col)), channels(j));
           if sum(match)
               anat(j) = subj_info(match, strcmp(varnames,"brain_area"));
           else
               anat(j) = "missing in excel";
           end
       end

       lim_max = max(max(mean(rest.aw_h2,3)));
       
       
       color = viridis;
       figure('Name', strcat(subj, "- amy enlargement: ", amy_enlarg))
    
   
       imagesc(mean(rest.aw_h2,3))
       colormap(color)
       title("FC matrix")
       colorbar
       xticks(1:length(channels));
       xticklabels(channels);
       xtickangle(60)
       yticks(1:length(channels));
       yticklabels(strcat(channels, "-", anat));
       ylabel = "channel";
       xlabel = "channel";
       clim([0, lim_max])
  
  
       print(fullfile(dir_data, newfold, subj), '-dpng');
       
    
   end   
    clearvars  subj channels check linksrest OUTrest TOTrest subj_info min max rest h2_array anat
    
end

%% save the table
% switch category
%    case "ipsi/contra"
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_ipsi_ipsi','TOTrest_ipsi_ipsi','OUTrest_contra_contra', 'TOTrest_contra_contra', 'OUTrest_ipsi_contra', 'TOTrest_ipsi_contra', 'OUTrest_contra_ipsi', 'TOTrest_contra_ipsi'};
%    case "inv/ni"
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_inv_inv','TOTrest_inv_inv','OUTrest_ni_ni', 'TOTrest_ni_ni', 'OUTrest_inv_ni', 'TOTrest_inv_ni', 'OUTrest_ni_inv', 'TOTrest_ni_inv'};
% end
    
varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all'};
FCtable = array2table(FCtable, 'VariableNames', varnames);


 % if also inv-ni or ipsi-contra
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_ipsi_ipsi','TOTrest_ipsi_ipsi','OUTrest_contra_contra', 'TOTrest_contra_contra', 'OUTrest_ipsi_contra', 'TOTrest_ipsi_contra', 'OUTrest_contra_ipsi', 'TOTrest_contra_ipsi'};
%        varnames = {varnames{:}, 'OUTrest_all', 'TOTrest_all', 'OUTrest_inv_inv','TOTrest_inv_inv','OUTrest_ni_ni', 'TOTrest_ni_ni', 'OUTrest_inv_ni', 'TOTrest_inv_ni', 'OUTrest_ni_inv', 'TOTrest_ni_inv'};


end   
