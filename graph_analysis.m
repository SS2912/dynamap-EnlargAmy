function [topo_table, topology] = graph_analysis(dir_data, varargin)
%
% Inputs:
%   dir_data    - input adress of folder containing h2 results from bids pipeline with sub-code at beginning 1 file per subject with subj code, all subjects in the same folder 
%   dir_info    - dir of info table with info on each channel
% 
% Varargin:
%   cfg         - optional parameters (fieldtrip config style) for analysis (threshold etc)
%   symmetric   - symmetrize (1, default) or not (0) the connectivity matrix to make it Undirected (WU,BU)
%   graphs      - if 1, plots the matrix after symmetrization 
%
%
% Output:
%   topo_table  - table with mean of some graph measures combined with subj_info table if present (saved as excel)
%   topology    - cell containing raw various graph measures calculated on the input h2 matrix, not averaged across windows 
%
% Required functions: 
%   - apply_topology_measures_on_matrix.m (SMV)
%
% Authors: Sara Simula (original: March 2024. Last version: )
% baseline

% 0. initialize default values and variables
cfg = struct();
symmetric = 1; 
graphs = 0;
add_info = 0;
cols = 3;
topo_table = [];
subj_info= [];

for ii = 1:2:nargin-1
        if strcmp('dir_info', varargin{ii})
            dir_info = varargin{ii+1}; 
            add_info = 1;
            cols = cols+1;
%         elseif strcmp('thr_lag', varargin{ii})
%             thr_lag = varargin{ii+1};
        end
end

%% 1. set directory for h2 matrices
cd(dir_data)
myfiles = dir('*.mat');
topology = cell(length(myfiles),cols);

% add info if necessary
if add_info 
    chan_infoAll = readtable(dir_info);
    varnames = chan_infoAll.Properties.VariableNames;
end

%% 2. Loop on subject and calculate graph measures

% % debug:
% keep = myfiles(1);
% clearvars myfiles
% myfiles = keep;

for sub_idx=1:length(myfiles)
%     sub_idx = 31;
    clearvars h2 rest
    h2 = load(myfiles(sub_idx).name);
    subj = string(extractBetween(myfiles(sub_idx).name, 'sub-', '_ses'));
    channels = upper(string(h2.electrode_names));
    conn_mat = h2.aw_h2;
    topology{sub_idx,1} = subj;
    topology{sub_idx,2} = channels;
    topology{sub_idx,3} = apply_topology_measures_on_matrix(conn_mat, cfg, symmetric, graphs);

    if add_info
        chan_infoSub = chan_infoAll(chan_infoAll.subject == subj, :);
        b = table2cell(chan_infoSub);
        subj_info = string(b);
        topology{sub_idx, 4} = subj_info;
        % calculate mean for 2 dimensional topological measures and put in info table
        all_measures        = topology{sub_idx,3};
        fnames              = fieldnames(all_measures);
        idx_mean            = 0;
        for i =1:length(fnames)
            if length(size(all_measures.(fnames{i}))) ==2
                if size(all_measures.(fnames{i}),1) == length(channels)
                    idx_mean = idx_mean +1;
                    result_col = length(subj_info)+1;
                    var_topo_mean{idx_mean} = fnames{i};
                    subj_info(:, result_col) = mean(all_measures.(fnames{i}), 2);
                elseif size(all_measures.(fnames{i}),1) == 1
                    idx_mean = idx_mean +1;
                    result_col = length(subj_info)+1;
                    var_topo_mean{idx_mean} = fnames{i};
                    subj_info(:, result_col) = repelem(mean(all_measures.(fnames{i}), 2), length(channels))';
                end
            end
        end
       
    end
    
topo_table = [topo_table; subj_info];

clearvars fnames subj_info b subj channels conn_mat all_measures result_col 
end

varnames = {varnames{:}, var_topo_mean{:}};
topo_table = array2table(topo_table, 'VariableNames', varnames);

end
