%% script to compute basic metrics on all networks

clear all; close all;

% data paths
data_path = '/Users/jk1/temp/stroke_resilience/BNA_240_flipped_N32_retroicor_SBB4_prop_bin/HC/randmio_connected_bin_pos_HC_perm01-1000_06-03-2020 16-17.mat';
output_path = '/Users/jk1/temp/stroke_resilience/output/stroke_resilience';

% SmallWorldNess toolbox by Mark Humphries needs to be declared here
path_to_SmallWorldNess_toolbox = '/Users/jk1/matlab/toolboxes/SmallWorldNess_toolbox';
data_field_name = 'CM_thresh_HC_bin';

secondary_sub_field_name = '';
group_name = 'HC'; % one of HC / ST01-03


%% preprare paths and directories

script_path = mfilename('fullpath');
script_folder = script_path(1 : end - size(mfilename, 2));
addpath(genpath(script_folder));
addpath(genpath(path_to_SmallWorldNess_toolbox));

if ~(exist(data_path))
    fprintf('Data file does not exist. Please enter a valid file.')
end


if ~exist(fullfile(output_path, group_name))
   mkdir(output_path, group_name);
else
   warning('Output directory already exists, old results will be overwritten.'); 
end
output_path = fullfile(output_path, group_name);

%% load the adjacency matrix for the network network
raw_data = load(data_path);  % loads struct of data 

thresh_data = raw_data.(data_field_name);
if ~isempty(secondary_sub_field_name)
    thresh_data = thresh_data.(secondary_sub_field_name);
end
thresholds = fieldnames(thresh_data);

adjacency_matrix_shape = size(thresh_data.(string(thresholds(1))));
n_subjects = adjacency_matrix_shape(3);

all_results = cell2table({});

parfor subject_idx=1:n_subjects
    fprintf('Subject %d of %d. \n',subject_idx, n_subjects);
    subject_results = cell2table({});
    for threshold_idx = 1:numel(thresholds)
        threshold = string(thresholds(threshold_idx));
        fprintf('Subject %d: %0.1f %% of thresholds processed. \n',subject_idx,(threshold_idx/numel(thresholds) * 100));
        
        if isempty(secondary_sub_field_name) % check for subgroups 
            A = raw_data.(data_field_name).(threshold)(:,:,subject_idx);
        else
            A = raw_data.(data_field_name).(secondary_sub_field_name).(threshold)(:,:,subject_idx);
        end

        [n_nodes, n_edges, mean_degree, degree_distribution] = basic_graph_properties(A);
        headers = {'subject', 'threshold', 'n_nodes', 'n_edges', 'mean_degree', 'degree_distribution'};
        metrics = {subject_idx, threshold, n_nodes, n_edges, mean_degree, degree_distribution};
       
        metrics = cell2table(metrics);
        metrics.Properties.VariableNames = headers;

        all_results = [all_results; metrics];
    end 
end

writetable(all_results, fullfile(output_path, 'basic_graph_metrics.csv'))

