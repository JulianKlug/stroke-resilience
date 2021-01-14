%% script to compute small-world-ness and do statistical testing on all networks

clear all; close all;

% data paths
data_path = '/home/klug/working_data/stroke_resilience/randmio_connected_bin_pos_ST_perm01-1000_05-28-2020 23-28.mat';
output_path = '/home/klug/output/stroke_resilience';

% SmallWorldNess toolbox by Mark Humphries needs to be declared here
path_to_SmallWorldNess_toolbox = '/home/klug/utils/SmallWorldNess';
data_field_name = 'CM_thresh_HC_bin';

% timepoints = {'ST01', 'ST02', 'ST03'};
% for tp_idx = 1:length(timepoints)
% tp = string(timepoints(tp_idx));

secondary_sub_field_name = '';
group_name = 'HC'; % one of HC / TP1-3


% analysis parameters
compute_significance = false; % significance is calculated via montecarlo simulation for all data points (takes a lot of time)
Num_ER_repeats = 100;  % to estimate C and L numerically for E-R random graph
Num_S_repeats = 1000; % to get P-value for S; min P = 0.001 for 1000 samples
I = 0.95; % confidence interval

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

        small_worldness_results = SmallWorldNess_Method_comparison(A, Num_ER_repeats, Num_S_repeats, I, compute_significance);

        small_worldness_results.threshold = repmat(threshold, height(small_worldness_results), 1);
        
        subject_results = [subject_results; small_worldness_results];
    end 
    
    subject_results.subject = repmat(subject_idx, height(subject_results), 1);
    all_results = [ all_results; subject_results ];
end

writetable(all_results, fullfile(output_path, 'small_worldness_comparison.csv'))

