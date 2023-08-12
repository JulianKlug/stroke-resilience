clear()
clc()

 data_dir = '/home/users/k/klug/data/resilience/RS';
 bct_path = '/home/users/k/klug/utils/2019_03_03_BCT';
 sw_toolbox_path = '/home/users/k/klug/utils/SmallWorldNess';
n_rand = 100;
flag = 1; %1 for clustering coefficient using Cws, 2 for clustering coefficient using transitivity
recompute = 0; %1 to recompute, 0 to skip recomputing

 addpath(genpath(bct_path));
 addpath(genpath(sw_toolbox_path));


%iterate through all subdirectories (subjects)
subdirectories = dir(data_dir);
hw = waitbar(0,'Computing small worldness...');
for i = 1:length(subdirectories)
    if subdirectories(i).isdir && ~strcmp(subdirectories(i).name,'.') && ~strcmp(subdirectories(i).name,'..')
	%  display progress percentage
        fprintf('Processing folder progress: %d/%d\n',i,length(subdirectories));
        sub_dir = [data_dir '/' subdirectories(i).name '/'];
        %find undirected_thresholded_graphs folder
        subject_subdirectories = dir(sub_dir);
        for j = 1:length(subject_subdirectories)
            if subject_subdirectories(j).isdir && strcmp(subject_subdirectories(j).name,'undirected_thresholded_graphs')
                graphs_dir = [sub_dir '/' subject_subdirectories(j).name '/'];
                %find all .mat files (that do not end in _sw_sigmas.mat)
                graphs_files = dir([graphs_dir '*.mat']);
                graphs_files = graphs_files(~contains({graphs_files.name},'_sw_sigmas.mat'));
                graphs_files_paths = fullfile(graphs_dir, {graphs_files.name});
                parfor k = 1:length(graphs_files_paths)
                    graph_file = graphs_files_paths{k};
                    % check if graph_sigmas already exists
                    if exist([graph_file(1:end-4) '_sw_sigmas.mat'],'file') == 2 && recompute == 0
                        fprintf('Already computed %s\n',[graph_file(1:end-4) '_sw_sigmas.mat']);
                        continue
                    end

                    %load .mat file
                    graph = load(graph_file);
                    %iterate through keys
                    threshold_keys = fieldnames(graph);
                    graph_sigmas = {};
                    for l = 1:length(threshold_keys)
                        %get threshold
                        threshold = threshold_keys{l};
                        %get graph
                        graph_at_t = graph.(threshold);
                        small_worldness_sigma = monte_carlo_small_worldness(graph_at_t,n_rand, flag);
                        graph_sigmas.(threshold) = small_worldness_sigma;
                    end
                    %save graph_sigmas for specific graph file
                    % save([graphs_dir graphs_files(k).name(1:end-4) '_sw_sigmas.mat'],'graph_sigmas');
                    save_graph_sigmas(graph_sigmas, [graph_file(1:end-4) '_sw_sigmas.mat']); 
                end
                waitbar(j/length(subdirectories),hw);
            end
        end

    end
end

function save_graph_sigmas(graph_sigmas, fullfilename)
    save(fullfilename, 'graph_sigmas');
end
