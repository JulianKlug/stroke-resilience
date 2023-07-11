clear()
clc()

data_dir = '';
n_rand = 100;
flag = 1; %1 for clustering coefficient using Cws, 2 for clustering coefficient using transitivity

%iterate through all subdirectories (subjects)
subdirectories = dir(data_dir);
hw = waitbar(0,'Computing small worldness...');
for i = 1:length(subdirectories)
    if subdirectories(i).isdir && ~strcmp(subdirectories(i).name,'.') && ~strcmp(subdirectories(i).name,'..')
        sub_dir = [data_dir '/' subdirectories(i).name '/'];
        %find undirected_thresholded_graphs folder
        subject_subdirectories = dir(sub_dir);
        for j = 1:length(subject_subdirectories)
            if subject_subdirectories(j).isdir && strcmp(subject_subdirectories(j).name,'undirected_thresholded_graphs')
                graphs_dir = [sub_dir '/' subject_subdirectories(j).name '/'];
                %find all .mat files
                graphs_files = dir([graphs_dir '*.mat']);
                for k = 1:length(graphs_files)
                    %load .mat file
                    graph = load([graphs_dir graphs_files(k).name]);
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
                        waitbar(l/length(threshold_keys)/length(graphs_files)/length(subdirectories),hw);
                    end
                    %save graph_sigmas for specific graph file
                    save([graphs_dir graphs_files(k).name(1:end-4) '_sw_sigmas.mat'],'graph_sigmas');
                end
            end
        end

    end
end
