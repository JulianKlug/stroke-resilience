clear all
data_path = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/attacks240/attack_HC/HC_connectivity_matrix.mat';
output_path = '/Users/jk1/temp/stroke_resilience/output/repeated_random_attack';
addpath('/Users/jk1/matlab/toolboxes/SmallWorldNess_toolbox')
addpath('/Users/jk1/matlab/toolboxes/2019_03_03_BCT')

%set param
group = 'HC';
node_nb = 240;
%totalnod2attack = node_nb - 1;
attack = 'rnd_attack'; % 'trg_attack' , 'rnd_attack'
propthr=0:0.1:1;
max_random_deletion_iterations = 100;

%take binarized tresholded matrices 
load(data_path)

if strcmp(group, 'HC')
    CM_thresh_bin = CM_thresh_HC_bin;
end
    
nsubj = size(CM_thresh_bin.top100,3);

tic
%ground truth glob eff
for s=1:nsubj
    for p=1:numel(propthr)
        propname=['top' num2str(propthr(p)*100)];
        myCM=CM_thresh_bin.(propname)(:,:,s);
        L = myCM;
        [R,D] = reachdist(L);%binary(directed/undirected) connection matrix
        [lambda,efficiency] = charpath(D);
        GlobEff_bin.(propname)(s)=efficiency;
    end
end


%% calculate basal BC
all_CM = CM_thresh_bin;
cd(output_path)

%% Do multiple iterations of random deletion
for random_deletion_iteration = 1:max_random_deletion_iterations
    iter = ['random_deletion_iteration' num2str(random_deletion_iteration)]
%     setting seed so that same nodes are deleted for all groups
    rng(random_deletion_iteration)
    nod2delorder = randperm(node_nb);
    temp_CM = all_CM;
    %% loop iteratively on nodes to delete
    for attacknb = 1:node_nb
        n_nodes_attacked = ['n_nodes_attacked' num2str(attacknb)];
        for p = 1:numel(propthr)
            propname = ['top' num2str(propthr(p)*100)];
            nod2del = nod2delorder(attacknb);
            for s = 1:nsubj
                myCM = temp_CM.(propname)(:,:,s);
                myCM(:,nod2del) = 0;
                myCM(nod2del,:) = 0;
                temp_CM.(propname)(:,:,s) = myCM;
            end
        end

        %recalculate efficiency
        for p = 2:numel(propthr)
            propname = ['top' num2str(propthr(p)*100)];
            for s = 1:nsubj
                mynewCM = temp_CM.(propname)(:,:,s);
                L = mynewCM;
                [R,D] = reachdist(L);
                [lambda,efficiency] = charpath(D); % global efficiency
                GlobEff_bin_new.(iter).(n_nodes_attacked)(s,p) = efficiency;
            end
        end

    end
end
fname = fullfile(output_path, [group '_rep' num2str(max_random_deletion_iterations) '_rng_attack_' datestr(now,'mm-dd-yyyy HH-MM')  '.mat'])
save(fname, GlobEff_bin_new)

time = toc









        
        
